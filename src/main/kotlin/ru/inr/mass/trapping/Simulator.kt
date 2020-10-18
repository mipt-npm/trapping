/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ru.inr.mass.trapping

import org.apache.commons.math3.distribution.ExponentialDistribution
import org.apache.commons.math3.geometry.euclidean.threed.Rotation
import org.apache.commons.math3.geometry.euclidean.threed.SphericalCoordinates
import org.apache.commons.math3.random.RandomGenerator
import org.apache.commons.math3.util.FastMath.*
import org.apache.commons.math3.util.Pair
import org.apache.commons.math3.util.Precision
import ru.inr.mass.trapping.Scatter.counter
import ru.inr.mass.trapping.Scatter.debug

/**
 * @author Alexander Nozik
 * @property magneticField Longitudal magnetic field distribution
 * @property gasDensity gas density in 1/m^3
 */
class Simulator(
        val eLow: Double,
        val thetaTransport: Double,
        val thetaPinch: Double,
        val gasDensity: Double,// m^-3
        val bSource: Double,
        val magneticField: ((Double) -> Double)?
) {

    enum class EndState {

        ACCEPTED, //трэппинговый электрон попал в аксептанс
        REJECTED, //трэппинговый электрон вылетел через заднюю пробку
        LOWENERGY, //потерял слишком много энергии
        PASS, //электрон никогда не запирался и прошел напрямую, нужно для нормировки
        NONE
    }


    /**
     * Perform scattering in the given position
     *
     * @param pos
     * @return
     */
    private fun scatter(pos: State): State {
        //Вычисляем сечения и нормируем их на полное сечение
        var sigmaIon = Scatter.sigmaion(pos.e)
        var sigmaEl = Scatter.sigmael(pos.e)
        var sigmaExc = Scatter.sigmaexc(pos.e)
        val sigmaTotal = sigmaEl + sigmaIon + sigmaExc
        sigmaIon /= sigmaTotal
        sigmaEl /= sigmaTotal
        sigmaExc /= sigmaTotal

        //проверяем нормировку
        assert(Precision.equals(sigmaEl + sigmaExc + sigmaIon, 1.0, 1e-2))

        val alpha = pos.generator.nextDouble()

        val delta: Pair<Double, Double> = when {
            alpha < sigmaEl -> Scatter.randomel(pos.e, pos.generator)// elastic
            alpha > sigmaEl + sigmaExc -> Scatter.randomion(pos.e, pos.generator)//ionization case
            else -> Scatter.randomexc(pos.e, pos.generator)//excitation case
        }

        //Обновляем значени угла и энергии независимо ни от чего
        pos.substractE(delta.first)
        //Изменение угла
        pos.addTheta(delta.second / 180 * Math.PI)

        return pos
    }

    /**
     * Calculate distance to the next scattering
     *
     * @param e
     * @return
     */
    private fun freePath(generator: RandomGenerator, e: Double): Double {
        //FIXME redundant cross-section calculation
        //All cross-sections are in m^2
        return ExponentialDistribution(generator, 1.0 / Scatter.sigmaTotal(e) / gasDensity).sample()
    }

    /**
     * Calculate propagated position before scattering
     *
     * @param deltaL
     * @return z shift and reflection counter
     */
    private fun propagate(pos: State, deltaL: Double): State {
        // if magnetic field not defined, consider it to be uniform and equal bSource
        if (magneticField == null) {
            val deltaZ = deltaL * cos(pos.theta) // direction already included in cos(theta)
            //            double z0 = pos.z;
            pos.addZ(deltaZ)
            pos.l += deltaL

            //            //if we are crossing source boundary, check for end condition
            //            while (abs(deltaZ + pos.z) > SOURCE_LENGTH / 2d && !pos.isFinished()) {
            //
            //                pos.checkEndState();
            //            }
            //
            //            // if track is finished apply boundary position
            //            if (pos.isFinished()) {
            //                // remembering old z to correctly calculate total l
            //                double oldz = pos.z;
            //                pos.z = pos.direction() * SOURCE_LENGTH / 2d;
            //                pos.l += (pos.z - oldz) / cos(pos.theta);
            //            } else {
            //                //else just add z
            //                pos.l += deltaL;
            //                pos.addZ(deltaZ);
            //            }

            return pos
        } else {
            var curL = 0.0
            val sin2 = sin(pos.theta) * sin(pos.theta)

            while (curL <= deltaL - 0.01 && !pos.isFinished) {
                //an l step
                val delta = min(deltaL - curL, DELTA_L)

                val b = field(pos.z)

                val root = 1 - sin2 * b / bSource
                //preliminary reflection
                if (root < 0) {
                    pos.flip()
                }
                pos.addZ(pos.direction() * delta * sqrt(abs(root)))
                //                // change direction in case of reflection. Loss of precision here?
                //                if (root < 0) {
                //                    // check if end state occurred. seem to never happen since it is reflection case
                //                    pos.checkEndState();
                //                    // finish if it does
                //                    if (pos.isFinished()) {
                //                        //TODO check if it happens
                //                        return pos;
                //                    } else {
                //                        //flip direction
                //                        pos.flip();
                //                        // move in reversed direction
                //                        pos.z += pos.direction() * delta * sqrt(-root);
                //                    }
                //
                //                } else {
                //                    // move forward
                //                    pos.z += pos.direction() * delta * sqrt(root);
                //                    //check if it is exit case
                //                    if (abs(pos.z) > SOURCE_LENGTH / 2d) {
                //                        // check if electron is out
                //                        pos.checkEndState();
                //                        // finish if it is
                //                        if (pos.isFinished()) {
                //                            return pos;
                //                        }
                //                        // PENDING no need to apply reflection, it is applied automatically when root < 0
                //                        pos.z = signum(pos.z) * SOURCE_LENGTH / 2d;
                //                        if (signum(pos.z) == pos.direction()) {
                //                            pos.flip();
                //                        }
                //                    }
                //                }


                curL += delta
                pos.l += delta
            }
            return pos
        }
    }

    /**
     * Magnetic field in the z point
     *
     * @param z
     * @return
     */
    private fun field(z: Double): Double {
        return magneticField?.invoke(z) ?: bSource
    }

    /**
     * Симулируем один пробег электрона от начального значения и до вылетания из
     * иточника или до того момента, как энергия становится меньше eLow.
     */
    fun simulate(id: Long, generator: RandomGenerator, initEnergy: Double, initTheta: Double, initZ: Double): SimulationResult {
        assert(initEnergy > 0)
        assert(initTheta > 0 && initTheta < Math.PI)
        assert(abs(initZ) <= SOURCE_LENGTH / 2.0)

        val pos = State(generator, initEnergy, initTheta, initZ)

        while (!pos.isFinished) {
            val dl = freePath(generator, pos.e) // path to next scattering
            // propagate to next scattering position
            propagate(pos, dl)

            if (!pos.isFinished) {
                // perform scatter
                scatter(pos)
                // increase collision number
                pos.colNum++
                if (pos.e < eLow) {
                    //Если энергия стала слишком маленькой
                    pos.setEndState(EndState.LOWENERGY)
                }
            }
        }

        return SimulationResult(id, pos.endState, pos.e, pos.theta, initTheta, pos.colNum, pos.l)
    }

    fun resetDebugCounters() {
        debug = true
        counter.resetAll()
    }

    fun printDebugCounters() {
        if (debug) {
            counter.print(System.out)
        } else {
            throw RuntimeException("Debug not initiated")
        }
    }

    data class SimulationResult(val id: Long, val state: EndState, val E: Double, val theta: Double, val initTheta: Double, val collisionNumber: Int, val l: Double)

    /**
     * Current electron position in simulation. Not thread safe!
     *
     * @property e Current energy
     * @property theta Current theta recalculated to the field in the center of the source
     * @property z current z. Zero is the center of the source
     */
    private inner class State(val generator: RandomGenerator, internal var e: Double, internal var theta: Double, internal var z: Double) {
        /**
         * Current total path
         */
        var l = 0.0

        /**
         * Number of scatterings
         */
        var colNum = 0

        var endState = EndState.NONE

        internal val isForward: Boolean
            get() = theta <= Math.PI / 2

        internal val isFinished: Boolean
            get() = this.endState != EndState.NONE

        /**
         * @param dE
         * @return resulting E
         */
        internal fun substractE(dE: Double): Double {
            this.e -= dE
            return e
        }

        internal fun direction(): Double {
            return (if (isForward) 1 else -1).toDouble()
        }

        internal fun setEndState(state: EndState) {
            this.endState = state
        }

        /**
         * add Z and calculate direction change
         *
         * @param dZ
         * @return
         */
        internal fun addZ(dZ: Double): Double {
            this.z += dZ
            while (abs(this.z) > SOURCE_LENGTH / 2.0 && !isFinished) {
                // reflecting from back wall
                if (z < 0) {
                    if (theta >= PI - thetaTransport) {
                        setEndState(EndState.REJECTED)
                    }
                    z = if (isFinished) {
                        -SOURCE_LENGTH / 2.0
                    } else {
                        // reflecting from rear pinch
                        -SOURCE_LENGTH - z
                    }
                } else {
                    if (theta < thetaPinch) {
                        if (colNum == 0) {
                            //counting pass electrons
                            setEndState(EndState.PASS)
                        } else {
                            setEndState(EndState.ACCEPTED)
                        }
                    }
                    z = if (isFinished) {
                        SOURCE_LENGTH / 2.0
                    } else {
                        // reflecting from forward transport magnet
                        SOURCE_LENGTH - z
                    }
                }
                if (!isFinished) {
                    flip()
                }
            }
            return z
        }

        //        /**
        //         * Check if this position is an end state and apply it if necessary. Does not check z position.
        //         *
        //         * @return
        //         */
        //        private void checkEndState() {
        //            //accepted by spectrometer
        //            if (theta < thetaPinch) {
        //                if (colNum == 0) {
        //                    //counting pass electrons
        //                    setEndState(EndState.PASS);
        //                } else {
        //                    setEndState(EndState.ACCEPTED);
        //                }
        //            }
        //
        //            //through the rear magnetic pinch
        //            if (theta >= PI - thetaTransport) {
        //                setEndState(EndState.REJECTED);
        //            }
        //        }

        /**
         * Reverse electron direction
         */
        fun flip() {
            if (theta < 0 || theta > PI) {
                throw Error()
            }
            theta = PI - theta
        }

        /**
         * Magnetic field in the current point
         *
         * @return
         */
        fun field(): Double {
            return this@Simulator.field(z)
        }

        /**
         * Real theta angle in current point
         *
         * @return
         */
        fun realTheta(): Double {
            if (magneticField == null) {
                return theta
            } else {
                var newTheta = asin(min(abs(sin(theta)) * sqrt(field() / bSource), 1.0))
                if (theta > PI / 2) {
                    newTheta = PI - newTheta
                }

                assert(!java.lang.Double.isNaN(newTheta))
                return newTheta
            }
        }


        /**
         * Сложение вектора с учетом случайно распределения по фи
         *
         * @param dTheta
         * @return resulting angle
         */
        fun addTheta(dTheta: Double): Double {
            //Генерируем случайный фи
            val phi = generator.nextDouble() * 2.0 * Math.PI

            //change to real angles
            val realTheta = realTheta()
            //Создаем начальный вектор в сферических координатах
            val init = SphericalCoordinates(1.0, 0.0, realTheta + dTheta)
            // Задаем вращение относительно оси, перпендикулярной исходному вектору
            val rotate = SphericalCoordinates(1.0, 0.0, realTheta)
            // поворачиваем исходный вектор на dTheta
            val rot = Rotation(rotate.cartesian, phi, null)

            val result = rot.applyTo(init.cartesian)

            val newTheta = acos(result.z)

            //            //следим чтобы угол был от 0 до Pi
            //            if (newtheta < 0) {
            //                newtheta = -newtheta;
            //            }
            //            if (newtheta > Math.PI) {
            //                newtheta = 2 * Math.PI - newtheta;
            //            }

            //change back to virtual angles
            if (magneticField == null) {
                theta = newTheta
            } else {
                theta = asin(sin(newTheta) * sqrt(bSource / field()))
                if (newTheta > PI / 2) {
                    theta = PI - theta
                }
            }

            require(!theta.isNaN())

            return theta
        }

    }

    companion object {

        const val SOURCE_LENGTH = 3.0
        private const val DELTA_L = 0.1 //step for propagate calculation
    }


}
