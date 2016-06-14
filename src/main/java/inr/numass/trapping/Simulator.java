/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package inr.numass.trapping;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.SphericalCoordinates;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.Pair;
import org.apache.commons.math3.util.Precision;

import static org.apache.commons.math3.util.FastMath.*;

/**
 * @author Darksnake
 */
public class Simulator {

    public static final double SOURCE_LENGTH = 3d;
    private static final double DELTA_L = 0.1; //step for propagate calculation

    private RandomGenerator generator;
    private Scatter scatter;
    double eLow = 14000d;//default value
    double thetaTransport = 24.107064 / 180 * Math.PI;// default value
    double thetaPinch = 19.481097 / 180 * Math.PI;// default value
    double gasDensity = 5e18;// m^-3
    double bSource = 0.6;
    private UnivariateFunction magneticField;


    public Simulator() {
        generator = new JDKRandomGenerator();// Be careful here since it is not synchronized
        scatter = new Scatter(generator);
    }

    public Simulator(double Bsource, double Btransport, double Bpinch, double Elow) {
        this();
        setFields(Bsource, Btransport, Bpinch);
        setELow(Elow);
    }

    public enum EndState {

        ACCEPTED,//трэппинговый электрон попал в аксептанс
        REJECTED,//трэппинговый электрон вылетел через заднюю пробку
        LOWENERGY,//потерял слишком много энергии
        PASS,//электрон никогда не запирался и прошел напрямую, нужно для нормировки
        NONE
    }

    public final void setELow(double eLow) {
        this.eLow = eLow;
    }

    public final void setFields(double Bsource, double Btransport, double Bpinch) {
        this.bSource = Bsource;
        this.thetaTransport = Math.asin(Math.sqrt(Bsource / Btransport));
        this.thetaPinch = Math.asin(Math.sqrt(Bsource / Bpinch));
    }

    /**
     * Set gas density in 1/m^3
     *
     * @param gasDensity
     */
    public void setGasDensity(double gasDensity) {
        this.gasDensity = gasDensity;
    }

    /**
     * Longitudal magnetic field distribution
     *
     * @param magneticField
     */
    public void setFieldFunc(UnivariateFunction magneticField) {
        this.magneticField = magneticField;
    }

    /**
     * Perform scattering in the given position
     *
     * @param pos
     * @return
     */
    private State scatter(State pos) {
        //Вычисляем сечения и нормируем их на полное сечение
        double sigmaIon = scatter.sigmaion(pos.e);
        double sigmaEl = scatter.sigmael(pos.e);
        double sigmaExc = scatter.sigmaexc(pos.e);
        double sigmaTotal = sigmaEl + sigmaIon + sigmaExc;
        sigmaIon /= sigmaTotal;
        sigmaEl /= sigmaTotal;
        sigmaExc /= sigmaTotal;

        //проверяем нормировку
        assert Precision.equals(sigmaEl + sigmaExc + sigmaIon, 1, 1e-2);

        double alpha = generator.nextDouble();

        Pair<Double, Double> delta;
        if (alpha > sigmaEl) {
            if (alpha > sigmaEl + sigmaExc) {
                //ionization case
                delta = scatter.randomion(pos.e);
            } else {
                //excitation case
                delta = scatter.randomexc(pos.e);
            }
        } else {
            // elastic
            delta = scatter.randomel(pos.e);
        }

        //Обновляем значени угла и энергии независимо ни от чего
        pos.substractE(delta.getFirst());
        //Изменение угла
        pos.addTheta(delta.getSecond() / 180 * Math.PI);

        return pos;
    }

    /**
     * Calculate distance to the next scattering
     *
     * @param e
     * @return
     */
    private double freePath(double e) {
        //FIXME redundant cross-section calculation
        //All cross-sections are in m^2
        return new ExponentialDistribution(generator, 1d / scatter.sigmaTotal(e) / gasDensity).sample();
    }

    /**
     * Calculate propagated position before scattering
     *
     * @param deltaL
     * @return z shift and reflection counter
     */
    private State propagate(State pos, double deltaL) {
        // if magnetic field not defined, consider it to be uniform and equal bSource
        if (magneticField == null) {
            double deltaZ = deltaL * cos(pos.theta); // direction already included in cos(theta)
            double z0 = pos.z;
            pos.addZ(deltaZ);
            pos.l += abs(pos.z - z0) / cos(pos.theta);

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

            return pos;
        } else {
            double curL = 0;
            double sin2 = sin(pos.theta) * sin(pos.theta);

            while (curL <= deltaL - 0.01 && !pos.isFinished()) {
                //an l step
                double delta = min(deltaL - curL, DELTA_L);

                double b = field(pos.z);

                double root = 1 - sin2 * b / bSource;
                //preliminary reflection
                if (root < 0) {
                    pos.flip();
                }
                pos.addZ(pos.direction() * delta * sqrt(abs(root)));
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


                curL += delta;
                pos.l += delta;
            }
            return pos;
        }
    }

    /**
     * Magnetic field in the z point
     *
     * @param z
     * @return
     */
    private double field(double z) {
        if (magneticField == null) {
            return bSource;
        } else {
            return magneticField.value(z);
        }
    }

    /**
     * Симулируем один пробег электрона от начального значения и до вылетания из
     * иточника или до того момента, как энергия становится меньше eLow.
     */
    public SimulationResult simulate(double initEnergy, double initTheta, double initZ) {
        assert initEnergy > 0;
        assert initTheta > 0 && initTheta < Math.PI;
        assert abs(initZ) <= SOURCE_LENGTH / 2d;

        State pos = new State(initEnergy, initTheta, initZ);

        while (!pos.isFinished()) {
            double dl = freePath(pos.e); // path to next scattering
            // propagate to next scattering position
            propagate(pos, dl);

            if (!pos.isFinished()) {
                // perform scatter
                scatter(pos);
                // increase collision number
                pos.colNum++;
                if (pos.e < eLow) {
                    //Если энергия стала слишком маленькой
                    pos.setEndState(EndState.LOWENERGY);
                }
            }
        }

        SimulationResult res = new SimulationResult(pos.endState, pos.e, pos.theta, initTheta, pos.colNum, pos.l);
        return res;
    }

    public void resetDebugCounters() {
        scatter.debug = true;
        scatter.counter.resetAll();
    }

    public void printDebugCounters() {
        if (scatter.debug) {
            scatter.counter.print(System.out);
        } else {
            throw new RuntimeException("Debug not initiated");
        }
    }

    public static class SimulationResult {

        public SimulationResult(EndState state, double E, double theta, double initTheta, int collisionNumber, double l) {
            this.state = state;
            this.E = E;
            this.theta = theta;
            this.initTheta = initTheta;
            this.collisionNumber = collisionNumber;
            this.l = l;
        }

        public EndState state;
        public double E;
        public double initTheta;
        public double theta;
        public int collisionNumber;
        public double l;
    }

    /**
     * Current electron position in simulation. Not thread safe!
     */
    private class State {
        public State(double e, double theta, double z) {
            this.e = e;
            this.theta = theta;
            this.z = z;
        }

        /**
         * Current energy
         */
        double e;
        /**
         * Current theta recalculated to the field in the center of the source
         */
        double theta;
        /**
         * Current total path
         */
        double l = 0;

        /**
         * Number of scatterings
         */
        int colNum = 0;

        /**
         * current z. Zero is the center of the source
         */
        double z;

        EndState endState = EndState.NONE;

        /**
         * @param dE
         * @return resulting E
         */
        double substractE(double dE) {
            this.e -= dE;
            return e;
        }

        boolean isForward() {
            return theta <= Math.PI / 2;
        }

        double direction() {
            return isForward() ? 1 : -1;
        }

        boolean isFinished() {
            return this.endState != EndState.NONE;
        }

        void setEndState(EndState state) {
            this.endState = state;
        }

        /**
         * add Z and calculate direction change
         *
         * @param dZ
         * @return
         */
        double addZ(double dZ) {
            this.z += dZ;
            while (abs(this.z) > SOURCE_LENGTH / 2d && !isFinished()) {
                checkEndState();
                if (!isFinished()) {
                    flip();
                }
                // reflecting from back wall
                if (z < 0) {
                    if (isFinished()) {
                        z = -SOURCE_LENGTH / 2d;
                    } else {
                        // reflecting from rear pinch
                        z = -SOURCE_LENGTH - z;
                    }
                } else {
                    if (isFinished()) {
                        z = SOURCE_LENGTH / 2d;
                    } else {
                        // reflecting from forward transport magnet
                        z = SOURCE_LENGTH - z;
                    }
                }
            }
            return z;
        }

        /**
         * Check if this position is an end state and apply it if necessary. Does not check z position.
         *
         * @return
         */
        private void checkEndState() {
            //accepted by spectrometer
            if (theta < thetaPinch) {
                if (colNum == 0) {
                    //counting pass electrons
                    setEndState(EndState.PASS);
                } else {
                    setEndState(EndState.ACCEPTED);
                }
            }

            //through the rear magnetic pinch
            if (theta >= PI - thetaTransport) {
                setEndState(EndState.REJECTED);
            }
        }

        /**
         * Reverse electron direction
         */
        void flip() {
            theta = PI - theta;
        }

        /**
         * Magnetic field in the current point
         *
         * @return
         */
        double field() {
            return Simulator.this.field(z);
        }

        /**
         * Real theta angle in current point
         *
         * @return
         */
        double realTheta() {
            if (magneticField == null) {
                return theta;
            } else {
                double newTheta = asin(min(abs(sin(theta)) * sqrt(field() / bSource), 1));
                if (theta > PI / 2) {
                    newTheta = PI - newTheta;
                }

                assert !Double.isNaN(newTheta);
                return newTheta;
            }
        }


        /**
         * Сложение вектора с учетом случайно распределения по фи
         *
         * @param dTheta
         * @return resulting angle
         */
        double addTheta(double dTheta) {
            //Генерируем случайный фи
            double phi = generator.nextDouble() * 2 * Math.PI;

            //change to real angles
            double realTheta = realTheta();
            //Создаем начальный вектор в сферических координатах
            SphericalCoordinates init = new SphericalCoordinates(1, 0, realTheta + dTheta);
            // Задаем вращение относительно оси, перпендикулярной исходному вектору
            SphericalCoordinates rotate = new SphericalCoordinates(1, 0, realTheta);
            // поворачиваем исходный вектор на dTheta
            Rotation rot = new Rotation(rotate.getCartesian(), phi, null);

            Vector3D result = rot.applyTo(init.getCartesian());

            double newtheta = acos(result.getZ());

//            //следим чтобы угол был от 0 до Pi
//            if (newtheta < 0) {
//                newtheta = -newtheta;
//            }
//            if (newtheta > Math.PI) {
//                newtheta = 2 * Math.PI - newtheta;
//            }

            //change back to virtual angles
            if (magneticField == null) {
                theta = newtheta;
            } else {
                theta = asin(sin(newtheta) * sqrt(bSource / field()));
                if (newtheta > PI / 2) {
                    theta = PI - theta;
                }
            }

            assert !Double.isNaN(theta);

            return theta;
        }

    }


}
