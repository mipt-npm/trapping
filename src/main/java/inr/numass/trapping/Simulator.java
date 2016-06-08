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

import java.util.function.Function;

import static org.apache.commons.math3.util.FastMath.*;

/**
 * @author Darksnake
 */
public class Simulator {

    public static final double SOURCE_LENGTH = 3d;
    private static final double DELTA_L = 0.1; //step for dZ calculation

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

    public static enum EndState {

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

    public void setGasDensity(double gasDensity) {
        this.gasDensity = gasDensity;
    }

    public void setFieldFunc(UnivariateFunction magneticField) {
        this.magneticField = magneticField;
    }

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
        //FIXME double cross-section calculation
        //All cross-sections are in m^2
        return new ExponentialDistribution(generator, 1 / scatter.sigmaTotal(e) / gasDensity).sample();
    }

    /**
     * @param position
     * @param dl
     * @return
     */
    private State propagate(State position, double dl, double dz) {
        //increase l
        position.l += dl;
        position.addZ(dz);
        return position;
    }

    /**
     * calculate z coordinate change with known path length. Does not change position.
     *
     * @param dl
     * @return
     */
    private double dZ(State position, double dl) {
        // if magnetic field not defined, consider it to be uniform and equal bSource
        if (magneticField == null) {
            return dl / cos(position.theta);
        } else {
            double dz = 0;
            double curL = 0;
            while (curL <= dl) {
                double delta = min(dl - curL, DELTA_L);
                double b = field(position.z + dz);

                dz += delta / sqrt(1 - pow(sin(position.theta), 2) * b / bSource);
                curL += DELTA_L;
            }
            return dz;
        }
    }

    /**
     * Magnetic field with reflection taken into account
     *
     * @param z
     * @return
     */
    private double field(double z) {
        if (magneticField == null) {
            return bSource;
        } else {
            while (!(z > 0 && z < SOURCE_LENGTH)) {
                // reflecting from back wall
                if (z < 0) {
                    z = -z;
                } else {
                    z -= SOURCE_LENGTH;
                }
            }
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
        assert initZ >= 0 && initZ < SOURCE_LENGTH;

//        if (initTheta < this.thetaPinch) {
//            if (generator.nextBoolean()) {
//                return new SimulationResult(EndState.PASS, initEnergy, initTheta, initTheta, 0);
//            } else {
//                return new SimulationResult(EndState.REJECTED, initEnergy, initTheta, initTheta, 0);
//            }
//        } else if (initTheta < this.thetaTransport) {
//            return new SimulationResult(EndState.REJECTED, initEnergy, initTheta, initTheta, 0);
//        }

        State pos = new State(initEnergy, initTheta, initZ);
        EndState endState = EndState.NONE;

        boolean stopflag = false;

        while (!stopflag) {

            double dl = freePath(pos.e); // path to next scattering
            double dz = dZ(pos, dl); // z coordinate to next scattering
            double expectedZ = pos.z + dz; // expected z position of next scattering

            //if no scattering on current source pass
            if (expectedZ < 0 || expectedZ > SOURCE_LENGTH) {
                //accepted by spectrometer
                if (pos.theta < thetaPinch) {
                    stopflag = true;
                    //Учитываем тот факт, что электрон мог вылететь в правильный угол, но назад
                    if (pos.colNum == 0) {
                        //counting pass electrons
                        endState = EndState.PASS;
                    } else {
                        endState = EndState.ACCEPTED;
                    }
                }


                //through the rear magnetic trap
                if (pos.theta >= PI - thetaTransport) {
                    stopflag = true;
                    endState = EndState.REJECTED;
                }

                if (pos.e < eLow) {
                    //Если энергия стала слишком маленькой
                    stopflag = true;
                    endState = EndState.LOWENERGY;
                }
            }
            if (!stopflag) {
                // perform scatter
                propagate(pos, dl, dz);
                pos.colNum++;
                scatter(pos);
            }

        }

        SimulationResult res = new SimulationResult(endState, pos.e, pos.theta, initTheta, pos.colNum, pos.l);
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
         * current z;
         */
        double z;

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

        /**
         * add Z and calculate direction change
         *
         * @param dZ
         * @return
         */
        double addZ(double dZ) {
            this.z += dZ;
            while (!(this.z > 0 && this.z < SOURCE_LENGTH)) {
                // reflecting from back wall
                if (z < 0) {
                    z = -z;
                    flip();
                } else {
                    z -= SOURCE_LENGTH;
                    flip();
                }
            }
            return z;
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
                return asin(sin(theta) * sqrt(field() / bSource));
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

            //следим чтобы угол был от 0 до Pi
            if (newtheta < 0) {
                newtheta = -newtheta;
            }
            if (newtheta > Math.PI && newtheta <= Math.PI * 3 / 2) {
                newtheta = 2 * Math.PI - newtheta;
            }

            //change back to virtual angles
            if (magneticField == null) {
                theta = newtheta;
            } else {
                theta = asin(sin(newtheta) * sqrt(bSource / field()));
            }
            return theta;
        }

    }


}
