/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package inr.numass.trapping;

import java.io.PrintStream;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.SphericalCoordinates;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.SynchronizedRandomGenerator;
import org.apache.commons.math3.util.Pair;
import org.apache.commons.math3.util.Precision;

/**
 *
 * @author Darksnake
 */
public class ElectronTrappingSimulator {

    private RandomGenerator generator;
    Scatter scatter;
    double Elow = 14000d;
    double thetaTransport = 24.107064 / 180 * Math.PI;
    double thetaPinch = 19.481097 / 180 * Math.PI;


    public ElectronTrappingSimulator() {
        generator = new SynchronizedRandomGenerator(new MersenneTwister());
        scatter = new Scatter(generator);
    }

    public static enum EndState {

        ACCEPTED,//трэппинговый электрон попал в аксептанс
        REJECTED,//трэппинговый электрон вылетел через заднюю пробку
        LOWENERGY,//потерял слишком много энергии
        PASS,//электрон никогда не запирался и прошел напрямую, нужно для нормировки
        NONE
    }

    public void setFields(double Bsource, double Btransport, double Bpinch) {
        this.thetaTransport = Math.asin(Math.sqrt(Bsource / Btransport));
        this.thetaPinch = Math.asin(Math.sqrt(Bsource / Bpinch));
    }

    /**
     * Симулируем один пробег электрона от начального значения и до вылетания из
     * иточника или до того момента, как энергия становится меньше Elow.
     * Считаем, что за один проход источника происходит в среднем существенно
     * меньше одного столкновения, поэтому выход вперед или назад совершенно
     * симментричны.
     *
     */
    public SimulaionResult simulateOne(double initEnergy, double initTheta) {
        assert initEnergy > 0;
        assert initTheta > 0 && initTheta < Math.PI / 2;

        if (initTheta < this.thetaPinch) {
            if (generator.nextBoolean()) {
                return new SimulaionResult(EndState.PASS, initEnergy, initTheta, initTheta, 0);
            } else {
                return new SimulaionResult(EndState.REJECTED, initEnergy, initTheta, initTheta, 0);
            }
        } else if (initTheta < this.thetaTransport) {
            return new SimulaionResult(EndState.REJECTED, initEnergy, initTheta, initTheta, 0);
        }

        double E = initEnergy;
        double theta = initTheta;
        int colNum = 0;
        EndState state = EndState.NONE;

        boolean stopflag = false;

        while (!stopflag) {
            colNum++;
            Pair<Double,Double> delta;

            //Вычисляем сечения и нормируем их на полное сечение
            double sigmaIon = scatter.sigmaion(E);
            double sigmaEl = scatter.sigmael(E);
            double sigmaexc = scatter.sigmaexc(E);
            double sigmaTotal = sigmaEl + sigmaIon + sigmaexc;
            sigmaIon /= sigmaTotal;
            sigmaEl /= sigmaTotal;
            sigmaexc /= sigmaTotal;
            //проверяем нормировку
            assert Precision.equals(sigmaEl + sigmaexc + sigmaIon, 1, 1e-2);

            double alpha = generator.nextDouble();

            if (alpha > sigmaEl) {
                if (alpha > sigmaEl + sigmaexc) {
                    //ionization case
                    delta = scatter.randomion(E);
                } else {
                    //excitation case
                    delta = scatter.randomexc(E);
                }
            } else {
                // elastic
                delta = scatter.randomel(E);
            }

            //Обновляем значени угла и энергии независимо ни от чего
            E -= delta.getFirst();
            //Изменение угла
            theta = addTheta(theta, delta.getSecond() / 180 * Math.PI);
            //следим чтобы угол был от 0 до 90, если он перекинется через границу, считаем что электрон остается в потоке
            theta = normalizeTheta(theta);

            if (theta < thetaPinch) {
                stopflag = true;
                if (generator.nextBoolean()) {
                    //Учитываем тот факт, что электрон мог вылететь в правильный угол, но назад
                    state = EndState.ACCEPTED;
                } else {
                    state = EndState.REJECTED;
                }
            }

            if (theta >= thetaPinch && theta <= thetaTransport) {
                //Вылет через заднюю пробку эквивалентен отражению от передней.
                //это верно при низких концентрациях
                stopflag = true;
                state = EndState.REJECTED;
            }

            if (E < Elow) {
                //Если энергия стала слишком маленькой
                stopflag = true;
                state = EndState.LOWENERGY;
            }
        }
        SimulaionResult res = new SimulaionResult(state, E, theta, initTheta, colNum);
        if (state == EndState.ACCEPTED) {
            printOne(System.out, res);
        }
        return res;

    }

    /**
     * Сложение вектора с учетом случайно распределения по фи
     *
     * @param theta
     * @param dTheta
     * @return
     */
    double addTheta(double theta, double dTheta) {
        //Генерируем случайный фи
        double phi = generator.nextDouble()* 2 * Math.PI;
        //Создаем начальный вектор в сферических координатах
        SphericalCoordinates init = new SphericalCoordinates(1, 0, theta + dTheta);
        // Задаем вращение относительно оси, перпендикулярной исходному вектору
        SphericalCoordinates rotate = new SphericalCoordinates(1, 0, theta);
        // поворачиваем исходный вектор на dTheta
        Rotation rot = new Rotation(rotate.getCartesian(), phi);

        Vector3D result = rot.applyTo(init.getCartesian());

        //      assert Vector3D.angle(result, rotate.getCartesian()) == dTheta;
        return Math.acos(result.getZ());
    }

    /**
     * Симулируем пролет num электронов.
     *
     * @param E
     * @param num
     * @return
     */
    public List<SimulaionResult> simulateAll(double E, int num) {
        System.out.printf("%nStarting sumulation with initial energy %g and %d electrons.%n%n", E, num);
        return Stream.generate(() -> getRandomTheta()).limit(num).parallel()
                .map(theta -> simulateOne(E, theta))
                .collect(Collectors.toList());
    }

    public static void printOne(PrintStream out, SimulaionResult res) {
        out.printf("%g\t%g\t%g\t%d\t%s%n", res.E, res.theta * 180 / Math.PI, res.initTheta * 180 / Math.PI, res.collisionNumber, res.state.toString());
    }

    private double normalizeTheta(double theta) {
        double res = theta;
        if (theta < 0) {
            res = -theta;
        }
        if (res >= 0 && res <= Math.PI / 2) {
            return res;
        } else if (res > Math.PI / 2 && res <= Math.PI) {
            return Math.PI - res;
        } else if (res > Math.PI && res <= Math.PI * 3 / 2) {
            return res - Math.PI;
        } else {
            throw new IllegalStateException();
        }
    }

    public double getRandomTheta() {
        double x = generator.nextDouble();
        return Math.acos(x);
    }

    public class SimulaionResult {

        public SimulaionResult(EndState state, double E, double theta, double initTheta, int collisionNumber) {
            this.state = state;
            this.E = E;
            this.theta = theta;
            this.initTheta = initTheta;
            this.collisionNumber = collisionNumber;
        }
        public EndState state;
        public double E;
        public double initTheta;
        public double theta;
        public int collisionNumber;
    }
}
