/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hep.dataforge.trapping;

import java.util.ArrayList;
import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.SphericalCoordinates;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.util.Precision;

/**
 *
 * @author Darksnake
 */
public class ElectronTrappingSimulator {

    private TrappingRandomGenerator generator = new TrappingRandomGenerator();
    double Elow = 14000d;
    double thetaTransport = 24.107064 / 180 * Math.PI;
    double thetaPinch = 19.481097 / 180 * Math.PI;

    public ElectronTrappingSimulator() {
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
            if (generator.heads()) {
                return new SimulaionResult(EndState.PASS, initEnergy, initTheta, 0);
            } else {
                return new SimulaionResult(EndState.REJECTED, initEnergy, initTheta, 0);
            }
        } else if(initTheta<this.thetaTransport){
                return new SimulaionResult(EndState.REJECTED, initEnergy, initTheta, 0);            
        }


        double E = initEnergy;
        double theta = initTheta;
        int colNum = 0;
        EndState state = EndState.NONE;

        boolean stopflag = false;

        while (!stopflag) {
            colNum++;
            DoubleValue dE = new DoubleValue(0);
            DoubleValue dTheta = new DoubleValue(0);


            //Вычисляем сечения и нормируем их на полное сечение
            double sigmaIon = Scatter.sigmaion(E);
            double sigmaEl = Scatter.sigmael(E);
            double sigmaexc = Scatter.sigmaexc(E);
            double sigmaTotal = sigmaEl + sigmaIon + sigmaexc;
            sigmaIon /= sigmaTotal;
            sigmaEl /= sigmaTotal;
            sigmaexc /= sigmaTotal;
            //проверяем нормировку
            assert Precision.equals(sigmaEl + sigmaexc + sigmaIon, 1, 1e-2);

            double alpha = generator.next();

            if (alpha > sigmaEl) {
                if (alpha > sigmaEl + sigmaexc) {
                    //ionization case
                    Scatter.randomion(E, dE, dTheta);
                } else {
                    //excitation case
                    Scatter.randomexc(E, dE, dTheta);
                }
            } else {
                // elastic
                Scatter.randomel(E, dE, dTheta);
            }

            //Обновляем значени угла и энергии независимо ни от чего
            E -= dE.getValue();
            //Изменение угла
            theta = addTheta(theta, dTheta.getValue() / 180 * Math.PI);
            //следим чтобы угол был от 0 до 90, если он перекинется через границу, считаем что электрон остается в потоке
            theta = normalizeTheta(theta);

            if (theta < thetaPinch) {
                stopflag = true;
                if (generator.heads()) {
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
        return new SimulaionResult(state, E, theta, colNum);

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
        double phi = generator.next() * 2 * Math.PI;
        //Создаем начальный вектор в сферических координатах
        SphericalCoordinates init = new SphericalCoordinates(1, 0, theta+dTheta);
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
    public ArrayList<SimulaionResult> simulateAll(double E, int num) {
        ArrayList<SimulaionResult> res = new ArrayList();
        double theta;

        System.out.printf("%nStarting sumulation with initial energy %g and %d electrons.%n%n", E, num);
        for (int i = 0; i < num; i++) {
//            System.out.printf("Running electron number %d", i);
            theta = this.getRandomTheta();
            res.add(this.simulateOne(E, theta));
        }
        return res;
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
        double x = generator.next();
        return Math.acos(x);
    }

    public class SimulaionResult {

        public SimulaionResult(EndState state, double E, double theta, int collisionNumber) {
            this.state = state;
            this.E = E;
            this.theta = theta;
            this.collisionNumber = collisionNumber;
        }
        public EndState state;
        public double E;
        public double theta;
        public int collisionNumber;
    }
}
