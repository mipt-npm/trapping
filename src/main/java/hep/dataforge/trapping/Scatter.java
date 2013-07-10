/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hep.dataforge.trapping;

import com.sun.jna.Library;
import com.sun.jna.Native;
import com.sun.jna.ptr.DoubleByReference;

/**
 *
 * @author Darksnake
 */
public class Scatter {

    public interface LibScatter extends Library {

        public void randomel(double E, DoubleByReference Eloss, DoubleByReference theta);

        public void randomexc(double E, DoubleByReference Eloss, DoubleByReference theta);

        public void randomion(double E, DoubleByReference Eloss, DoubleByReference theta);

        public double sigmael(double E);

        public double sigmaexc(double E);

        public double sigmaion(double E);
    }

    /**
     * PENDING переделать, чтобы возвращались нормальные значения
     * @param E
     * @param Eloss
     * @param theta 
     */
    static void randomel(double E, DoubleValue Eloss, DoubleValue theta) {
        LibScatter lib = (LibScatter) Native.loadLibrary("libScatter", LibScatter.class);

        DoubleByReference ElossPointer = new DoubleByReference(Eloss.getValue());
        DoubleByReference thetaPointer = new DoubleByReference(theta.getValue());
        lib.randomel(E, ElossPointer, thetaPointer);
        Eloss.setValue(ElossPointer.getValue());
        theta.setValue(thetaPointer.getValue());

    }

    static void randomexc(double E, DoubleValue Eloss, DoubleValue theta) {
        LibScatter lib = (LibScatter) Native.loadLibrary("libScatter", LibScatter.class);

        DoubleByReference ElossPointer = new DoubleByReference(Eloss.getValue());
        DoubleByReference thetaPointer = new DoubleByReference(theta.getValue());
        lib.randomexc(E, ElossPointer, thetaPointer);
        Eloss.setValue(ElossPointer.getValue());
        theta.setValue(thetaPointer.getValue());

    }

    static void randomion(double E, DoubleValue Eloss, DoubleValue theta) {
        LibScatter lib = (LibScatter) Native.loadLibrary("libScatter", LibScatter.class);

        DoubleByReference ElossPointer = new DoubleByReference(Eloss.getValue());
        DoubleByReference thetaPointer = new DoubleByReference(theta.getValue());
        lib.randomion(E, ElossPointer, thetaPointer);
        Eloss.setValue(ElossPointer.getValue());
        theta.setValue(thetaPointer.getValue());

    }

    /**
     * Все сечения в м^2
     * @param E
     * @return 
     */
    public static double sigmael(double E) {
        LibScatter lib = (LibScatter) Native.loadLibrary("libScatter", LibScatter.class);
        return lib.sigmael(E);
    }

    public static double sigmaexc(double E) {
        LibScatter lib = (LibScatter) Native.loadLibrary("libScatter", LibScatter.class);
        return lib.sigmaexc(E);
    }

    public static double sigmaion(double E) {
        LibScatter lib = (LibScatter) Native.loadLibrary("libScatter", LibScatter.class);
        return lib.sigmaion(E);
    }
    
    /**
     * Полное сечение с учетом квазиупругих столкновений
     * @param E
     * @return 
     */
    public static double sigmaTotal(double E){
        LibScatter lib = (LibScatter) Native.loadLibrary("libScatter", LibScatter.class);
        return lib.sigmael(E)+lib.sigmaexc(E)+lib.sigmaion(E);
    }
}
