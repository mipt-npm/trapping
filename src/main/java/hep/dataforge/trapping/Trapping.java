package hep.dataforge.trapping;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;

/**
 * Hello world!
 *
 */
public class Trapping {

    public static void main(String[] args) throws FileNotFoundException {
        PrintWriter out = null;
        if ((args != null) && (args[0] != null)) {
            File file = new File(args[0]);
            out = new PrintWriter(file);
        }

        double E = 18000d;
        System.out.println();
//        System.setProperty("jna.library.path", "d:\\projects\\Trapping\\target\\classes\\win32-amd64\\libScatter.dll");
        ElectronTrappingSimulator simulator = new ElectronTrappingSimulator();
        simulator.setFields(0.6, 3.6, 7.2);
//        ElectronTrappingSimulator.SimulaionResult result;

        int accepted = 0;
        int pass = 0;
        int rejected = 0;
        int lowE = 0;

        ArrayList<ElectronTrappingSimulator.SimulaionResult> results = simulator.simulateAll(E, 100000);

        for (Iterator<ElectronTrappingSimulator.SimulaionResult> it = results.iterator(); it.hasNext();) {
            ElectronTrappingSimulator.SimulaionResult res = it.next();

            if (out != null) {
                out.printf("%g\t%g\t%d\t%s%n", res.E, res.theta * 180 / Math.PI, res.collisionNumber, res.state.toString());
            }
            switch (res.state) {

                case ACCEPTED:
                    System.out.printf("%g\t%g\t%d\t%s%n", res.E, res.theta * 180 / Math.PI, res.collisionNumber, res.state.toString());
                    accepted++;
                    break;
                case LOWENERGY:
                    lowE++;
                    break;
                case PASS:
                    pass++;
                    break;
                case REJECTED:
                    rejected++;
                    break;
                case NONE:
                    throw new IllegalStateException();
            }
        }

        System.out.printf("%nThe total number of events is %d.%n%n", results.size());

        System.out.printf("The spectrometer acceptance angle is %g.%n", simulator.thetaPinch * 180 / Math.PI);
        System.out.printf("The transport mirroring angle is %g.%n", simulator.thetaTransport * 180 / Math.PI);
        System.out.printf("The spectrometer acceptance angle is %g.%n", simulator.thetaPinch * 180 / Math.PI);
        System.out.printf("The starting energy is %g.%n", E);
        System.out.printf("The lower energy boundary is %g.%n%n", simulator.Elow);

        System.out.printf("The total number of ACCEPTED events is %d.%n", accepted);
        System.out.printf("The total number of PASS events is %d.%n", pass);
        System.out.printf("The total number of REJECTED events is %d.%n", rejected);
        System.out.printf("The total number of LOWENERGY events is %d.%n%n", lowE);

        if (out != null) {
            out.printf("The total number of events is %d.%n%n", results.size());

            out.printf("The spectrometer acceptance angle is %g.%n", simulator.thetaPinch * 180 / Math.PI);
            out.printf("The transport mirroring angle is %g.%n", simulator.thetaTransport * 180 / Math.PI);
            out.printf("The spectrometer acceptance angle is %g.%n", simulator.thetaPinch * 180 / Math.PI);
            out.printf("The starting energy is %g.%n", E);
            out.printf("The lower energy boundary is %g.%n%n", simulator.Elow);

            out.printf("The total number of ACCEPTED events is %d.%n", accepted);
            out.printf("The total number of PASS events is %d.%n", pass);
            out.printf("The total number of REJECTED events is %d.%n", rejected);
            out.printf("The total number of LOWENERGY events is %d.%n%n", lowE);
            out.close();
        }
        

    }
}
