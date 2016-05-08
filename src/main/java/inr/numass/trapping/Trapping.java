package inr.numass.trapping;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.List;

/**
 * Hello world!
 *
 */
public class Trapping {

    public static void main(String[] args) throws FileNotFoundException {
        PrintStream out = null;
        if ((args.length > 0) && (args[0] != null)) {
            File file = new File(args[0]);
            out = new PrintStream(file);
        } else {

        }

        double E = 18000d;
        System.out.println();
        ElectronTrappingSimulator simulator = new ElectronTrappingSimulator();
        simulator.setFields(0.6, 3.6, 7.2);

        int accepted = 0;
        int pass = 0;
        int rejected = 0;
        int lowE = 0;

        simulator.scatter.counter.resetAll();
        List<ElectronTrappingSimulator.SimulaionResult> results = simulator.simulateAll(E, (int) 1e6);
        simulator.scatter.counter.print(System.out);

        System.out.printf("%nSimulation complete.%n%n");

        for (ElectronTrappingSimulator.SimulaionResult res : results) {
            switch (res.state) {
                case ACCEPTED:
//                    ElectronTrappingSimulator.printOne(System.out, res);
                    if (out != null) {
                        ElectronTrappingSimulator.printOne(out, res);
                    }
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
        System.out.printf("The transport reflection angle is %g.%n", simulator.thetaTransport * 180 / Math.PI);
        System.out.printf("The starting energy is %g.%n", E);
        System.out.printf("The lower energy boundary is %g.%n%n", simulator.Elow);

        System.out.printf("The total number of ACCEPTED events is %d.%n", accepted);
        System.out.printf("The total number of PASS events is %d.%n", pass);
        System.out.printf("The total number of REJECTED events is %d.%n", rejected);
        System.out.printf("The total number of LOWENERGY events is %d.%n%n", lowE);

        if (out != null) {
            out.println();
            out.printf("The total number of events is %d.%n%n", results.size());

            out.printf("The spectrometer acceptance angle is %g.%n", simulator.thetaPinch * 180 / Math.PI);
            out.printf("The transport reflection angle is %g.%n", simulator.thetaTransport * 180 / Math.PI);
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