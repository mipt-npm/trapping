package inr.numass.trapping;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.time.Duration;
import java.time.Instant;

public class Trapping {
    public static void main(String[] args) throws IOException {
//        new SimulationManager().withParameters(0.6, 3.7, 7.2, 18000d, 4000).simulateAll((int) 1e6);

//        -0.236	3.70754
//        0.23	0.62786
//        0.746	0.60474
//        1.262	0.60325
//        1.778	0.60333
//        2.294	0.60503
//        2.81	0.6285
//        3.276	3.70478
        double[] z = {-0.236, 0.23, 0.746, 1.262, 1.778, 2.294, 2.81, 3.276};
        double[] b = {3.70754, 0.62786, 0.60474, 0.60325, 0.60333, 0.60503, 0.6285, 3.70478};
        System.out.println("Press any key to start...");
        System.in.read();
        Instant startTime = Instant.now();
        System.out.printf("Starting at %s%n%n", startTime.toString());
        new SimulationManager()
                .withParameters(0.6, 3.7, 4.84, 18000d, 4000)
//                .withFieldMap(z, b)
//                .withDensity(5e20)
                .simulateAll((int) 1e5);
        Instant finishTime = Instant.now();
        System.out.printf("%nFinished at %s%n", finishTime.toString());
        System.out.printf("Calculation took %s%n", Duration.between(startTime,finishTime).toString());
    }
}
