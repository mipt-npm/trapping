package inr.numass.trapping;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.time.Duration;
import java.time.Instant;

public class Trapping {
    public static void main(String[] args) throws IOException {
//        new SimulationManager().withParameters(0.6, 3.7, 7.2, 18000d, 4000).simulateAll((int) 1e6);
        double[] z = {-1.736, -1.27, -0.754, -0.238, 0.278, 0.794, 1.31, 1.776};
        double[] b = {3.70754, 0.62786, 0.60474, 0.60325, 0.60333, 0.60503, 0.6285, 3.70478};
        System.out.println("Press any key to start...");
        System.in.read();
        Instant startTime = Instant.now();
        System.out.printf("Starting at %s%n%n", startTime.toString());
        new SimulationManager()
                .withParameters(0.6, 3.7, 4.84, 18000d, 4000)
//                .withOutputFile("D:\\Work\\Numass\\trapping\\trap 18, pinch 100A, fields.out")
//                .withFieldMap(z, b)
//                .withDensity(5e20)
                .simulateAll((int) 1e4);
        Instant finishTime = Instant.now();
        System.out.printf("%nFinished at %s%n", finishTime.toString());
        System.out.printf("Calculation took %s%n", Duration.between(startTime,finishTime).toString());
    }
}
