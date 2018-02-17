package inr.numass.trapping;

import java.io.IOException;
import java.time.Duration;
import java.time.Instant;
import java.util.function.Predicate;

public class Trapping {
    public static void main(String[] args) throws IOException {
        double[] z = {-1.736, -1.27, -0.754, -0.238, 0.278, 0.794, 1.31, 1.776};
        double[] b = {3.70754, 0.62786, 0.60474, 0.60325, 0.60333, 0.60503, 0.6285, 3.70478};
//        System.out.println("Press any key to start...");
//        System.in.read();
        Instant startTime = Instant.now();
        System.out.printf("Starting at %s%n%n", startTime.toString());
        new SimulationManager()
                .withOutputFile("D:\\Work\\Numass\\trapping\\test1.out")
                .withFields(0.6, 3.7, 7.2)
//                .withFieldMap(z, b)
                .withGasDensity(1e19) // per m^3
                .withReportFilter(res -> true)
                .simulateAll((int) 1e4);
        Instant finishTime = Instant.now();
        System.out.printf("%nFinished at %s%n", finishTime.toString());
        System.out.printf("Calculation took %s%n", Duration.between(startTime, finishTime).toString());
    }
}
