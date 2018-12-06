package inr.numass.trapping

import org.apache.commons.rng.simple.RandomSource
import java.io.File
import java.time.Duration
import java.time.Instant


fun main(args: Array<String>) {
    //val z = doubleArrayOf(-1.736, -1.27, -0.754, -0.238, 0.278, 0.794, 1.31, 1.776)
    //val b = doubleArrayOf(3.70754, 0.62786, 0.60474, 0.60325, 0.60333, 0.60503, 0.6285, 3.70478)
    //        System.out.println("Press any key to start...");
    //        System.in.read();

    val es = listOf(12000.0, 14000.0, 16000.0, 18000.0)

    val startTime = Instant.now()
    System.out.printf("Starting at %s%n%n", startTime.toString())

    for (e in es) {
        SimulationManager().apply {
            comment = "Out of the box cross-sections"
            fileName = "trap[$e]"
            generator = RandomSource.create(RandomSource.MWC_256)
            setFields(0.6, 3.7, 7.2)
            gasDensity = 1e19
            initialE = e
            range = 4000.0
        }.simulateAll(1e7)
    }

    val finishTime = Instant.now()
    System.out.printf("%nFinished at %s%n", finishTime.toString())
    System.out.printf("Calculation took %s%n", Duration.between(startTime, finishTime).toString())
}

