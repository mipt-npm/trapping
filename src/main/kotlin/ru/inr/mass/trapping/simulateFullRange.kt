package ru.inr.mass.trapping

import java.time.Duration
import java.time.Instant


fun main() {
    //val z = doubleArrayOf(-1.736, -1.27, -0.754, -0.238, 0.278, 0.794, 1.31, 1.776)
    //val b = doubleArrayOf(3.70754, 0.62786, 0.60474, 0.60325, 0.60333, 0.60503, 0.6285, 3.70478)
    //        System.out.println("Press any key to start...");
    //        System.in.read();

    val startTime = Instant.now()
    System.out.printf("Starting at %s%n%n", startTime.toString())

    SimulationManager().apply {
        comment = "Out of the box cross-sections"
        fileName = "trap-full-range"
        setFields(0.6, 3.6, 7.2)
        gasDensity = 1e19 // m^-3
        initialE = 18500.0
        range = 16000.0
    }.simulateAll(1_000_001)

    val finishTime = Instant.now()
    System.out.printf("%nFinished at %s%n", finishTime.toString())
    System.out.printf("Calculation took %s%n", Duration.between(startTime, finishTime).toString())
}

