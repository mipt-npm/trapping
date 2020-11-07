package ru.inr.mass.trapping

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator
import org.apache.commons.math3.random.JDKRandomGenerator
import org.apache.commons.rng.UniformRandomProvider
import org.apache.commons.rng.simple.RandomSource
import java.io.File
import java.io.PrintStream
import java.nio.file.Files
import java.nio.file.StandardOpenOption
import java.util.stream.LongStream

private val seedGenerator = JDKRandomGenerator()


/**
 * Created by darksnake on 04-Jun-16.
 */
class SimulationManager() {

    var outputDirectory: File = File("./output")
    var fileName = "trap[test]"

    var comment = ""

    /**
     * A supplier for random generator. Each track has its own generator
     */
    var generatorFactory: (Long) -> UniformRandomProvider =
        { RandomSource.create(RandomSource.MT_64, seedGenerator.nextInt()) }
    var reportFilter: (Simulator.SimulationResult) -> Boolean = { it.state == Simulator.EndState.ACCEPTED }

    var initialE = 18000.0
    var eLow: Double = 14000.0
    var thetaTransport: Double = 24.107064 / 180 * Math.PI
    var thetaPinch: Double = 19.481097 / 180 * Math.PI
    var gasDensity: Double = 5e20// m^-3
    var bSource: Double = 0.6
    var magneticField: ((Double) -> Double)? = null

    /**
     * A syntentic property which allows to set lower energy by specifying range instead of energy itself
     */
    var range: Double
        get() = initialE - eLow
        set(value) {
            eLow = initialE - value
        }


//    fun setOutputFile(fileName: String) {
//        val outputFile = File(fileName)
//        if (!outputFile.exists()) {
//            outputFile.parentFile.mkdirs()
//            outputFile.createNewFile()
//        }
//        this.output = PrintStream(FileOutputStream(outputFile))
//    }


    fun withReportFilter(filter: (Simulator.SimulationResult) -> Boolean): SimulationManager {
        this.reportFilter = filter
        return this
    }

    fun setFields(Bsource: Double, Btransport: Double, Bpinch: Double) {
        this.bSource = Bsource
        this.thetaTransport = Math.asin(Math.sqrt(Bsource / Btransport))
        this.thetaPinch = Math.asin(Math.sqrt(Bsource / Bpinch))
    }

    /**
     * Set field map from values
     *
     * @param z
     * @param b
     * @return
     */
    fun setFieldMap(z: DoubleArray, b: DoubleArray) {
        this.magneticField = { LinearInterpolator().interpolate(z, b).value(it) }
    }

    /**
     * Симулируем пролет num электронов.
     *
     * @param num
     * @return
     */
    @Synchronized
    fun simulateAll(num: Number): Counter {
        if (!outputDirectory.exists()) {
            outputDirectory.mkdirs()
        }


        val simulator = Simulator(eLow, thetaTransport, thetaPinch, gasDensity, bSource, magneticField)

        val counter = Counter()

        val header = """
                E_init = $initialE;
                E_low = $eLow;
                theta_pinch = $thetaPinch;
                theta_transport = $thetaTransport;
                density = $gasDensity;
            """.trimIndent() + "\n\n" + comment


        val outputPath = outputDirectory.toPath().resolve("$fileName.out")

        PrintStream(
            Files.newOutputStream(
                outputPath,
                StandardOpenOption.CREATE,
                StandardOpenOption.TRUNCATE_EXISTING,
                StandardOpenOption.WRITE
            )
        ).use { output ->
            output.println("# " + header.replace("\n", "\n# "))//adding comment symbols

            System.out.printf("%nStarting simulation with initial energy %g and %d electrons.%n%n",
                initialE,
                num.toLong())
            output.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s%n", "id", "E", "theta", "theta_start", "colNum", "L", "state")
            LongStream.rangeClosed(1, num.toLong()).parallel().mapToObj { id ->
                val generator = RandomGeneratorBridge(generatorFactory(id))
                val theta = Math.acos(1 - 2 * generator.nextDouble())// from 0 to Pi
                val z = (generator.nextDouble() - 0.5) * Simulator.SOURCE_LENGTH
                simulator.simulate(id, generator, initialE, theta, z).also { counter.count(it) }
            }.filter(reportFilter).forEach { res ->
                printOne(output, res)
            }
        }

        val statisticsPath = outputDirectory.toPath().resolve("$fileName.stat")
        PrintStream(Files.newOutputStream(statisticsPath,
            StandardOpenOption.CREATE,
            StandardOpenOption.TRUNCATE_EXISTING,
            StandardOpenOption.WRITE)).use {
            it.println(header + "\n")
            printStatistics(it, simulator, counter)

        }
        return counter
    }

    private fun printStatistics(output: PrintStream, simulator: Simulator, counter: Counter) {
        output.apply {
            println()
            println("***RESULT***")
            printf("The total number of events is %d.%n%n", counter.count)
            printf("The spectrometer acceptance angle is %g.%n", simulator.thetaPinch * 180 / Math.PI)
            printf("The transport reflection angle is %g.%n", simulator.thetaTransport * 180 / Math.PI)
            printf("The starting energy is %g.%n", initialE)
            printf("The lower energy boundary is %g.%n", simulator.eLow)
            printf("The source density is %g.%n", simulator.gasDensity)

            println()
            printf("The total number of ACCEPTED events is %d.%n", counter.accepted)
            printf("The total number of PASS events is %d.%n", counter.pass)
            printf("The total number of REJECTED events is %d.%n", counter.rejected)
            printf("The total number of LOWENERGY events is %d.%n%n", counter.lowE)
        }

    }

    private fun printOne(out: PrintStream, res: Simulator.SimulationResult) {
        out.printf("%d\t%g\t%g\t%g\t%d\t%g\t%s%n",
            res.id,
            res.E,
            res.theta * 180 / Math.PI,
            res.initTheta * 180 / Math.PI,
            res.collisionNumber,
            res.l,
            res.state.toString())
        out.flush()
    }


    /**
     * Statistic counter
     */
    class Counter {
        internal var accepted = 0
        internal var pass = 0
        internal var rejected = 0
        internal var lowE = 0
        internal var count = 0

        fun count(res: Simulator.SimulationResult): Simulator.SimulationResult {
            synchronized(this) {
                count++
                when (res.state) {
                    Simulator.EndState.ACCEPTED -> accepted++
                    Simulator.EndState.LOWENERGY -> lowE++
                    Simulator.EndState.PASS -> pass++
                    Simulator.EndState.REJECTED -> rejected++
                    Simulator.EndState.NONE -> throw IllegalStateException()
                }
            }
            return res
        }

        @Synchronized
        fun reset() {
            count = 0
            accepted = 0
            pass = 0
            rejected = 0
            lowE = 0
        }

    }
}
