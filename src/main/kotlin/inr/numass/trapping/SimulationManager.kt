package inr.numass.trapping

import org.apache.commons.math3.analysis.UnivariateFunction
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator
import org.apache.commons.math3.random.JDKRandomGenerator
import org.apache.commons.math3.random.RandomGenerator

import java.io.*
import java.util.stream.Stream

/**
 * Created by darksnake on 04-Jun-16.
 */
class SimulationManager {

    var generator: RandomGenerator = JDKRandomGenerator()

    /**
     *  output for accepted events
     */
    var output: PrintStream = System.out
    /**
     * output for statistics
     */
    var statisticOutput = System.out
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


    // from 0 to Pi
    private fun getRandomTheta(): Double {
        val x = generator.nextDouble()
        return Math.acos(1 - 2 * x)
    }

    fun setOutputFile(fileName: String) {
        val outputFile = File(fileName)
        if (!outputFile.exists()) {
            outputFile.parentFile.mkdirs()
            outputFile.createNewFile()
        }
        this.output = PrintStream(FileOutputStream(outputFile))
    }


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
    fun simulateAll(num: Int): Counter {
        val counter = Counter()
        val simulator = Simulator(eLow, thetaTransport, thetaPinch, gasDensity, bSource, magneticField, generator)

        System.out.printf("%nStarting sumulation with initial energy %g and %d electrons.%n%n", initialE, num)
        output.printf("%s\t%s\t%s\t%s\t%s\t%s%n", "E", "theta", "theta_start", "colNum", "L", "state")
        Stream.generate { getRandomTheta() }.limit(num.toLong()).parallel()
                .forEach { theta ->
                    val initZ = (generator.nextDouble() - 0.5) * Simulator.SOURCE_LENGTH
                    val res = simulator.simulate(initialE, theta, initZ)
                    if (reportFilter(res)) {
                        printOne(output, res)
                    }
                    counter.count(res)
                }
        printStatistics(statisticOutput, simulator, counter)
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
        out.printf("%g\t%g\t%g\t%d\t%g\t%s%n", res.E, res.theta * 180 / Math.PI, res.initTheta * 180 / Math.PI, res.collisionNumber, res.l, res.state.toString())
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
