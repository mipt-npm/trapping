package inr.numass.trapping;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.random.ISAACRandom;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.SynchronizedRandomGenerator;

import java.io.PrintStream;
import java.util.Map;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Stream;

/**
 * Created by darksnake on 04-Jun-16.
 */
public class SimulationManager {

    RandomGenerator generator = new JDKRandomGenerator();
    Simulator simulator = new Simulator();
    private double initialE = 18000;
    private PrintStream output = System.out;
    private PrintStream statisticOutput = System.out;


    public SimulationManager withParameters(double bSource, double bTransport, double bPinch, double initialE, double energyRange) {
        this.simulator = new Simulator(bSource, bTransport, bPinch, initialE - energyRange);
        this.initialE = initialE;
        return this;
    }

    public SimulationManager withOutput(PrintStream output) {
        this.output = output;
        return this;
    }

    public SimulationManager withStatisticOutput(PrintStream statisticOutput) {
        this.statisticOutput = statisticOutput;
        return this;
    }

    public SimulationManager withFieldMap(UnivariateFunction fieldMap) {
        this.simulator.setFieldFunc(fieldMap);
        return this;
    }

    public SimulationManager withFieldMap(double[] z, double[] b) {
        this.simulator.setFieldFunc(new LinearInterpolator().interpolate(z, b));
        return this;
    }

    public SimulationManager withDensity(double density){
        this.simulator.setGasDensity(density);
        return this;
    }

    /**
     * Симулируем пролет num электронов.
     *
     * @param num
     * @return
     */
    public synchronized Counter simulateAll(int num) {
        Counter counter = new Counter();
        Predicate<Simulator.SimulationResult> reportIf = (res) -> res.state == Simulator.EndState.ACCEPTED;
        System.out.printf("%nStarting sumulation with initial energy %g and %d electrons.%n%n", initialE, num);
        Stream.generate(() -> getRandomTheta()).limit(num).parallel()
                .forEach((theta) -> {
                    double initZ = generator.nextDouble() * Simulator.SOURCE_LENGTH;
                    Simulator.SimulationResult res = simulator.simulate(initialE, theta, initZ);
                    if (reportIf.test(res)) {
                        if (output != null) {
                            printOne(output, res);
                        }
                    }
                    counter.count(res);
                });
        printStatistics(counter);
        return counter;
    }

    private double getRandomTheta() {
        double x = generator.nextDouble();
        // from 0 to 2 Pi
        return Math.acos(1 - 2 * x);
    }

    private void printStatistics(Counter counter) {
        output.printf("The total number of events is %d.%n%n", counter.count);

        output.printf("The spectrometer acceptance angle is %g.%n", simulator.thetaPinch * 180 / Math.PI);
        output.printf("The transport reflection angle is %g.%n", simulator.thetaTransport * 180 / Math.PI);
        output.printf("The starting energy is %g.%n", initialE);
        output.printf("The lower energy boundary is %g.%n%n", simulator.eLow);

        output.printf("The total number of ACCEPTED events is %d.%n", counter.accepted);
        output.printf("The total number of PASS events is %d.%n", counter.pass);
        output.printf("The total number of REJECTED events is %d.%n", counter.rejected);
        output.printf("The total number of LOWENERGY events is %d.%n%n", counter.lowE);
    }

    private void printOne(PrintStream out, Simulator.SimulationResult res) {
        out.printf("%g\t%g\t%g\t%d\t%s%n", res.E, res.theta * 180 / Math.PI, res.initTheta * 180 / Math.PI, res.collisionNumber, res.state.toString());
    }


    public static class Counter {
        int accepted = 0;
        int pass = 0;
        int rejected = 0;
        int lowE = 0;
        int count = 0;

        public synchronized Simulator.SimulationResult count(Simulator.SimulationResult res) {
            count++;
            switch (res.state) {
                case ACCEPTED:
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
            return res;
        }

        public synchronized void reset() {
            count = 0;
            accepted = 0;
            pass = 0;
            rejected = 0;
            lowE = 0;
        }

    }
}
