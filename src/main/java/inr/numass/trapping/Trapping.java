package inr.numass.trapping;

import java.io.FileNotFoundException;

/**
 * Hello world!
 */
public class Trapping {
    public static void main(String[] args) throws FileNotFoundException {
//        new SimulationManager().withParameters(0.6, 3.7, 7.2, 18000d, 4000).simulateAll((int) 1e6);
        new SimulationManager().withParameters(0.6, 3.7, 4.84, 14000d, 4000).simulateAll((int) 1e6);
    }
}
