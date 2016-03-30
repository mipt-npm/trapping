/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hep.dataforge.trapping;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.SynchronizedRandomGenerator;

/**
 *
 * @author Darksnake
 */
public class TrappingRandomGenerator {
    RandomGenerator generator;

    public TrappingRandomGenerator() {
        this.generator = new SynchronizedRandomGenerator(new MersenneTwister());
    }

    public TrappingRandomGenerator(RandomGenerator generator) {
        this.generator = generator;
    }
    
    /**
     * heads-tails random.
     * @return 
     */
    public boolean heads(){
        return generator.nextBoolean();
    }
    
    /**
     * next uniform in [0;1] 
     * @return 
     */
    public double next(){
        return generator.nextDouble();
    }
    
    
    
}
