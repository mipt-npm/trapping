/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hep.dataforge.trapping;

import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomGenerator;

/**
 *
 * @author Darksnake
 */
public class TrappingRandomGenerator {
    RandomGenerator generator;

    public TrappingRandomGenerator() {
        this.generator = new JDKRandomGenerator();
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
