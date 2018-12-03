package inr.numass.trapping

import org.apache.commons.math3.random.RandomGenerator
import org.apache.commons.rng.UniformRandomProvider

class RandomGeneratorBridge(val provider: UniformRandomProvider) : RandomGenerator {
    override fun nextBoolean(): Boolean = provider.nextBoolean()

    override fun nextFloat(): Float = provider.nextFloat()

    override fun setSeed(seed: Int) = error("Not supported")

    override fun setSeed(seed: IntArray?) = error("Not supported")

    override fun setSeed(seed: Long) = error("Not supported")

    override fun nextBytes(bytes: ByteArray) = provider.nextBytes(bytes)

    override fun nextInt(): Int = provider.nextInt()

    override fun nextInt(n: Int): Int = provider.nextInt(n)

    override fun nextGaussian(): Double = error("Not supported")

    override fun nextDouble(): Double = provider.nextDouble()

    override fun nextLong(): Long = provider.nextLong()
}