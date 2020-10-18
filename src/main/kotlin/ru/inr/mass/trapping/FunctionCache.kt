package ru.inr.mass.trapping

import java.util.*

class FunctionCache(private val xPrecision: Double, val function: (Double) -> Double) {
    private val values = TreeMap<Double, Double>()

    operator fun invoke(x: Double): Double {
        val floor: MutableMap.MutableEntry<Double, Double>? = values.floorEntry(x)
        val ceil: MutableMap.MutableEntry<Double, Double>? = values.ceilingEntry(x)
        return if (floor == null || ceil == null || ceil.key - floor.key <= xPrecision) {
            val newValue = function(x)
            values[x] = newValue
            newValue
        } else {
            val x0 = floor.key
            val x1 = ceil.key
            val y0 = floor.value
            val y1 = ceil.value
            y0 + (x - x0) * (y1 - y0) / (x1 - x0)
        }
    }
}