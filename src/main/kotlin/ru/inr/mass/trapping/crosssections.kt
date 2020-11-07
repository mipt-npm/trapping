package ru.inr.mass.trapping

import de.m3y.kformat.Table
import de.m3y.kformat.table
import kscience.plotly.Plotly
import kscience.plotly.layout
import kscience.plotly.makeFile
import kscience.plotly.models.AxisType
import kscience.plotly.scatter

infix fun ClosedFloatingPointRange<Double>.step(step: Double): Sequence<Double> = sequence {
    var current = start

    while (current <= endInclusive) {
        yield(current)
        current += step
    }
}


/**
 * Extract cross-sections from the code
 */
fun main() {

    val energies = ((1.0..20.0) step 0.2).toList() // energy in keV

    val barn = 1e-22
    val sigmaEl = energies.map { Scatter.sigmael(it * 1000) / barn }
    val sigmaIon = energies.map { Scatter.sigmaion(it * 1000) / barn }
    val sigmaExc = energies.map { Scatter.sigmaexc(it * 1000) / barn }

    Plotly.plot {
        scatter {
            x.numbers = energies
            y.numbers = sigmaEl
            name = "Elastic"
        }
        scatter {
            x.numbers = energies
            y.numbers = sigmaIon
            name = "Ionization"
        }
        scatter {
            x.numbers = energies
            y.numbers = sigmaExc
            name = "Excitation"
        }
        layout{
            yaxis {
                title = "Cross-section (Barn)"
                type = AxisType.log
            }
            xaxis{
                title = "Electron energy (KeV)"
            }

        }
    }.makeFile()

    val table = table {
        header("E(keV)", "elastic", "ion", "exc")

        for (e in (12.0..18.6) step 0.1) {
            row(
                e,
                Scatter.sigmael(e * 1000) / barn,
                Scatter.sigmaion(e * 1000) / barn,
                Scatter.sigmaexc(e * 1000) / barn
            )
        }
        hints {
            borderStyle = Table.BorderStyle.SINGLE_LINE // or NONE
        }
    }.render(StringBuilder())

    println(table)
}