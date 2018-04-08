/*
 * written by Sebastian Voecking <seb.voeck@uni-muenster.de>
 *
 * See scatter.h for details
 *
 * Included in this file are function from Ferenc Glueck for calculation of
 * cross sections.
 */
package inr.numass.trapping

import org.apache.commons.math3.random.RandomGenerator
import org.apache.commons.math3.util.FastMath.*
import org.apache.commons.math3.util.Pair


//var generator: RandomGenerator = JDKRandomGenerator()

var debug = false

var counter = MultiCounter("Accept-reject calls")

/**
 * Using top level object instead of class
 * @author Darksnake
 */
object Scatter {

    private var fmax = 0.0

    private val a02 = 28e-22 // Bohr radius squared
    private val clight = 137.0 // velocity of light in atomic units
    private val emass = 18780.0 // Electron mass in atomic units
    private val R = 13.6 // Ryberg energy in eV

    private fun count(counterName: String) {
        if (debug) {
            counter.increase(counterName)
        }
    }

    fun randomel(E: Double, generator: RandomGenerator): Pair<Double, Double> {
        // This subroutine generates  energy loss and polar scatt. angle according to
        // electron elastic scattering in molecular hydrogen.
        // Input:
        //    E: incident electron energy in eV.
        // Output:
        //   *Eloss: energy loss in eV
        //   *theta: change of polar angle in degrees
        val H2molmass = 69e6
        val T: Double = E / 27.2
        var c = 1.0
        val b: Double
        var G: Double
        var a: Double
        val gam: Double
        var K2: Double
        val Gmax: Double = when {
            E >= 250.0 -> 1e-19
            E < 250.0 && E >= 150.0 -> 2.5e-19
            else -> 1e-18
        }
        var i: Int = 1

        count("randomel-calls")
        gam = 1.0 + T / (clight * clight) // relativistic correction factor
        b = 2.0 / (1.0 + gam) / T
        while (i < 5000) {
            count("randomel")
            c = 1.0 + b - b * (2.0 + b) / (b + 2.0 * generator.nextDouble())
            K2 = 2.0 * T * (1.0 + gam) * abs(1.0 - c) // momentum transfer squared
            a = (4.0 + K2) * (4.0 + K2) / (gam * gam)
            G = a * Del(E, c)
            if (G > Gmax * generator.nextDouble()) {
                break
            }
            i++
        }
        return Pair(2.0 * emass / H2molmass * (1.0 - c) * E, acos(c) * 180.0 / Math.PI)
    }

    fun randomexc(E: Double, generator: RandomGenerator): Pair<Double, Double> {
        // This subroutine generates  energy loss and polar scatt. angle according to
        // electron excitation scattering in molecular hydrogen.
        // Input:
        //    E: incident electron energy in eV.
        // Output:
        //   *Eloss: energy loss in eV
        //   *theta: change of polar angle in degrees

        val Ecen = 12.6 / 27.21
        val sum = DoubleArray(1001)
        var T: Double = 20000.0 / 27.2
        var c = 0.0
        var K: Double
        var xmin: Double
        var ymin: Double
        var ymax: Double
        val x: Double
        var y: Double
        var fy: Double
        val dy: Double
        var pmax: Double
        var D: Double
        val Dmax: Double
        var i: Int
        var j: Int
        var n = 0
        var N: Int
        var v = 0
        // Energy values of the excited electronic states:
        //  (from Mol. Phys. 41 (1980) 1501, in Hartree atomic units)
        val En = doubleArrayOf(12.73 / 27.2, 13.2 / 27.2, 14.77 / 27.2, 15.3 / 27.2, 14.93 / 27.2, 15.4 / 27.2, 13.06 / 27.2)
        // Probability numbers of the electronic states:
        //  (from testelectron7.c calculation )
        val p = doubleArrayOf(35.86, 40.05, 6.58, 2.26, 9.61, 4.08, 1.54)
        // Energy values of the B vibrational states:
        //   (from: Phys. Rev. A51 (1995) 3745 , in Hartree atomic units)
        val EB = doubleArrayOf(0.411, 0.417, 0.423, 0.428, 0.434, 0.439, 0.444, 0.449, 0.454, 0.459, 0.464, 0.468, 0.473, 0.477, 0.481, 0.485, 0.489, 0.493, 0.496, 0.500, 0.503, 0.507, 0.510, 0.513, 0.516, 0.519, 0.521, 0.524)
        // Energy values of the C vibrational states:
        //   (from: Phys. Rev. A51 (1995) 3745 , in Hartree atomic units)
        val EC = doubleArrayOf(0.452, 0.462, 0.472, 0.481, 0.490, 0.498, 0.506, 0.513, 0.519, 0.525, 0.530, 0.534, 0.537, 0.539)
        // Franck-Condon factors of the B vibrational states:
        //   (from: Phys. Rev. A51 (1995) 3745 )
        val pB = doubleArrayOf(4.2e-3, 1.5e-2, 3.0e-2, 4.7e-2, 6.3e-2, 7.3e-2, 7.9e-2, 8.0e-2, 7.8e-2, 7.3e-2, 6.6e-2, 5.8e-2, 5.1e-2, 4.4e-2, 3.7e-2, 3.1e-2, 2.6e-2, 2.2e-2, 1.8e-2, 1.5e-2, 1.3e-2, 1.1e-2, 8.9e-3, 7.4e-3, 6.2e-3, 5.2e-3, 4.3e-3, 3.6e-3)
        // Franck-Condon factors of the C vibrational states:
        //   (from: Phys. Rev. A51 (1995) 3745 )
        val pC = doubleArrayOf(1.2e-1, 1.9e-1, 1.9e-1, 1.5e-1, 1.1e-1, 7.5e-2, 5.0e-2, 3.3e-2, 2.2e-2, 1.4e-2, 9.3e-3, 6.0e-3, 3.7e-3, 1.8e-3)
        count("randomexc-calls")
        //
        xmin = Ecen * Ecen / (2.0 * T)
        ymin = log(xmin)
        ymax = log(8.0 * T + xmin)
        dy = (ymax - ymin) / 1000.0

        // Initialization of the sum[] vector, and fmax calculation:
        synchronized(this) {
            if (fmax == 0.0) {
                i = 0
                while (i <= 1000) {
                    y = ymin + dy * i
                    K = exp(y / 2.0)
                    sum[i] = sumexc(K)
                    if (sum[i] > fmax) {
                        fmax = sum[i]
                    }
                    i++
                }
                fmax *= 1.05
            }
        }

        //
        //  Scattering angle *theta generation:
        //
        T = E / 27.2
        val theta: Double
        if (E >= 100.0) {
            xmin = Ecen * Ecen / (2.0 * T)
            ymin = log(xmin)
            ymax = log(8.0 * T + xmin)
            //            dy = (ymax - ymin) / 1000.;
            // Generation of y values with the Neumann acceptance-rejection method:
            y = ymin
            j = 1
            while (j < 5000) {
                count("randomexc1")
                y = ymin + (ymax - ymin) * generator.nextDouble()
                K = exp(y / 2.0)
                fy = sumexc(K)
                if (fmax * generator.nextDouble() < fy) {
                    break
                }
                j++
            }
            // Calculation of c=cos(theta) and theta:
            x = exp(y)
            c = 1.0 - (x - xmin) / (4.0 * T)
            theta = acos(c) * 180.0 / Math.PI
        } else {
            Dmax = when {
                E <= 25.0 -> 60.0
                E > 25.0 && E <= 35.0 -> 95.0
                E > 35.0 && E <= 50.0 -> 150.0
                else -> 400.0
            }
            j = 1

            while (j < 5000) {
                count("randomexc2")
                c = -1.0 + 2.0 * generator.nextDouble()
                D = Dexc(E, c) * 1e22
                if (Dmax * generator.nextDouble() < D) {
                    break
                }
                j++
            }
            theta = acos(c) * 180.0 / Math.PI
        }
        // Energy loss *Eloss generation:

        // First we generate the electronic state, using the Neumann
        // acceptance-rejection method for discrete distribution:
        N = 7 // the number of electronic states in our calculation
        pmax = p[1] // the maximum of the p[] values
        j = 1
        while (j < 5000) {
            count("randomexc3")
            n = (N * generator.nextDouble()).toInt()
            if (generator.nextDouble() * pmax < p[n]) {
                break
            }
            j++
        }
        if (n < 0) {
            n = 0
        }
        if (n > 6) {
            n = 6
        }
        if (n > 1) {

        }
        val Eloss: Double
        when (n) {
            0 -> {
                // B state; we generate now a vibrational state,
                // using the Frank-Condon factors
                N = 28 // the number of B vibrational states in our calculation
                pmax = pB[7] // maximum of the pB[] values
                j = 1
                while (j < 5000) {
                    count("randomexc4")
                    v = (N * generator.nextDouble()).toInt()
                    if (generator.nextDouble() * pmax < pB[v]) {
                        break
                    }
                    j++
                }
                if (v < 0) {
                    v = 0
                }
                if (v > 27) {
                    v = 27
                }
                Eloss = EB[v] * 27.2
            }
            1 -> {
                // C state; we generate now a vibrational state,
                // using the Franck-Condon factors
                N = 14 // the number of C vibrational states in our calculation
                pmax = pC[1] // maximum of the pC[] values
                j = 1
                while (j < 5000) {
                    count("randomexc4")
                    v = (N * generator.nextDouble()).toInt()
                    if (generator.nextDouble() * pmax < pC[v]) {
                        break
                    }
                    j++
                }
                if (v < 0) {
                    v = 0
                }
                if (v > 13) {
                    v = 13
                }
                Eloss = EC[v] * 27.2
            }
            else ->
                // Bp, Bpp, D, Dp, EF states
                Eloss = En[n] * 27.2
        }
        return Pair(Eloss, theta)
    }

    fun randomion(E: Double, generator: RandomGenerator): Pair<Double, Double> {
        // This subroutine generates  energy loss and polar scatt. angle according to
        // electron ionization scattering in molecular hydrogen.
        // Input:
        //    E: incident electron energy in eV.
        // Output:
        //   *Eloss: energy loss in eV
        //   *theta: change of polar angle in degrees
        // The kinetic energy of the secondary electron is: Eloss-15.4 eV
        //
        val Ei = 15.45 / 27.21
        var c: Double
        val b: Double
        var K: Double
        val xmin: Double
        val ymin: Double
        val ymax: Double
        val x: Double
        var y: Double
        val T: Double = E / 27.2
        var G: Double
        var Gmax: Double
        var q: Double
        val h: Double
        var F: Double
        val Fmin: Double
        val Fmax: Double
        var Gp: Double
        val Elmin: Double
        val Elmax: Double = (E + 15.45) / 2.0 / 27.2
        val qmin: Double
        val qmax: Double
        var El: Double
        val wmax: Double
        var WcE: Double
        var Jstarq: Double
        var WcstarE: Double
        var w: Double
        var D2ion: Double
        var j: Int = 1
        var K2: Double
        var KK: Double
        var fE: Double
        var kej: Double
        var ki: Double
        var kf: Double
        var Rex: Double
        var arg: Double
        var arctg: Double
        var i: Int
        var st1: Double
        var st2: Double
        count("randomion-calls")
        //
        // I. Generation of theta
        // -----------------------
        Gmax = 1e-20
        if (E < 200.0) {
            Gmax = 2e-20
        }
        xmin = Ei * Ei / (2.0 * T)
        b = xmin / (4.0 * T)
        ymin = log(xmin)
        ymax = log(8.0 * T + xmin)
        // Generation of y values with the Neumann acceptance-rejection method:
        y = ymin
        while (j < 5000) {
            count("randomion1")
            y = ymin + (ymax - ymin) * generator.nextDouble()
            K = exp(y / 2.0)
            c = 1.0 + b - K * K / (4.0 * T)
            G = K * K * (Dinel(E, c) - Dexc(E, c))
            if (Gmax * generator.nextDouble() < G) {
                break
            }
            j++
        }
        // y --> x --> c --> theta
        x = exp(y)
        c = 1.0 - (x - xmin) / (4.0 * T)
        val theta = acos(c) * 180.0 / Math.PI
        //
        // II. Generation of Eloss, for fixed theta
        // ----------------------------------------
        //
        // For E<=100 eV we use subr. gensecelen
        //   (in this case no correlation between theta and Eloss)
        if (E <= 100.0) {
            return Pair(15.45 + gensecelen(E, generator), theta)
        }
        // For theta>=20 the free electron model is used
        //   (with full correlation between theta and Eloss)
        if (theta >= 20.0) {
            return Pair(E * (1.0 - c * c), theta)
        }
        // For E>100 eV and theta<20: analytical first Born approximation
        //   formula of Bethe for H atom (with modification for H2)
        //
        // Calc. of wmax:
        if (theta >= 0.7) {
            wmax = 1.1
        } else if (theta <= 0.7 && theta > 0.2) {
            wmax = 2.0
        } else if (theta <= 0.2 && theta > 0.05) {
            wmax = 4.0
        } else {
            wmax = 8.0
        }
        // We generate the q value according to the Jstarq pdf. We have to
        // define the qmin and qmax limits for this generation:
        K = sqrt(4.0 * T * (1.0 - Ei / (2.0 * T) - sqrt(1.0 - Ei / T) * c))
        Elmin = Ei
        qmin = Elmin / K - K / 2.0
        qmax = Elmax / K - K / 2.0
        //
        q = qmax
        Fmax = 1.0 / 2.0 + 1.0 / Math.PI * (q / (1.0 + q * q) + atan(q))
        q = qmin
        Fmin = 1.0 / 2.0 + 1.0 / Math.PI * (q / (1.0 + q * q) + atan(q))
        h = Fmax - Fmin
        // Generation of Eloss with the Neumann acceptance-rejection method:
        El = 0.0
        j = 1
        while (j < 5000) {
            // Generation of q with inverse transform method
            // (we use the Newton-Raphson method in order to solve the nonlinear eq.
            // for the inversion) :
            count("randomion2")
            F = Fmin + h * generator.nextDouble()
            y = 0.0
            i = 1
            while (i <= 30) {
                G = 1.0 / 2.0 + (y + sin(2.0 * y) / 2.0) / Math.PI
                Gp = (1.0 + cos(2.0 * y)) / Math.PI
                y = y - (G - F) / Gp
                if (abs(G - F) < 1e-8) {
                    break
                }
                i++
            }
            q = tan(y)
            // We have the q value, so we can define El, and calculate the weight:
            El = q * K + K * K / 2.0
            // First Born approximation formula of Bethe for e-H ionization:
            KK = K
            ki = sqrt(2.0 * T)
            kf = sqrt(2.0 * (T - El))
            K2 = 4.0 * T * (1.0 - El / (2.0 * T) - sqrt(1.0 - El / T) * c)
            if (K2 < 1e-9) {
                K2 = 1e-9
            }
            K = sqrt(K2) // momentum transfer
            Rex = 1.0 - K * K / (kf * kf) + K2 * K2 / (kf * kf * kf * kf)
            kej = sqrt(2.0 * abs(El - Ei) + 1e-8)
            st1 = K2 - 2.0 * El + 2.0
            if (abs(st1) < 1e-9) {
                st1 = 1e-9
            }
            arg = 2.0 * kej / st1
            arctg = if (arg >= 0.0) {
                atan(arg)
            } else {
                atan(arg) + Math.PI
            }
            st1 = (K + kej) * (K + kej) + 1.0
            st2 = (K - kej) * (K - kej) + 1.0
            fE = 1024.0 * El * (K2 + 2.0 / 3.0 * El) / (st1 * st1 * st1 * st2 * st2 * st2) * exp(-2.0 / kej * arctg) / (1.0 - exp(-2.0 * Math.PI / kej))
            D2ion = 2.0 * kf / ki * Rex / (El * K2) * fE
            K = KK
            //
            WcE = D2ion
            Jstarq = 16.0 / (3.0 * Math.PI * (1.0 + q * q) * (1.0 + q * q))
            WcstarE = 4.0 / (K * K * K * K * K) * Jstarq
            w = WcE / WcstarE
            if (wmax * generator.nextDouble() < w) {
                break
            }
            j++
        }

        return Pair(El * 27.2, theta)
    }

    fun gensecelen(E: Double, generator: RandomGenerator): Double {
        // This subroutine generates secondary electron energy W
        // from ionization of incident electron energy E, by using
        // the Lorentzian of Aseev  et al. (Eq. 8).
        // E and W in eV.
        val Ei = 15.45
        val eps2 = 14.3
        val b = 6.25
        val B: Double = atan((Ei - eps2) / b)
        val D: Double
        val A: Double
        val eps: Double
        val a: Double
        val u: Double = generator.nextDouble()
        val epsmax: Double

        epsmax = (E + Ei) / 2.0
        A = atan((epsmax - eps2) / b)
        D = b / (A - B)
        a = b / D * (u + D / b * B)
        eps = eps2 + b * tan(a)
        return eps - Ei
    }

    fun Del(E: Double, c: Double): Double {
        // This subroutine computes the differential cross section
        // Del= d sigma/d Omega  of  elastic electron scattering
        // on molecular hydrogen.
        // See: Nishimura et al., J. Phys. Soc. Jpn. 54 (1985) 1757.
        // Input:  E= electron kinetic energy in eV
        //         c= cos(theta), where theta is the polar scatt. angle
        // Del: in m^2/steradian
        val Cel = doubleArrayOf(-0.512, -0.512, -0.509, -0.505, -0.499, -0.491, -0.476, -0.473, -0.462, -0.452, -0.438, -0.422, -0.406, -0.388, -0.370, -0.352, -0.333, -0.314, -0.296, -0.277, -0.258, -0.239, -0.221, -0.202, -0.185, -0.167, -0.151, -0.135, -0.120, -0.105, -0.092, -0.070, -0.053, -0.039, -0.030, -0.024, -0.019, -0.016, -0.014, -0.013, -0.012, -0.009, -0.008, -0.006, -0.005, -0.004, -0.003, -0.002, -0.002, -0.001)
        val e = doubleArrayOf(0.0, 3.0, 6.0, 12.0, 20.0, 32.0, 55.0, 85.0, 150.0, 250.0)
        val t = doubleArrayOf(0.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 140.0, 180.0)
        val D = arrayOf(doubleArrayOf(2.9, 2.7, 2.5, 2.1, 1.8, 1.2, 0.9, 1.0, 1.6, 1.9), doubleArrayOf(4.2, 3.6, 3.1, 2.5, 1.9, 1.1, 0.8, 0.9, 1.3, 1.4), doubleArrayOf(6.0, 4.4, 3.2, 2.3, 1.8, 1.1, 0.7, 0.54, 0.5, 0.6), doubleArrayOf(6.0, 4.1, 2.8, 1.9, 1.3, 0.6, 0.3, 0.17, 0.16, 0.23), doubleArrayOf(4.9, 3.2, 2.0, 1.2, 0.8, 0.3, 0.15, 0.09, 0.05, 0.05), doubleArrayOf(5.2, 2.5, 1.2, 0.64, 0.36, 0.13, 0.05, 0.03, 0.016, 0.02), doubleArrayOf(4.0, 1.7, 0.7, 0.3, 0.16, 0.05, 0.02, 0.013, 0.01, 0.01), doubleArrayOf(2.8, 1.1, 0.4, 0.15, 0.07, 0.02, 0.01, 0.007, 0.004, 0.003), doubleArrayOf(1.2, 0.53, 0.2, 0.08, 0.03, 0.0074, 0.003, 0.0016, 0.001, 0.0008))
        val T: Double = E / 27.2
        var K2: Double
        var K: Double
        val d: Double
        val st1: Double
        val st2: Double
        val DH: Double
        val gam: Double
        var Delreturn = 0.0
        val CelK: Double
        val Ki: Double
        val theta: Double
        var i: Int
        var j: Int
        if (E >= 250.0) {
            gam = 1.0 + T / (clight * clight) // relativistic correction factor
            K2 = 2.0 * T * (1.0 + gam) * (1.0 - c)
            if (K2 < 0.0) {
                K2 = 1e-30
            }
            K = sqrt(K2)
            if (K < 1e-9) {
                K = 1e-9 // momentum transfer
            }
            d = 1.4009 // distance of protons in H2
            st1 = 8.0 + K2
            st2 = 4.0 + K2
            // DH is the diff. cross section for elastic electron scatt.
            // on atomic hydrogen within the first Born approximation :
            DH = 4.0 * st1 * st1 / (st2 * st2 * st2 * st2) * a02
            // CelK calculation with linear interpolation.
            // CelK is the correction of the elastic electron
            // scatt. on molecular hydrogen compared to the independent atom
            // model.
            when {
                K < 3.0 -> {
                    i = (K / 0.1).toInt()
                    Ki = i * 0.1
                    CelK = Cel[i] + (K - Ki) / 0.1 * (Cel[i + 1] - Cel[i])
                }
                K >= 3.0 && K < 5.0 -> {
                    i = (30 + (K - 3.0) / 0.2).toInt()
                    Ki = 3.0 + (i - 30) * 0.2
                    CelK = Cel[i] + (K - Ki) / 0.2 * (Cel[i + 1] - Cel[i])
                }
                K >= 5.0 && K < 9.49 -> {
                    i = (40 + (K - 5.0) / 0.5).toInt()
                    Ki = 5.0 + (i - 40) * 0.5
                    CelK = Cel[i] + (K - Ki) / 0.5 * (Cel[i + 1] - Cel[i])
                }
                else -> CelK = 0.0
            }
            Delreturn = 2.0 * gam * gam * DH * (1.0 + sin(K * d) / (K * d)) * (1.0 + CelK)
        } else {
            theta = acos(c) * 180.0 / Math.PI
            i = 0
            while (i <= 8) {
                if (E >= e[i] && E < e[i + 1]) {
                    j = 0
                    while (j <= 8) {
                        if (theta >= t[j] && theta < t[j + 1]) {
                            Delreturn = 1e-20 * (D[i][j] + (D[i][j + 1] - D[i][j]) * (theta - t[j]) / (t[j + 1] - t[j]))
                        }
                        j++
                    }
                }
                i++
            }
        }
        return Delreturn
    }

    fun Dexc(E: Double, c: Double): Double {
        // This subroutine computes the differential cross section
        // Del= d sigma/d Omega  of excitation electron scattering
        // on molecular hydrogen.
        // Input:  E= electron kinetic energy in eV
        //         c= cos(theta), where theta is the polar scatt. angle
        // Dexc: in m^2/steradian
        var K2: Double
        val K: Double
        var sigma = 0.0
        val T: Double = E / 27.2
        val theta: Double
        val EE = 12.6 / 27.2
        val e = doubleArrayOf(0.0, 25.0, 35.0, 50.0, 100.0)
        val t = doubleArrayOf(0.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 180.0)
        val D = arrayOf(doubleArrayOf(60.0, 43.0, 27.0, 18.0, 13.0, 8.0, 6.0, 6.0, 6.0), doubleArrayOf(95.0, 70.0, 21.0, 9.0, 6.0, 3.0, 2.0, 2.0, 2.0), doubleArrayOf(150.0, 120.0, 32.0, 8.0, 3.7, 1.9, 1.2, 0.8, 0.8), doubleArrayOf(400.0, 200.0, 12.0, 2.0, 1.4, 0.7, 0.3, 0.2, 0.2))
        var i: Int
        var j: Int
        //
        if (E >= 100.0) {
            K2 = 4.0 * T * (1.0 - EE / (2.0 * T) - sqrt(1.0 - EE / T) * c)
            if (K2 < 1e-9) {
                K2 = 1e-9
            }
            K = sqrt(K2) // momentum transfer
            sigma = 2.0 / K2 * sumexc(K) * a02
        } else if (E <= 10.0) {
            sigma = 0.0
        } else {
            theta = acos(c) * 180.0 / Math.PI
            i = 0
            while (i <= 3) {
                if (E >= e[i] && E < e[i + 1]) {
                    j = 0
                    while (j <= 7) {
                        if (theta >= t[j] && theta < t[j + 1]) {
                            sigma = 1e-22 * (D[i][j] + (D[i][j + 1] - D[i][j]) * (theta - t[j]) / (t[j + 1] - t[j]))
                        }
                        j++
                    }
                }
                i++
            }
        }
        return sigma
    }

    fun Dinel(E: Double, c: Double): Double {
        // This subroutine computes the differential cross section
        // Dinel= d sigma/d Omega  of  inelastic electron scattering
        // on molecular hydrogen, within the first Born approximation.
        // Input:  E= electron kinetic energy in eV
        //         c= cos(theta), where theta is the polar scatt. angle
        // Dinel: in m2/steradian
        val Cinel = doubleArrayOf(-0.246, -0.244, -0.239, -0.234, -0.227, -0.219, -0.211, -0.201, -0.190, -0.179, -0.167, -0.155, -0.142, -0.130, -0.118, -0.107, -0.096, -0.085, -0.076, -0.067, -0.059, -0.051, -0.045, -0.039, -0.034, -0.029, -0.025, -0.022, -0.019, -0.016, -0.014, -0.010, -0.008, -0.006, -0.004, -0.003, -0.003, -0.002, -0.002, -0.001, -0.001, -0.001, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000)
        val Ei = 0.568
        val T: Double = E / 27.2
        var K2: Double
        val K: Double
        val st1: Double
        val F: Double
        val DH: Double
        val Dinelreturn: Double
        val CinelK: Double
        val Ki: Double
        val i: Int
        if (E < 16.0) {
            return Dexc(E, c)
        }
        K2 = 4.0 * T * (1.0 - Ei / (2.0 * T) - sqrt(1.0 - Ei / T) * c)
        if (K2 < 1e-9) {
            K2 = 1e-9
        }
        K = sqrt(K2) // momentum transfer
        st1 = 1.0 + K2 / 4.0
        F = 1.0 / (st1 * st1) // scatt. formfactor of hydrogen atom
        // DH is the diff. cross section for inelastic electron scatt.
        // on atomic hydrogen within the first Born approximation :
        DH = 4.0 / (K2 * K2) * (1.0 - F * F) * a02
        // CinelK calculation with linear interpolation.
        // CinelK is the correction of the inelastic electron
        // scatt. on molecular hydrogen compared to the independent atom
        // model.
        if (K < 3.0) {
            i = (K / 0.1).toInt()
            Ki = i * 0.1
            CinelK = Cinel[i] + (K - Ki) / 0.1 * (Cinel[i + 1] - Cinel[i])
        } else if (K >= 3.0 && K < 5.0) {
            i = (30 + (K - 3.0) / 0.2).toInt()
            Ki = 3.0 + (i - 30) * 0.2
            CinelK = Cinel[i] + (K - Ki) / 0.2 * (Cinel[i + 1] - Cinel[i])
        } else if (K >= 5.0 && K < 9.49) {
            i = (40 + (K - 5.0) / 0.5).toInt()
            Ki = 5.0 + (i - 40) * 0.5
            CinelK = Cinel[i] + (K - Ki) / 0.5 * (Cinel[i + 1] - Cinel[i])
        } else {
            CinelK = 0.0
        }
        Dinelreturn = 2.0 * DH * (1.0 + CinelK)
        return Dinelreturn
    }

    private fun sumexc(K: Double): Double {
        val Kvec = doubleArrayOf(0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0)
        val fvec = arrayOf(doubleArrayOf(2.907e-1, 2.845e-1, 2.665e-1, 2.072e-1, 1.389e-1, // B
                8.238e-2, 4.454e-2, 2.269e-2, 7.789e-3, 2.619e-3, 1.273e-3, 2.218e-4, 4.372e-5, 2.889e-6, 4.247e-7), doubleArrayOf(3.492e-1, 3.367e-1, 3.124e-1, 2.351e-1, 1.507e-1, // C
                8.406e-2, 4.214e-2, 1.966e-2, 5.799e-3, 1.632e-3, 6.929e-4, 8.082e-5, 9.574e-6, 1.526e-7, 7.058e-9), doubleArrayOf(6.112e-2, 5.945e-2, 5.830e-2, 5.072e-2, 3.821e-2, // Bp
                2.579e-2, 1.567e-2, 8.737e-3, 3.305e-3, 1.191e-3, 6.011e-4, 1.132e-4, 2.362e-5, 1.603e-6, 2.215e-7), doubleArrayOf(2.066e-2, 2.127e-2, 2.137e-2, 1.928e-2, 1.552e-2, // Bpp
                1.108e-2, 7.058e-3, 4.069e-3, 1.590e-3, 5.900e-4, 3.046e-4, 6.142e-5, 1.369e-5, 9.650e-7, 1.244e-7), doubleArrayOf(9.405e-2, 9.049e-2, 8.613e-2, 7.301e-2, 5.144e-2, // D
                3.201e-2, 1.775e-2, 8.952e-3, 2.855e-3, 8.429e-4, 3.655e-4, 4.389e-5, 5.252e-6, 9.010e-8, 7.130e-9), doubleArrayOf(4.273e-2, 3.862e-2, 3.985e-2, 3.362e-2, 2.486e-2, // Dp
                1.612e-2, 9.309e-3, 4.856e-3, 1.602e-3, 4.811e-4, 2.096e-4, 2.498e-5, 2.905e-6, 5.077e-8, 6.583e-9), doubleArrayOf(0.000e-3, 2.042e-3, 7.439e-3, 2.200e-2, 3.164e-2, // EF
                3.161e-2, 2.486e-2, 1.664e-2, 7.562e-3, 3.044e-3, 1.608e-3, 3.225e-4, 7.120e-5, 6.290e-6, 1.066e-6))
        val EeV = doubleArrayOf(12.73, 13.20, 14.77, 15.3, 14.93, 15.4, 13.06)
        var n: Int = 0
        var j: Int
        var jmin = 0
        val nmax: Int = 6
        var En: Double
        var sum: Double = 0.0
        val f = DoubleArray(7)
        val x4 = DoubleArray(4)
        val f4 = DoubleArray(4)

        //
        while (n <= nmax) {
            En = EeV[n] / 27.21 // En is the excitation energy in Hartree atomic units
            if (K >= 5.0) {
                f[n] = 0.0
            } else if (K >= 3.0 && K <= 4.0) {
                f[n] = fvec[n][12] + (K - 3.0) * (fvec[n][13] - fvec[n][12])
            } else if (K >= 4.0 && K <= 5.0) {
                f[n] = fvec[n][13] + (K - 4.0) * (fvec[n][14] - fvec[n][13])
            } else {
                j = 0
                while (j < 14) {
                    if (K >= Kvec[j] && K <= Kvec[j + 1]) {
                        jmin = j - 1
                    }
                    j++
                }
                if (jmin < 0) {
                    jmin = 0
                }
                if (jmin > 11) {
                    jmin = 11
                }
                j = 0
                while (j <= 3) {
                    x4[j] = Kvec[jmin + j]
                    f4[j] = fvec[n][jmin + j]
                    j++
                }
                f[n] = lagrange(4, x4, f4, K)
            }
            sum += f[n] / En
            n++
        }
        return sum
    }

    fun lagrange(n: Int, xn: DoubleArray, fn: DoubleArray, x: Double): Double {
        var i: Int
        var j: Int
        var f: Double = 0.0
        var aa: Double
        var bb: Double
        val a = DoubleArray(100)
        val b = DoubleArray(100)
        j = 0
        while (j < n) {
            i = 0
            while (i < n) {
                a[i] = x - xn[i]
                b[i] = xn[j] - xn[i]
                i++
            }
            bb = 1.0
            aa = bb
            b[j] = aa
            a[j] = b[j]
            i = 0
            while (i < n) {
                aa *= a[i]
                bb *= b[i]
                i++
            }
            f += fn[j] * aa / bb
            j++
        }
        return f
    }

    fun sigmael(E: Double): Double {
        // This function computes the total elastic cross section of
        // electron scatt. on molecular hydrogen.
        // See: Liu, Phys. Rev. A35 (1987) 591,
        //      Trajmar, Phys Reports 97 (1983) 221.
        // E: incident electron energy in eV
        // sigmael: cross section in m^2
        val e = doubleArrayOf(0.0, 1.5, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 60.0, 100.0, 150.0, 200.0, 300.0, 400.0)
        val s = doubleArrayOf(9.6, 13.0, 15.0, 12.0, 10.0, 7.0, 5.6, 3.3, 1.1, 0.9, 0.5, 0.36, 0.23, 0.15)
        val gam: Double
        var sigma = 0.0
        val T: Double = E / 27.2
        var i: Int
        if (E >= 400.0) {
            gam = (emass + T) / emass
            sigma = gam * gam * Math.PI / (2.0 * T) * (4.2106 - 1.0 / T) * a02
        } else {
            i = 0
            while (i <= 12) {
                if (E >= e[i] && E < e[i + 1]) {
                    sigma = 1e-20 * (s[i] + (s[i + 1] - s[i]) * (E - e[i]) / (e[i + 1] - e[i]))
                }
                i++
            }
        }
        return sigma
    }

    fun sigmaexc(E: Double): Double {
        // This function computes the electronic excitation cross section of
        // electron scatt. on molecular hydrogen.
        // E: incident electron energy in eV,
        // sigmaexc: cross section in m^2
        val sigma: Double
        if (E < 9.8) {
            sigma = 1e-40
        } else if (E in 9.8..250.0) {
            sigma = sigmaBC(E) + sigmadiss10(E) + sigmadiss15(E)
        } else {
            sigma = 4.0 * Math.PI * a02 * R / E * (0.80 * log(E / R) + 0.28)
        }
        //    sigma=sigmainel(E)-sigmaion(E);
        return sigma
    }

    fun sigmaion(E: Double): Double {
        // This function computes the total ionization cross section of
        // electron scatt. on molecular hydrogen.
        // E: incident electron energy in eV,
        // sigmaion: total ionization cross section of
        //   e+H2 --> e+e+H2^+  or  e+e+H^+ +H
        // process in m^2.
        //
        // E<250 eV: Eq. 5 of J. Chem. Phys. 104 (1996) 2956
        // E>250: sigma_i formula on page 107 in
        //   Phys. Rev. A7 (1973) 103.
        // Good agreement with measured results of
        // PR A 54 (1996) 2146, and
        //   Physica 31 (1965) 94.
        //
        val B = 15.43
        val U = 15.98
        val sigma: Double
        val t: Double
        val u: Double
        val S: Double
        val r: Double
        val lnt: Double
        if (E < 16.0) {
            sigma = 1e-40
        } else if (E >= 16.0 && E <= 250.0) {
            t = E / B
            u = U / B
            r = R / B
            S = 4.0 * Math.PI * a02 * 2.0 * r * r
            lnt = log(t)
            sigma = S / (t + u + 1.0) * (lnt / 2.0 * (1.0 - 1.0 / (t * t)) + 1.0 - 1.0 / t - lnt / (t + 1.0))
        } else {
            sigma = 4.0 * Math.PI * a02 * R / E * (0.82 * log(E / R) + 1.3)
        }
        return sigma
    }

    fun sigmaBC(E: Double): Double {
        // This function computes the sigmaexc electronic excitation
        // cross section to the B and C states, with energy loss
        // about 12.5 eV.
        // E is incident electron energy in eV,
        // sigmaexc in m^2
        val aB = doubleArrayOf(-4.2935194e2, 5.1122109e2, -2.8481279e2, 8.8310338e1, -1.6659591e1, 1.9579609, -1.4012824e-1, 5.5911348e-3, -9.5370103e-5)
        val aC = doubleArrayOf(-8.1942684e2, 9.8705099e2, -5.3095543e2, 1.5917023e2, -2.9121036e1, 3.3321027, -2.3305961e-1, 9.1191781e-3, -1.5298950e-4)
        var lnsigma: Double = 0.0
        var lnE: Double = log(E)
        var lnEn: Double = 1.0
        val sigmaB: Double
        var Emin: Double = 12.5
        var sigma: Double = 0.0
        val sigmaC: Double
        var n: Int
        if (E < Emin) {
            sigmaB = 0.0
        } else {
            n = 0
            while (n <= 8) {
                lnsigma += aB[n] * lnEn
                lnEn *= lnE
                n++
            }
            sigmaB = exp(lnsigma)
        }
        sigma += sigmaB
        //  sigma=0.;
        // C state:
        Emin = 15.8
        lnE = log(E)
        lnEn = 1.0
        lnsigma = 0.0
        if (E < Emin) {
            sigmaC = 0.0
        } else {
            n = 0
            while (n <= 8) {
                lnsigma += aC[n] * lnEn
                lnEn *= lnE
                n++
            }
            sigmaC = exp(lnsigma)
        }
        sigma += sigmaC
        return sigma * 1e-4
    }

    //////////////////////////////////////////////////////////////////
    fun sigmadiss10(E: Double): Double {
        // This function computes the sigmadiss10 electronic
        // dissociative excitation
        // cross section, with energy loss
        // about 10 eV.
        // E is incident electron energy in eV,
        // sigmadiss10 in m^2
        val a = doubleArrayOf(-2.297914361e5, 5.303988579e5, -5.316636672e5, 3.022690779e5, -1.066224144e5, 2.389841369e4, -3.324526406e3, 2.624761592e2, -9.006246604)
        var lnsigma: Double = 0.0
        val lnE: Double = log(E)
        var lnEn: Double = 1.0
        val Emin = 9.8
        var n: Int
        //  E is in eV
        val sigma = if (E < Emin) {
            0.0
        } else {
            n = 0
            while (n <= 8) {
                lnsigma += a[n] * lnEn
                lnEn *= lnE
                n++
            }
            exp(lnsigma)
        }
        return sigma * 1e-4
    }

    //////////////////////////////////////////////////////////////////
    fun sigmadiss15(E: Double): Double {
        // This function computes the sigmadiss15 electronic
        // dissociative excitation
        // cross section, with energy loss
        // about 15 eV.
        // E is incident electron energy in eV,
        // sigmadiss15 in m^2
        val a = doubleArrayOf(-1.157041752e3, 1.501936271e3, -8.6119387e2, 2.754926257e2, -5.380465012e1, 6.573972423, -4.912318139e-1, 2.054926773e-2, -3.689035889e-4)
        var lnsigma: Double = 0.0
        val lnE: Double = log(E)
        var lnEn: Double
        val Emin: Double = 16.5
        var n: Int
        //  E is in eV
        lnEn = 1.0
        val sigma = if (E < Emin) {
            0.0
        } else {
            n = 0
            while (n <= 8) {
                lnsigma += a[n] * lnEn
                lnEn *= lnE
                n++
            }
            exp(lnsigma)
        }
        return sigma * 1e-4
    }

    /**
     * Полное сечение с учетом квазиупругих столкновений
     *
     * @param E
     * @return
     */
    fun sigmaTotal(E: Double): Double {
        return sigmael(E) + sigmaexc(E) + sigmaion(E)

    }
}
