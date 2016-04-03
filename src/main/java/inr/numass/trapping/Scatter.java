/*
 * written by Sebastian Voecking <seb.voeck@uni-muenster.de>
 *
 * See scatter.h for details
 *
 * Included in this file are function from Ferenc Glueck for calculation of
 * cross sections.
 */
package inr.numass.trapping;

import static java.lang.Math.abs;
import static java.lang.Math.acos;
import static java.lang.Math.atan;
import static java.lang.Math.cos;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.tan;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.Pair;

/**
 *
 * @author Darksnake
 */
public class Scatter {

    static final double a02 = 28e-22; // Bohr radius squared
    static final double clight = 137; // velocity of light in atomic units    
    static final double emass = 18780; // Electron mass in atomic units
    static final double R = 13.6; // Ryberg energy in eV
    private final RandomGenerator generator;
    MultiCounter counter = new MultiCounter("Accept-reject calls");

    private double fmax = 0;

    public Scatter(RandomGenerator generator) {
        this.generator = generator;
    }

    public Pair<Double, Double> randomel(double E) {
        // This subroutine generates  energy loss and polar scatt. angle according to
        // electron elastic scattering in molecular hydrogen.
        // Input:
        //    E: incident electron energy in eV.
        // Output:
        //   *Eloss: energy loss in eV
        //   *theta: change of polar angle in degrees
        double H2molmass = 69.e6;
        double T, c = 1, b, G, a, gam, K2, Gmax;
        double[] u = new double[3];
        int i;

        counter.increase("randomel-calls");
        if (E >= 250.) {
            Gmax = 1.e-19;
        } else if (E < 250. && E >= 150.) {
            Gmax = 2.5e-19;
        } else {
            Gmax = 1.e-18;
        }
        T = E / 27.2;
        gam = 1. + T / (clight * clight); // relativistic correction factor
        b = 2. / (1. + gam) / T;
        for (i = 1; i < 5000; i++) {
            counter.increase("randomel");
            c = 1. + b - b * (2. + b) / (b + 2. * generator.nextDouble());
            K2 = 2. * T * (1. + gam) * abs(1d - c); // momentum transfer squared
            a = (4. + K2) * (4. + K2) / (gam * gam);
            G = a * Del(E, c);
            if (G > Gmax * generator.nextDouble()) {
                break;
            }
        }
        return new Pair<>(2d * emass / H2molmass * (1d - c) * E, acos(c) * 180d / Math.PI);
    }

    Pair<Double, Double> randomexc(double E) {
        // This subroutine generates  energy loss and polar scatt. angle according to
        // electron excitation scattering in molecular hydrogen.
        // Input:
        //    E: incident electron energy in eV.
        // Output:
        //   *Eloss: energy loss in eV
        //   *theta: change of polar angle in degrees

        double Ecen = 12.6 / 27.21;
        double[] sum = new double[1001];
        double T, c = 0., K, xmin, ymin, ymax, x, y, fy, dy, pmax;
        double D, Dmax;
        int i, j, n = 0, N, v = 0;
        // Energy values of the excited electronic states:
        //  (from Mol. Phys. 41 (1980) 1501, in Hartree atomic units)
        double[] En = {12.73 / 27.2, 13.2 / 27.2, 14.77 / 27.2, 15.3 / 27.2,
            14.93 / 27.2, 15.4 / 27.2, 13.06 / 27.2};
        // Probability numbers of the electronic states:
        //  (from testelectron7.c calculation )
        double[] p = {35.86, 40.05, 6.58, 2.26, 9.61, 4.08, 1.54};
        // Energy values of the B vibrational states:
        //   (from: Phys. Rev. A51 (1995) 3745 , in Hartree atomic units)
        double[] EB = {0.411, 0.417, 0.423, 0.428, 0.434, 0.439, 0.444, 0.449,
            0.454, 0.459, 0.464, 0.468, 0.473, 0.477, 0.481, 0.485,
            0.489, 0.493, 0.496, 0.500, 0.503, 0.507, 0.510, 0.513,
            0.516, 0.519, 0.521, 0.524};
        // Energy values of the C vibrational states:
        //   (from: Phys. Rev. A51 (1995) 3745 , in Hartree atomic units)
        double[] EC = {0.452, 0.462, 0.472, 0.481, 0.490, 0.498, 0.506, 0.513,
            0.519, 0.525, 0.530, 0.534, 0.537, 0.539};
        // Franck-Condon factors of the B vibrational states:
        //   (from: Phys. Rev. A51 (1995) 3745 )
        double[] pB = {4.2e-3, 1.5e-2, 3.0e-2, 4.7e-2, 6.3e-2, 7.3e-2, 7.9e-2,
            8.0e-2, 7.8e-2, 7.3e-2, 6.6e-2, 5.8e-2, 5.1e-2, 4.4e-2,
            3.7e-2, 3.1e-2, 2.6e-2, 2.2e-2, 1.8e-2, 1.5e-2, 1.3e-2,
            1.1e-2, 8.9e-3, 7.4e-3, 6.2e-3, 5.2e-3, 4.3e-3, 3.6e-3};
        // Franck-Condon factors of the C vibrational states:
        //   (from: Phys. Rev. A51 (1995) 3745 )
        double[] pC = {1.2e-1, 1.9e-1, 1.9e-1, 1.5e-1, 1.1e-1, 7.5e-2, 5.0e-2,
            3.3e-2, 2.2e-2, 1.4e-2, 9.3e-3, 6.0e-3, 3.7e-3, 1.8e-3};
        counter.increase("randomexc-calls");
        T = 20000. / 27.2;
        //
        xmin = Ecen * Ecen / (2. * T);
        ymin = log(xmin);
        ymax = log(8. * T + xmin);
        dy = (ymax - ymin) / 1000.;

        // Initialization of the sum[] vector, and fmax calculation:
        if (fmax == 0) {
            synchronized (this) {
                for (i = 0; i <= 1000; i++) {
                    y = ymin + dy * i;
                    K = exp(y / 2.);
                    sum[i] = sumexc(K);
                    if (sum[i] > fmax) {
                        fmax = sum[i];
                    }
                }
                fmax = 1.05 * fmax;
            }
        }

        //
        //  Scattering angle *theta generation:
        //
        T = E / 27.2;
        double theta;
        if (E >= 100.) {
            xmin = Ecen * Ecen / (2. * T);
            ymin = log(xmin);
            ymax = log(8. * T + xmin);
//            dy = (ymax - ymin) / 1000.;
            // Generation of y values with the Neumann acceptance-rejection method:
            y = ymin;
            for (j = 1; j < 5000; j++) {
                counter.increase("randomexc1");
                y = ymin + (ymax - ymin) * generator.nextDouble();
                K = exp(y / 2.);
                fy = sumexc(K);
                if (fmax * generator.nextDouble() < fy) {
                    break;
                }
            }
            // Calculation of c=cos(theta) and theta:
            x = exp(y);
            c = 1. - (x - xmin) / (4. * T);
            theta = acos(c) * 180. / Math.PI;
        } else {
            if (E <= 25.) {
                Dmax = 60.;
            } else if (E > 25. && E <= 35.) {
                Dmax = 95.;
            } else if (E > 35. && E <= 50.) {
                Dmax = 150.;
            } else {
                Dmax = 400.;
            }
            for (j = 1; j < 5000; j++) {
                counter.increase("randomexc2");
                c = -1. + 2. * generator.nextDouble();
                D = Dexc(E, c) * 1.e22;
                if (Dmax * generator.nextDouble() < D) {
                    break;
                }
            }
            theta = acos(c) * 180. / Math.PI;
        }
        // Energy loss *Eloss generation:

        // First we generate the electronic state, using the Neumann
        // acceptance-rejection method for discrete distribution:
        N = 7; // the number of electronic states in our calculation
        pmax = p[1]; // the maximum of the p[] values
        for (j = 1; j < 5000; j++) {
            counter.increase("randomexc3");
            n = (int) (N * generator.nextDouble());
            if (generator.nextDouble() * pmax < p[n]) {
                break;
            }
        }
        if (n < 0) {
            n = 0;
        }
        if (n > 6) {
            n = 6;
        }
        if (n > 1) {

        }
        double Eloss;
        switch (n) {
            case 0:
                // B state; we generate now a vibrational state,
                // using the Frank-Condon factors
                N = 28; // the number of B vibrational states in our calculation
                pmax = pB[7]; // maximum of the pB[] values
                for (j = 1; j < 5000; j++) {
                    counter.increase("randomexc4");
                    v = (int) (N * generator.nextDouble());
                    if (generator.nextDouble() * pmax < pB[v]) {
                        break;
                    }
                }
                if (v < 0) {
                    v = 0;
                }
                if (v > 27) {
                    v = 27;
                }
                Eloss = EB[v] * 27.2;
                break;
            case 1:
                // C state; we generate now a vibrational state,
                // using the Franck-Condon factors
                N = 14; // the number of C vibrational states in our calculation
                pmax = pC[1]; // maximum of the pC[] values
                for (j = 1; j < 5000; j++) {
                    counter.increase("randomexc4");
                    v = (int) (N * generator.nextDouble());
                    if (generator.nextDouble() * pmax < pC[v]) {
                        break;
                    }
                }
                if (v < 0) {
                    v = 0;
                }
                if (v > 13) {
                    v = 13;
                }
                Eloss = EC[v] * 27.2;
                break;
            default:
                // Bp, Bpp, D, Dp, EF states
                Eloss = En[n] * 27.2;
                break;
        }
        return new Pair<>(Eloss, theta);
    }

    Pair<Double, Double> randomion(double E) {
        // This subroutine generates  energy loss and polar scatt. angle according to
        // electron ionization scattering in molecular hydrogen.
        // Input:
        //    E: incident electron energy in eV.
        // Output:
        //   *Eloss: energy loss in eV
        //   *theta: change of polar angle in degrees
        // The kinetic energy of the secondary electron is: Eloss-15.4 eV
        //
        double Ei = 15.45 / 27.21;
        double c, b, K, xmin, ymin, ymax, x, y, T, G, W, Gmax;
        double q, h, F, Fmin, Fmax, Gp, Elmin, Elmax, qmin, qmax, El, wmax;
        double WcE, Jstarq, WcstarE, w, D2ion;
        int j;
        double K2, KK, fE, kej, ki, kf, Rex, arg, arctg;
        int i;
        double st1, st2;
        counter.increase("randomion-calls");
        //
        // I. Generation of theta
        // -----------------------
        Gmax = 1.e-20;
        if (E < 200.) {
            Gmax = 2.e-20;
        }
        T = E / 27.2;
        xmin = Ei * Ei / (2. * T);
        b = xmin / (4. * T);
        ymin = log(xmin);
        ymax = log(8. * T + xmin);
        // Generation of y values with the Neumann acceptance-rejection method:
        y = ymin;
        for (j = 1; j < 5000; j++) {
            counter.increase("randomion1");
            y = ymin + (ymax - ymin) * generator.nextDouble();
            K = exp(y / 2.);
            c = 1. + b - K * K / (4. * T);
            G = K * K * (Dinel(E, c) - Dexc(E, c));
            if (Gmax * generator.nextDouble() < G) {
                break;
            }
        }
        // y --> x --> c --> theta
        x = exp(y);
        c = 1. - (x - xmin) / (4. * T);
        double theta = acos(c) * 180. / Math.PI;
        //
        // II. Generation of Eloss, for fixed theta
        // ----------------------------------------
        //
        // For E<=100 eV we use subr. gensecelen
        //   (in this case no correlation between theta and Eloss)
        if (E <= 100.) {
            return new Pair<>(15.45 + gensecelen(E), theta);
        }
        // For theta>=20 the free electron model is used
        //   (with full correlation between theta and Eloss)
        if (theta >= 20.) {
            return new Pair<>(E * (1. - c * c), theta);
        }
        // For E>100 eV and theta<20: analytical first Born approximation
        //   formula of Bethe for H atom (with modification for H2)
        //
        // Calc. of wmax:
        if (theta >= 0.7) {
            wmax = 1.1;
        } else if (theta <= 0.7 && theta > 0.2) {
            wmax = 2.;
        } else if (theta <= 0.2 && theta > 0.05) {
            wmax = 4.;
        } else {
            wmax = 8.;
        }
        // We generate the q value according to the Jstarq pdf. We have to
        // define the qmin and qmax limits for this generation:
        K = sqrt(4. * T * (1. - Ei / (2. * T) - sqrt(1. - Ei / T) * c));
        Elmin = Ei;
        Elmax = (E + 15.45) / 2. / 27.2;
        qmin = Elmin / K - K / 2.;
        qmax = Elmax / K - K / 2.;
        //
        q = qmax;
        Fmax = 1. / 2. + 1. / Math.PI * (q / (1. + q * q) + atan(q));
        q = qmin;
        Fmin = 1. / 2. + 1. / Math.PI * (q / (1. + q * q) + atan(q));
        h = Fmax - Fmin;
        // Generation of Eloss with the Neumann acceptance-rejection method:
        El = 0;
        for (j = 1; j < 5000; j++) {
            // Generation of q with inverse transform method
            // (we use the Newton-Raphson method in order to solve the nonlinear eq.
            // for the inversion) :
            counter.increase("randomion2");
            F = Fmin + h * generator.nextDouble();
            y = 0.;
            for (i = 1; i <= 30; i++) {
                G = 1. / 2. + (y + sin(2. * y) / 2.) / Math.PI;
                Gp = (1. + cos(2. * y)) / Math.PI;
                y = y - (G - F) / Gp;
                if (abs(G - F) < 1.e-8) {
                    break;
                }
            }
            q = tan(y);
            // We have the q value, so we can define El, and calculate the weight:
            El = q * K + K * K / 2.;
            // First Born approximation formula of Bethe for e-H ionization:
            KK = K;
            ki = sqrt(2. * T);
            kf = sqrt(2. * (T - El));
            K2 = 4. * T * (1. - El / (2. * T) - sqrt(1. - El / T) * c);
            if (K2 < 1.e-9) {
                K2 = 1.e-9;
            }
            K = sqrt(K2); // momentum transfer
            Rex = 1. - K * K / (kf * kf) + K2 * K2 / (kf * kf * kf * kf);
            kej = sqrt(2. * abs(El - Ei) + 1.e-8);
            st1 = K2 - 2. * El + 2.;
            if (abs(st1) < 1.e-9) {
                st1 = 1.e-9;
            }
            arg = 2. * kej / st1;
            if (arg >= 0.) {
                arctg = atan(arg);
            } else {
                arctg = atan(arg) + Math.PI;
            }
            st1 = (K + kej) * (K + kej) + 1.;
            st2 = (K - kej) * (K - kej) + 1.;
            fE = 1024. * El * (K2 + 2. / 3. * El) / (st1 * st1 * st1 * st2 * st2 * st2)
                    * exp(-2. / kej * arctg) / (1. - exp(-2. * Math.PI / kej));
            D2ion = 2. * kf / ki * Rex / (El * K2) * fE;
            K = KK;
            //
            WcE = D2ion;
            Jstarq = 16. / (3. * Math.PI * (1. + q * q) * (1. + q * q));
            WcstarE = 4. / (K * K * K * K * K) * Jstarq;
            w = WcE / WcstarE;
            if (wmax * generator.nextDouble() < w) {
                break;
            }
        }

        return new Pair<>(El * 27.2, theta);
    }

    double gensecelen(double E) {
        // This subroutine generates secondary electron energy W
        // from ionization of incident electron energy E, by using
        // the Lorentzian of Aseev  et al. (Eq. 8).
        // E and W in eV.
        double Ei = 15.45, eps2 = 14.3, b = 6.25;
        double B;
        double D, A, eps, a, u, epsmax;
        int iff = 0;
        B = 0;
        if (iff == 0) {
            B = atan((Ei - eps2) / b);
            iff = 1;
        }
        epsmax = (E + Ei) / 2.;
        A = atan((epsmax - eps2) / b);
        D = b / (A - B);
        u = generator.nextDouble();
        a = b / D * (u + D / b * B);
        eps = eps2 + b * tan(a);
        return eps - Ei;
    }

    double Del(double E, double c) {
        // This subroutine computes the differential cross section
        // Del= d sigma/d Omega  of  elastic electron scattering
        // on molecular hydrogen.
        // See: Nishimura et al., J. Phys. Soc. Jpn. 54 (1985) 1757.
        // Input:  E= electron kinetic energy in eV
        //         c= cos(theta), where theta is the polar scatt. angle
        // Del: in m^2/steradian
        double[] Cel = {
            -0.512, -0.512, -0.509, -0.505, -0.499,
            -0.491, -0.476, -0.473, -0.462, -0.452,
            -0.438, -0.422, -0.406, -0.388, -0.370,
            -0.352, -0.333, -0.314, -0.296, -0.277,
            -0.258, -0.239, -0.221, -0.202, -0.185,
            -0.167, -0.151, -0.135, -0.120, -0.105,
            -0.092, -0.070, -0.053, -0.039, -0.030,
            -0.024, -0.019, -0.016, -0.014, -0.013,
            -0.012, -0.009, -0.008, -0.006, -0.005,
            -0.004, -0.003, -0.002, -0.002, -0.001
        };
        double[] e = {0., 3., 6., 12., 20., 32., 55., 85., 150., 250.};
        double[] t = {0., 10., 20., 30., 40., 60., 80., 100., 140., 180.};
        double[][] D = {
            {2.9, 2.7, 2.5, 2.1, 1.8, 1.2, 0.9, 1., 1.6, 1.9},
            {4.2, 3.6, 3.1, 2.5, 1.9, 1.1, 0.8, 0.9, 1.3, 1.4},
            {6., 4.4, 3.2, 2.3, 1.8, 1.1, 0.7, 0.54, 0.5, 0.6},
            {6., 4.1, 2.8, 1.9, 1.3, 0.6, 0.3, 0.17, 0.16, 0.23},
            {4.9, 3.2, 2., 1.2, 0.8, 0.3, 0.15, 0.09, 0.05, 0.05},
            {5.2, 2.5, 1.2, 0.64, 0.36, 0.13, 0.05, 0.03, 0.016, 0.02},
            {4., 1.7, 0.7, 0.3, 0.16, 0.05, 0.02, 0.013, 0.01, 0.01},
            {2.8, 1.1, 0.4, 0.15, 0.07, 0.02, 0.01, 0.007, 0.004, 0.003},
            {1.2, 0.53, 0.2, 0.08, 0.03, 0.0074, 0.003, 0.0016, 0.001, 0.0008}
        };
        double T, K2, K, d, st1, st2, DH, gam, Delreturn = 0., CelK, Ki, theta;
        int i, j;
        T = E / 27.2;
        if (E >= 250.) {
            gam = 1. + T / (clight * clight); // relativistic correction factor
            K2 = 2. * T * (1. + gam) * (1. - c);
            if (K2 < 0.) {
                K2 = 1.e-30;
            }
            K = sqrt(K2);
            if (K < 1.e-9) {
                K = 1.e-9; // momentum transfer
            }
            d = 1.4009; // distance of protons in H2
            st1 = 8. + K2;
            st2 = 4. + K2;
            // DH is the diff. cross section for elastic electron scatt.
            // on atomic hydrogen within the first Born approximation :
            DH = 4. * st1 * st1 / (st2 * st2 * st2 * st2) * a02;
            // CelK calculation with linear interpolation.
            // CelK is the correction of the elastic electron
            // scatt. on molecular hydrogen compared to the independent atom
            // model.
            if (K < 3.) {
                i = (int) (K / 0.1);
                Ki = i * 0.1;
                CelK = Cel[i] + (K - Ki) / 0.1 * (Cel[i + 1] - Cel[i]);
            } else if (K >= 3. && K < 5.) {
                i = (int) (30 + (K - 3.) / 0.2);
                Ki = 3. + (i - 30) * 0.2;
                CelK = Cel[i] + (K - Ki) / 0.2 * (Cel[i + 1] - Cel[i]);
            } else if (K >= 5. && K < 9.49) {
                i = (int) (40 + (K - 5.) / 0.5);
                Ki = 5. + (i - 40) * 0.5;
                CelK = Cel[i] + (K - Ki) / 0.5 * (Cel[i + 1] - Cel[i]);
            } else {
                CelK = 0.;
            }
            Delreturn = 2. * gam * gam * DH * (1. + sin(K * d) / (K * d)) * (1. + CelK);
        } else {
            theta = acos(c) * 180. / Math.PI;
            for (i = 0; i <= 8; i++) {
                if (E >= e[i] && E < e[i + 1]) {
                    for (j = 0; j <= 8; j++) {
                        if (theta >= t[j] && theta < t[j + 1]) {
                            Delreturn = 1.e-20 * (D[i][j] + (D[i][j + 1] - D[i][j])
                                    * (theta - t[j]) / (t[j + 1] - t[j]));
                        }
                    }
                }
            }
        }
        return Delreturn;
    }

    double Dexc(double E, double c) {
        // This subroutine computes the differential cross section
        // Del= d sigma/d Omega  of excitation electron scattering
        // on molecular hydrogen.
        // Input:  E= electron kinetic energy in eV
        //         c= cos(theta), where theta is the polar scatt. angle
        // Dexc: in m^2/steradian
        double K2, K, sigma = 0., T, theta;
        double EE = 12.6 / 27.2;
        double[] e = {0., 25., 35., 50., 100.};
        double[] t = {0., 10., 20., 30., 40., 60., 80., 100., 180.};
        double[][] D = {
            {60., 43., 27., 18., 13., 8., 6., 6., 6.},
            {95., 70., 21., 9., 6., 3., 2., 2., 2.,},
            {150., 120., 32., 8., 3.7, 1.9, 1.2, 0.8, 0.8},
            {400., 200., 12., 2., 1.4, 0.7, 0.3, 0.2, 0.2}
        };
        int i, j;
        //
        T = E / 27.2;
        if (E >= 100.) {
            K2 = 4. * T * (1. - EE / (2. * T) - sqrt(1. - EE / T) * c);
            if (K2 < 1.e-9) {
                K2 = 1.e-9;
            }
            K = sqrt(K2); // momentum transfer
            sigma = 2. / K2 * sumexc(K) * a02;
        } else if (E <= 10.) {
            sigma = 0.;
        } else {
            theta = acos(c) * 180. / Math.PI;
            for (i = 0; i <= 3; i++) {
                if (E >= e[i] && E < e[i + 1]) {
                    for (j = 0; j <= 7; j++) {
                        if (theta >= t[j] && theta < t[j + 1]) {
                            sigma = 1.e-22 * (D[i][j] + (D[i][j + 1] - D[i][j])
                                    * (theta - t[j]) / (t[j + 1] - t[j]));
                        }
                    }
                }
            }
        }
        return sigma;
    }

    double Dinel(double E, double c) {
        // This subroutine computes the differential cross section
        // Dinel= d sigma/d Omega  of  inelastic electron scattering
        // on molecular hydrogen, within the first Born approximation.
        // Input:  E= electron kinetic energy in eV
        //         c= cos(theta), where theta is the polar scatt. angle
        // Dinel: in m2/steradian
        double[] Cinel = {
            -0.246, -0.244, -0.239, -0.234, -0.227,
            -0.219, -0.211, -0.201, -0.190, -0.179,
            -0.167, -0.155, -0.142, -0.130, -0.118,
            -0.107, -0.096, -0.085, -0.076, -0.067,
            -0.059, -0.051, -0.045, -0.039, -0.034,
            -0.029, -0.025, -0.022, -0.019, -0.016,
            -0.014, -0.010, -0.008, -0.006, -0.004,
            -0.003, -0.003, -0.002, -0.002, -0.001,
            -0.001, -0.001, 0.000, 0.000, 0.000,
            0.000, 0.000, 0.000, 0.000, 0.000
        };
        double Ei = 0.568;
        double T, K2, K, st1, F, DH, Dinelreturn, CinelK, Ki;
        int i;
        if (E < 16.) {
            return Dexc(E, c);
        }
        T = E / 27.2;
        K2 = 4. * T * (1. - Ei / (2. * T) - sqrt(1. - Ei / T) * c);
        if (K2 < 1.e-9) {
            K2 = 1.e-9;
        }
        K = sqrt(K2); // momentum transfer
        st1 = 1. + K2 / 4.;
        F = 1. / (st1 * st1); // scatt. formfactor of hydrogen atom
        // DH is the diff. cross section for inelastic electron scatt.
        // on atomic hydrogen within the first Born approximation :
        DH = 4. / (K2 * K2) * (1. - F * F) * a02;
        // CinelK calculation with linear interpolation.
        // CinelK is the correction of the inelastic electron
        // scatt. on molecular hydrogen compared to the independent atom
        // model.
        if (K < 3.) {
            i = (int) (K / 0.1);
            Ki = i * 0.1;
            CinelK = Cinel[i] + (K - Ki) / 0.1 * (Cinel[i + 1] - Cinel[i]);
        } else if (K >= 3. && K < 5.) {
            i = (int) (30 + (K - 3.) / 0.2);
            Ki = 3. + (i - 30) * 0.2;
            CinelK = Cinel[i] + (K - Ki) / 0.2 * (Cinel[i + 1] - Cinel[i]);
        } else if (K >= 5. && K < 9.49) {
            i = (int) (40 + (K - 5.) / 0.5);
            Ki = 5. + (i - 40) * 0.5;
            CinelK = Cinel[i] + (K - Ki) / 0.5 * (Cinel[i + 1] - Cinel[i]);
        } else {
            CinelK = 0.;
        }
        Dinelreturn = 2. * DH * (1. + CinelK);
        return Dinelreturn;
    }

    double sumexc(double K) {
        double[] Kvec = {0., 0.1, 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.5, 1.8, 2., 2.5, 3., 4., 5.};
        double[][] fvec = {
            {2.907e-1, 2.845e-1, 2.665e-1, 2.072e-1, 1.389e-1, // B
                8.238e-2, 4.454e-2, 2.269e-2, 7.789e-3, 2.619e-3,
                1.273e-3, 2.218e-4, 4.372e-5, 2.889e-6, 4.247e-7},
            {3.492e-1, 3.367e-1, 3.124e-1, 2.351e-1, 1.507e-1, // C
                8.406e-2, 4.214e-2, 1.966e-2, 5.799e-3, 1.632e-3,
                6.929e-4, 8.082e-5, 9.574e-6, 1.526e-7, 7.058e-9},
            {6.112e-2, 5.945e-2, 5.830e-2, 5.072e-2, 3.821e-2, // Bp
                2.579e-2, 1.567e-2, 8.737e-3, 3.305e-3, 1.191e-3,
                6.011e-4, 1.132e-4, 2.362e-5, 1.603e-6, 2.215e-7},
            {2.066e-2, 2.127e-2, 2.137e-2, 1.928e-2, 1.552e-2, // Bpp
                1.108e-2, 7.058e-3, 4.069e-3, 1.590e-3, 5.900e-4,
                3.046e-4, 6.142e-5, 1.369e-5, 9.650e-7, 1.244e-7},
            {9.405e-2, 9.049e-2, 8.613e-2, 7.301e-2, 5.144e-2, // D
                3.201e-2, 1.775e-2, 8.952e-3, 2.855e-3, 8.429e-4,
                3.655e-4, 4.389e-5, 5.252e-6, 9.010e-8, 7.130e-9},
            {4.273e-2, 3.862e-2, 3.985e-2, 3.362e-2, 2.486e-2, // Dp
                1.612e-2, 9.309e-3, 4.856e-3, 1.602e-3, 4.811e-4,
                2.096e-4, 2.498e-5, 2.905e-6, 5.077e-8, 6.583e-9},
            {0.000e-3, 2.042e-3, 7.439e-3, 2.200e-2, 3.164e-2, // EF
                3.161e-2, 2.486e-2, 1.664e-2, 7.562e-3, 3.044e-3,
                1.608e-3, 3.225e-4, 7.120e-5, 6.290e-6, 1.066e-6}};
        double[] EeV = {12.73, 13.20, 14.77, 15.3, 14.93, 15.4, 13.06};
        int n, j, jmin = 0, nmax;
        double En, sum;
        double[] f = new double[7];
        double[] x4 = new double[4];
        double[] f4 = new double[4];

        //
        sum = 0.;
        nmax = 6;
        for (n = 0; n <= nmax; n++) {
            En = EeV[n] / 27.21; // En is the excitation energy in Hartree atomic units
            if (K >= 5.) {
                f[n] = 0.;
            } else if (K >= 3. && K <= 4.) {
                f[n] = fvec[n][12] + (K - 3.) * (fvec[n][13] - fvec[n][12]);
            } else if (K >= 4. && K <= 5.) {
                f[n] = fvec[n][13] + (K - 4.) * (fvec[n][14] - fvec[n][13]);
            } else {
                for (j = 0; j < 14; j++) {
                    if (K >= Kvec[j] && K <= Kvec[j + 1]) {
                        jmin = j - 1;
                    }
                }
                if (jmin < 0) {
                    jmin = 0;
                }
                if (jmin > 11) {
                    jmin = 11;
                }
                for (j = 0; j <= 3; j++) {
                    x4[j] = Kvec[jmin + j];
                    f4[j] = fvec[n][jmin + j];
                }
                f[n] = lagrange(4, x4, f4, K);
            }
            sum += f[n] / En;
        }
        return sum;
    }

    double lagrange(int n, double[] xn, double[] fn, double x) {
        int i, j;
        double f, aa, bb;
        double[] a = new double[100];
        double[] b = new double[100];
        f = 0.;
        for (j = 0; j < n; j++) {
            for (i = 0; i < n; i++) {
                a[i] = x - xn[i];
                b[i] = xn[j] - xn[i];
            }
            a[j] = b[j] = aa = bb = 1.;
            for (i = 0; i < n; i++) {
                aa = aa * a[i];
                bb = bb * b[i];
            }
            f += fn[j] * aa / bb;
        }
        return f;
    }

    double sigmael(double E) {
        // This function computes the total elastic cross section of
        // electron scatt. on molecular hydrogen.
        // See: Liu, Phys. Rev. A35 (1987) 591,
        //      Trajmar, Phys Reports 97 (1983) 221.
        // E: incident electron energy in eV
        // sigmael: cross section in m^2
        double[] e = {0., 1.5, 5., 7., 10., 15., 20., 30., 60., 100., 150., 200., 300., 400.};
        double[] s = {9.6, 13., 15., 12., 10., 7., 5.6, 3.3, 1.1, 0.9, 0.5, 0.36, 0.23, 0.15};
        double gam, sigma = 0., T;
        int i;
        T = E / 27.2;
        if (E >= 400.) {
            gam = (emass + T) / emass;
            sigma = gam * gam * Math.PI / (2. * T) * (4.2106 - 1. / T) * a02;
        } else {
            for (i = 0; i <= 12; i++) {
                if (E >= e[i] && E < e[i + 1]) {
                    sigma = 1.e-20 * (s[i] + (s[i + 1] - s[i]) * (E - e[i]) / (e[i + 1] - e[i]));
                }
            }
        }
        return sigma;
    }

    double sigmaexc(double E) {
        // This function computes the electronic excitation cross section of
        // electron scatt. on molecular hydrogen.
        // E: incident electron energy in eV,
        // sigmaexc: cross section in m^2
        double sigma;
        if (E < 9.8) {
            sigma = 1.e-40;
        } else if (E >= 9.8 && E <= 250.) {
            sigma = sigmaBC(E) + sigmadiss10(E) + sigmadiss15(E);
        } else {
            sigma = 4. * Math.PI * a02 * R / E * (0.80 * log(E / R) + 0.28);
        }
        //    sigma=sigmainel(E)-sigmaion(E);
        return sigma;
    }

    double sigmaion(double E) {
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
        double B = 15.43, U = 15.98;
        double sigma, t, u, S, r, lnt;
        if (E < 16.) {
            sigma = 1.e-40;
        } else if (E >= 16. && E
                <= 250.) {
            t = E / B;
            u = U / B;
            r = R / B;
            S = 4. * Math.PI * a02 * 2. * r * r;
            lnt = log(t);
            sigma = S / (t + u + 1.) * (lnt / 2. * (1. - 1. / (t * t)) + 1. - 1. / t - lnt / (t + 1.));
        } else {
            sigma = 4. * Math.PI * a02 * R / E * (0.82 * log(E / R) + 1.3);
        }
        return sigma;
    }

    double sigmaBC(double E) {
        // This function computes the sigmaexc electronic excitation
        // cross section to the B and C states, with energy loss
        // about 12.5 eV.
        // E is incident electron energy in eV,
        // sigmaexc in m^2
        double[] aB = {-4.2935194e2, 5.1122109e2, -2.8481279e2,
            8.8310338e1, -1.6659591e1, 1.9579609,
            -1.4012824e-1, 5.5911348e-3, -9.5370103e-5};
        double[] aC = {-8.1942684e2, 9.8705099e2, -5.3095543e2,
            1.5917023e2, -2.9121036e1, 3.3321027,
            -2.3305961e-1, 9.1191781e-3, -1.5298950e-4};
        double lnsigma, lnE, lnEn, sigmaB, Emin, sigma, sigmaC;
        int n;
        sigma = 0.;
        Emin = 12.5;
        lnE = log(E);
        lnEn = 1.;
        lnsigma = 0.;
        if (E < Emin) {
            sigmaB = 0.;
        } else {
            for (n = 0; n <= 8; n++) {
                lnsigma += aB[n] * lnEn;
                lnEn = lnEn * lnE;
            }
            sigmaB = exp(lnsigma);
        }
        sigma += sigmaB;
        //  sigma=0.;
        // C state:
        Emin = 15.8;
        lnE = log(E);
        lnEn = 1.;
        lnsigma = 0.;
        if (E < Emin) {
            sigmaC = 0.;
        } else {
            for (n = 0; n <= 8; n++) {
                lnsigma += aC[n] * lnEn;
                lnEn = lnEn * lnE;
            }
            sigmaC = exp(lnsigma);
        }
        sigma += sigmaC;
        return sigma * 1.e-4;
    }

//////////////////////////////////////////////////////////////////
    double sigmadiss10(double E) {
        // This function computes the sigmadiss10 electronic
        // dissociative excitation
        // cross section, with energy loss
        // about 10 eV.
        // E is incident electron energy in eV,
        // sigmadiss10 in m^2
        double[] a = {-2.297914361e5, 5.303988579e5, -5.316636672e5,
            3.022690779e5, -1.066224144e5, 2.389841369e4,
            -3.324526406e3, 2.624761592e2, -9.006246604};
        double lnsigma, lnE, lnEn, Emin, sigma;
        int n;
        //  E is in eV
        sigma = 0.;
        Emin = 9.8;
        lnE = log(E);
        lnEn = 1.;
        lnsigma = 0.;
        if (E < Emin) {
            sigma = 0.;
        } else {
            for (n = 0; n <= 8; n++) {
                lnsigma += a[n] * lnEn;
                lnEn = lnEn * lnE;
            }
            sigma = exp(lnsigma);
        }
        return sigma * 1.e-4;
    }

//////////////////////////////////////////////////////////////////
    double sigmadiss15(double E) {
        // This function computes the sigmadiss15 electronic
        // dissociative excitation
        // cross section, with energy loss
        // about 15 eV.
        // E is incident electron energy in eV,
        // sigmadiss15 in m^2
        double[] a = {-1.157041752e3, 1.501936271e3, -8.6119387e2,
            2.754926257e2, -5.380465012e1, 6.573972423,
            -4.912318139e-1, 2.054926773e-2, -3.689035889e-4};
        double lnsigma, lnE, lnEn, Emin, sigma;
        int n;
        //  E is in eV
        sigma = 0.;
        Emin = 16.5;
        lnE = log(E);
        lnEn = 1.;
        lnsigma = 0.;
        if (E < Emin) {
            sigma = 0.;
        } else {
            for (n = 0; n <= 8; n++) {
                lnsigma += a[n] * lnEn;
                lnEn = lnEn * lnE;
            }
            sigma = exp(lnsigma);
        }
        return sigma * 1.e-4;
    }

    /**
     * Полное сечение с учетом квазиупругих столкновений
     *
     * @param E
     * @return
     */
    public double sigmaTotal(double E) {
        return sigmael(E) + sigmaexc(E) + sigmaion(E);

    }
}
