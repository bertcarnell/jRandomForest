package randomForest;

import java.util.Arrays;
import java.util.Random;
import java.util.TreeMap;

/// <summary>
/// Utility methods for the random forest code (rf.cs)
/// </summary>
public class rfutils {

    /// <summary>
    /// Fill an array with zero
    /// </summary>
    /// <param name="x">the array to fill</param>
    /// <param name="length">the length of the locations to fill (not
    /// necessarily the full length)</param>
    public static void zeroInt(ScootArray<Integer> x, int length) {
        if (x.size() >= length) {
            for (int i = 0; i < length; i++)
                x.set(i, 0);
        }
    }

    /// <summary>
    /// Fill an array with zero
    /// </summary>
    /// <param name="x">the array to fill</param>
    /// <param name="length">the length of the locations to fill (not
    /// necessarily the full length)</param>
    public static void zeroInt(int[] x, int length) {
        if (x.length >= length) {
            for (int i = 0; i < length; i++)
                x[i] = 0;
        }
    }

    /// <summary>
    /// Fill an array with zero
    /// </summary>
    /// <param name="x">the array to fill</param>
    /// <param name="length">the length of the locations to fill (not
    /// necessarily the full length)</param>
    public static void zeroDouble(double[] x, int length) {
        if (x.length >= length) {
            for (int i = 0; i < length; i++)
                x[i] = 0.0;
        }
    }

    /// <summary>
    /// Fill an array with zero
    /// </summary>
    /// <param name="x">the array to fill</param>
    /// <param name="length">the length of the locations to fill (not
    /// necessarily the full length)</param>
    public static void zeroDouble(ScootArray<Double> x, int length) {
        if (x.size() >= length) {
            for (int i = 0; i < length; i++)
                x.set(i, 0.0);
        }
    }

    /// <summary>
    /// Normalize Class Weights
    /// </summary>
    /// <param name="cl"></param>
    /// <param name="nsample"></param>
    /// <param name="nclass"></param>
    /// <param name="useWt"></param>
    /// <param name="classwt"></param>
    /// <param name="classFreq"></param>
    public static void normClassWt(ScootArray<Integer> cl, int nsample,
            int nclass,
            int useWt, double[] classwt, int[] classFreq) {
        int i;
        double sumwt = 0.0;

        if (useWt > 0) {
            /* Normalize user-supplied weights so they sum to one. */
            for (i = 0; i < nclass; ++i)
                sumwt += classwt[i];
            for (i = 0; i < nclass; ++i)
                classwt[i] /= sumwt;
        } else {
            for (i = 0; i < nclass; ++i) {
                classwt[i] = ((double) classFreq[i]) / nsample;
            }
        }
        for (i = 0; i < nclass; ++i) {
            classwt[i] = (classFreq[i] > 0)
                    ? classwt[i] * nsample / classFreq[i]
                    : 0.0;
        }
    }

    /// <summary>
    /// makeA() constructs the <code>mdim by nsample</code> integer array a. For
    /// each
    /// numerical variable with values x(m, n), n=1, ...,nsample, the x-values
    /// are sorted from lowest to highest. Denote these by xs(m, n). Then
    /// a(m,n) is the case number in which xs(m, n) occurs. The b matrix is
    /// also contructed here. If the mth variable is categorical, then
    /// a(m, n) is the category of the nth case number.
    /// </summary>
    /// <param name="x"></param>
    /// <param name="mdim"></param>
    /// <param name="nsample"></param>
    /// <param name="cat"></param>
    /// <param name="a"><code>mdim x nsample</code> integer array</param>
    /// <param name="b"></param>
    public static void makeA(double[] x, int mdim, int nsample,
            ScootArray<Integer> cat, int[] a,
            ScootArray<Integer> b) {
        int i, j, n1, n2;
        int[] index;
        double[] v;

        v = new double[nsample]; // v = (double *) Calloc(nsample, double);
        index = new int[nsample]; // index = (int *) Calloc(nsample, int);

        for (i = 0; i < mdim; ++i) {
            if (cat.get(i) == 1) { /* numerical predictor */
                for (j = 0; j < nsample; ++j) {
                    v[j] = x[i + j * mdim];
                    index[j] = j + 1;
                }
                // R_qsort_I(v, index, 1, nsample);
                sortArrays(v, index); // sort the items (index) by the keys (v)
                                      // and the keys (v)

                /*
                 * this sorts the v(n) in ascending order. index(n) is the case
                 * number of that v(n) nth from the lowest (assume the original
                 * case numbers are 1,2,...).
                 */
                for (j = 0; j < nsample - 1; ++j) {
                    n1 = index[j];
                    n2 = index[j + 1];
                    a[i + j * mdim] = n1;
                    if (j == 0)
                        b.set(i + (n1 - 1) * mdim, 1);
                    if (v[j] < v[j + 1]) {
                        b.set(i + (n2 - 1) * mdim,
                                b.get(i + (n1 - 1) * mdim) + 1);
                    } else {
                        b.set(i + (n2 - 1) * mdim, b.get(i + (n1 - 1) * mdim));
                    }
                }
                a[i + (nsample - 1) * mdim] = index[nsample - 1];
            } else { /* categorical predictor */
                for (j = 0; j < nsample; ++j)
                    a[i + j * mdim] = (int) x[i + j * mdim];
            }
        }
        // Free(index);
        // Free(v);
    }

    /// <summary>
    ///
    /// </summary>
    /// <param name="x"></param>
    /// <param name="realN"></param>
    /// <param name="totalN"></param>
    /// <param name="mdim"></param>
    /// <param name="oRandom">Random number generator</param>
    public static void createClass(double[] x, int realN, int totalN, int mdim,
            Random oRandom) {
        /*
         * Create the second class by bootstrapping each variable independently.
         */
        int i, j, k;
        for (i = realN; i < totalN; ++i) {
            for (j = 0; j < mdim; ++j) {
                // k = (int)unif_rand() * realN;
                k = (int) (Math.floor(oRandom.nextDouble() * (double) (realN)));
                x[j + i * mdim] = x[j + k * mdim];
            }
        }
    }

    /// <summary>
    ///
    /// </summary>
    /// <param name="a"></param>
    /// <param name="nuse"></param>
    /// <param name="nsample"></param>
    /// <param name="mdim"></param>
    /// <param name="cat"></param>
    /// <param name="maxcat"></param>
    /// <param name="ncase"></param>
    /// <param name="jin"></param>
    public static void modA(ScootArray<Integer> a, int nuse, int nsample,
            int mdim,
            ScootArray<Integer> cat, int maxcat, ScootArray<Integer> ncase,
            int[] jin) {
        int i, j, k, m, nt;

        nuse = 0;
        for (i = 0; i < nsample; ++i)
            if (jin[i] > 0)
                nuse++;

        for (i = 0; i < mdim; ++i) {
            k = 0;
            nt = 0;
            if (cat.get(i) == 1) {
                for (j = 0; j < nsample; ++j) {
                    if (jin[a.get(i + k * mdim) - 1] > 0) {
                        a.set(i + nt * mdim, a.get(i + k * mdim));
                        k++;
                    } else {
                        for (m = 0; m < nsample - k; ++m) {
                            if (jin[a.get(i + (k + m) * mdim) - 1] > 0) {
                                a.set(i + nt * mdim, a.get(i + (k + m) * mdim));
                                k += m + 1;
                                break;
                            }
                        }
                    }
                    nt++;
                    if (nt >= nuse)
                        break;
                }
            }
        }
        if (maxcat > 1) {
            k = 0;
            nt = 0;
            for (i = 0; i < nsample; ++i) {
                if (jin[k] > 0) {
                    k++;
                    ncase.set(nt, k);
                } else {
                    for (j = 0; j < nsample - k; ++j) {
                        if (jin[k + j] > 0) {
                            ncase.set(nt, k + j + 1);
                            k += j + 1;
                            break;
                        }
                    }
                }
                nt++;
                if (nt >= nuse)
                    break;
            }
        }
    }

    /// <summary>
    /// This subroutine takes the splits on numerical variables and translates
    /// them
    /// back into x-values. It also unpacks each categorical split into a
    /// 32-dimensional vector with components of zero or one--a one indicates
    /// that
    /// the corresponding category goes left in the split.
    /// </summary>
    /// <param name="x"></param>
    /// <param name="mdim"></param>
    /// <param name="nrnodes"></param>
    /// <param name="nsample"></param>
    /// <param name="bestvar"></param>
    /// <param name="bestsplit"></param>
    /// <param name="bestsplitnext"></param>
    /// <param name="xbestsplit"></param>
    /// <param name="nodestatus"></param>
    /// <param name="cat"></param>
    /// <param name="treeSize"></param>
    public static void Xtranslate(double[] x, int mdim, int nrnodes,
            int nsample,
            ScootArray<Integer> bestvar, ScootArray<Integer> bestsplit,
            ScootArray<Integer> bestsplitnext,
            ScootArray<Double> xbestsplit, ScootArray<Integer> nodestatus,
            ScootArray<Integer> cat, int treeSize) {
        int i, m;

        for (i = 0; i < treeSize; ++i) {
            if (nodestatus.get(i) == 1) {
                m = bestvar.get(i) - 1;
                if (cat.get(m) == 1) {
                    xbestsplit.set(i,
                            (0.5 * (x[m + (bestsplit.get(i) - 1) * mdim] +
                                    x[m + (bestsplitnext.get(i) - 1) * mdim])));
                } else {
                    xbestsplit.set(i, (double) bestsplit.get(i));
                }
            }
        }
    }

    /* Compute proximity. */
    public static void computeProximity(double[] prox, int oobprox,
            ScootArray<Integer> node,
            int[] inbag, int[] oobpair, int n) {
        /*
         * Accumulate the number of times a pair of points fall in the same
         * node. prox: n x n proximity matrix oobprox: should the accumulation
         * only count OOB cases? (0=no, 1=yes) node: vector of terminal node
         * labels inbag: indicator of whether a case is in-bag oobpair: matrix
         * to accumulate the number of times a pair is OOB together n: total
         * number of cases
         */
        int i, j;
        for (i = 0; i < n; ++i) {
            for (j = i + 1; j < n; ++j) {
                if (oobprox > 0) {
                    if ((inbag[i] > 0) ^ (inbag[j] > 0)) {
                        oobpair[j * n + i]++;
                        oobpair[i * n + j]++;
                        if (node.get(i) == node.get(j)) {
                            prox[j * n + i] += 1.0;
                            prox[i * n + j] += 1.0;
                        }
                    }
                } else {
                    if (node.get(i) == node.get(j)) {
                        prox[j * n + i] += 1.0;
                        prox[i * n + j] += 1.0;
                    }
                }
            }
        }
    }

    /// <summary>
    /// Permute the OOB part of a variable in x.
    /// </summary>
    /// <param name="m">the variable to be permuted</param>
    /// <param name="x">the data matrix (variables in rows)</param>
    /// <param name="_in">vector indicating which case is OOB</param>
    /// <param name="nsample">number of cases in the data</param>
    /// <param name="mdim">number of variables in the data</param>
    /// <param name="oRandom">Random number generator</param>
    public static void permuteOOB(int m, double[] x, int[] _in, int nsample,
            int mdim, Random oRandom) {
        double tmp;
        double[] tp;
        int i, last, k, nOOB = 0;

        tp = new double[nsample]; // (double *) Calloc(nsample, double);

        for (i = 0; i < nsample; ++i) {
            /*
             * make a copy of the OOB part of the data into tp (for permuting)
             */
            if (_in[i] == 0) {
                tp[nOOB] = x[m + i * mdim];
                nOOB++;
            }
        }
        /* Permute tp */
        last = nOOB;
        for (i = 0; i < nOOB; ++i) {
            /*
             * #ifdef USER k = (int) last * unif_rand(); #else k = (int)
             * ((double)last * unif_rand()); #endif
             */
            k = (int) (Math.floor((double) (last) * oRandom.nextDouble()));
            tmp = tp[last - 1];
            tp[last - 1] = tp[k];
            tp[k] = tmp;
            last--;
        }

        /* Copy the permuted OOB data back into x. */
        nOOB = 0;
        for (i = 0; i < nsample; ++i) {
            if (_in[i] == 0) {
                x[m + i * mdim] = tp[nOOB];
                nOOB++;
            }
        }
        // Free(tp);
    }

    /// <summary>
    ///
    /// </summary>
    /// <param name="nBits"></param>
    /// <param name="bits"></param>
    /// <returns></returns>
    public static int pack(int nBits, int[] bits) {
        int i = nBits;
        int pack = 0;
        while (i-- >= 0)
            pack += (int) (bits[i] << i);
        return pack;
    }

    /// <summary>
    /// The subroutine returns icat, an integer array of zeroes and ones
    /// corresponding to the coefficients
    /// in the binary expansion of pack
    /// </summary>
    /// <param name="nBits"></param>
    /// <param name="pack">4-byte integer</param>
    /// <param name="bits"></param>
    public static void unpack(int nBits, int pack, int[] bits) {
        int i;
        for (i = 0; i < nBits; pack >>= 1, ++i)
            bits[i] = (int) (pack & 1);
    }

    public static void sortArrays(double[] key, int[] data) {
        TreeMap<Double, Integer> helper = new TreeMap<Double, Integer>();
        for (int i = 0; i < key.length; i++) {
            helper.put(key[i], data[i]);
        }
        Arrays.sort(key);
        for (int i = 0; i < data.length; i++) {
            data[i] = helper.get(key[i]);
        }
    }
}
