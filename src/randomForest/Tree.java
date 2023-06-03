package randomForest;

import java.util.Random;

public class Tree {

    private static int NODE_TERMINAL = -1;

    /// <summary>
    ///
    /// </summary>
    /// <param name="x"></param>
    /// <param name="n"></param>
    /// <param name="mdim"></param>
    /// <param name="treemap"></param>
    /// <param name="nodestatus"></param>
    /// <param name="xbestsplit"></param>
    /// <param name="bestvar"></param>
    /// <param name="nodeclass"></param>
    /// <param name="treeSize"></param>
    /// <param name="cat"></param>
    /// <param name="nclass"></param>
    /// <param name="jts"></param>
    /// <param name="nodex"></param>
    /// <param name="maxcat"></param>
    public static void predictClassTree(double[] x, int n, int mdim,
            ScootArray<Integer> treemap,
            ScootArray<Integer> nodestatus, ScootArray<Double> xbestsplit,
            ScootArray<Integer> bestvar, ScootArray<Integer> nodeclass,
            int treeSize, ScootArray<Integer> cat, int nclass,
            ScootArray<Integer> jts, ScootArray<Integer> nodex, int maxcat) {
        int m, i, j, k;
        int[] cbestsplit = null;
        int npack;

        /* decode the categorical splits */
        if (maxcat > 1) {
            cbestsplit = new int[maxcat * treeSize];
            for (i = 0; i < treeSize; i++) {
                if (nodestatus.get(i) != NODE_TERMINAL) {
                    if (cat.get(bestvar.get(i) - 1) > 1) {
                        npack = (int) Math.round(xbestsplit.get(i));
                        /* unpack `npack' into bits */
                        for (j = 0; npack > 0; npack >>= 1, j++) {
                            if (npack == 0)
                                cbestsplit[j + i * maxcat] = 0;
                            else
                                cbestsplit[j + i * maxcat] = 1;
                        }
                    }
                }
            }
        }
        for (i = 0; i < n; i++) {
            k = 0;
            while (nodestatus.get(k) != NODE_TERMINAL) {
                m = bestvar.get(k) - 1;
                if (cat.get(m) == 1) {
                    /* Split by a numerical predictor */
                    k = (x[m + i * mdim] <= xbestsplit.get(k))
                            ? treemap.get(k * 2) - 1
                            : treemap.get(1 + k * 2) - 1;
                } else {
                    /* Split by a categorical predictor */
                    k = (cbestsplit[(int) x[m + i * mdim] - 1 + k * maxcat] > 0)
                            ? treemap.get(k * 2) - 1
                            : treemap.get(1 + k * 2) - 1;
                }
            }
            /* Terminal node: assign class label */
            jts.set(i, nodeclass.get(k));
            nodex.set(i, k + 1);
        }
        // if (maxcat > 1)
        // Free(cbestsplit);
    }

    /// <summary>
    /// This finds the best split of a categorical variable with
    /// <code>lcat</code> categories
    /// and <code>nclass</code> classes The method uses an exhaustive search
    /// over all partitions
    /// of the category values if the number of categories is 10 or fewer.
    /// Otherwise ncsplit randomly
    /// selected splits are tested and best used.
    /// </summary>
    /// <param name="parentDen"></param>
    /// <param name="tclasscat"><code>tclasscat(j,k) is the number of cases in
    /// class <code>j</code> with category value <code>k</code></param>
    /// <param name="tclasspop"></param>
    /// <param name="nclass">number of classes</param>
    /// <param name="lcat">number of categories</param>
    /// <param name="ncatsp"></param>
    /// <param name="critmax"></param>
    /// <param name="nhit"></param>
    /// <param name="maxcat"></param>
    /// <param name="ncmax"></param>
    /// <param name="ncsplit">The number of radnomly selected splits</param>
    /// <param name="oRandom">The random number generator</param>
    public static void catmax(double parentDen, ScootArray<Double> tclasscat,
            ScootArray<Double> tclasspop, int nclass, int lcat,
                          int ncatsp, double critmax, int nhit,
                          int maxcat, int ncmax, int ncsplit,
                          Random oRandom) {
        int j, k, n, nsplit;
        int[] icat = new int[32];
        double leftNum, leftDen, rightNum, decGini;
        double[] leftCatClassCount = null;

        leftCatClassCount = new double[nclass]; //(double *) Calloc(*nclass, double);
        nhit = 0;
        nsplit = (lcat > ncmax) ?
            ncsplit : (int) (Math.pow(2.0, (double) (lcat - 1)) - 1.0);

        for (n = 0; n < nsplit; ++n) {
            if (lcat > ncmax) {
                /* Generate random split.
                   TODO: consider changing to generating random bits with more
                   efficient algorithm */
                for (j = 0; j < lcat; ++j) 
                    //icat[j] = unif_rand() > 0.5 ? 1 : 0;
                    icat[j] = oRandom.nextDouble() > 0.5 ? 1 : 0;
            } else {
                rfutils.unpack(lcat, n + 1, icat);
            }
            for (j = 0; j < nclass; ++j) {
                leftCatClassCount[j] = 0;
                for (k = 0; k < lcat; ++k) {
                    if (icat[k] > 0) {
                        leftCatClassCount[j] += tclasscat.get(j + k * nclass);
                    }
                }
            }
            leftNum = 0.0;
            leftDen = 0.0;
            for (j = 0; j < nclass; ++j) {
                leftNum += leftCatClassCount[j] * leftCatClassCount[j];
                leftDen += leftCatClassCount[j];
            }
            /* If either node is empty, try another split. */
            if (leftDen <= 1.0e-8 || parentDen - leftDen <= 1.0e-5) continue;
            rightNum = 0.0;
            for (j = 0; j < nclass; ++j) {
                leftCatClassCount[j] = tclasspop.get(j) - leftCatClassCount[j];
                rightNum += leftCatClassCount[j] * leftCatClassCount[j];
            }
            decGini = (leftNum / leftDen) + (rightNum / (parentDen - leftDen));
            if (decGini > critmax) {
                critmax = decGini;
                nhit = 1;
                ncatsp = (lcat > ncmax) ? rfutils.pack(lcat, icat) : n + 1;
            }
        }

    // Free(leftCatClassCount);
    }

    /// <summary>
    /// Find the best split of with categorical variable when there are two
    /// classes
    /// </summary>
    /// <param name="totalWt"></param>
    /// <param name="tclasscat"></param>
    /// <param name="classCount"></param>
    /// <param name="nclass"></param>
    /// <param name="nCat"></param>
    /// <param name="nbest"></param>
    /// <param name="critmax"></param>
    /// <param name="nhit"></param>
    /// <param name="catCount"></param>
    public static void catmaxb(double totalWt, ScootArray<Double> tclasscat, 
            ScootArray<Double> classCount, 
                           int nclass, int nCat, int nbest, double critmax,
                           int nhit, double[] catCount) {

        //double catProportion[32], cp[32], cm[32];
        double[] catProportion = new double[32];
        double[] cp = new double[32];
        double[] cm = new double[32];
        int[] kcat = new int[32];
        int i, j;
        double bestsplit=0.0, rightDen, leftDen, leftNum, rightNum, crit;

        nhit = 0;
        for (i = 0; i < nCat; ++i) {
            catProportion[i] = (catCount[i] > 0) ?
                tclasscat.get(i * nclass) / catCount[i] : 0.0;
            kcat[i] = i + 1;
        }
        rfutils.sortArrays(catProportion, kcat); // sort paired arrays by the keys (catProportion)
        for (i = 0; i < nclass; ++i)
        {
            cp[i] = 0;
            cm[i] = classCount.get(i);
        }
        rightDen = totalWt;
        leftDen = 0.0;
        for (i = 0; i < nCat - 1; ++i) {
            leftDen += catCount[kcat[i]-1];
            rightDen -= catCount[kcat[i]-1];
            leftNum = 0.0;
            rightNum = 0.0;
            for (j = 0; j < nclass; ++j) {
                cp[j] += tclasscat.get(j + (kcat[i]-1) * nclass);
                cm[j] -= tclasscat.get(j + (kcat[i]-1) * nclass);
                leftNum += cp[j] * cp[j];
                rightNum += cm[j] * cm[j];
            }
            if (catProportion[i] < catProportion[i + 1]) {
                /* If neither node is empty, check the split. */
                if (rightDen > 1.0e-5 && leftDen > 1.0e-5) {
                    crit = (leftNum / leftDen) + (rightNum / rightDen);
                    if (crit > critmax) {
                        critmax = crit;
                        bestsplit = .5 * (catProportion[i] + catProportion[i + 1]);
                        nhit = 1;
                    }
                }
            }
        }
        if (nhit == 1) {
            for (i = 0; i < nCat; ++i) {
                catProportion[i] = (catCount[i] > 0) ?
                    tclasscat.get(i * nclass) / catCount[i] : 0.0;
                kcat[i] = catProportion[i] < bestsplit ? 1 : 0;
                /* Rprintf("%i ", kcat[i]); */
            }
            nbest = rfutils.pack(nCat, kcat);
            /* Rprintf("\nnbest=%u\nnbest=%i\n", *nbest, *nbest); */
        }
}
}
