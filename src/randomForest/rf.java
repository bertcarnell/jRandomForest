package randomForest;

import java.util.ArrayList;
import java.util.Random;

public class rf {
    private static void swapInt(int a, int b) {
        int oldb = b;
        b = a;
        a = oldb;
    }

    private static void TestSetError(double[] countts, ScootArray<Integer> jts,
            int[] clts,
            int[] jet, int ntest,
            int nclass, int nvote, ScootArray<Double> errts,
            int labelts, int[] nclts, double[] cutoff, Random oRandom) {
        int j, n, ntie;
        double cmax, crit;

        for (n = 0; n < ntest; ++n)
            countts[jts.get(n) - 1 + n * nclass] += 1.0;

        /* Prediction is the class with the maximum votes */
        for (n = 0; n < ntest; ++n) {
            cmax = 0.0;
            ntie = 1;
            for (j = 0; j < nclass; ++j) {
                crit = (countts[j + n * nclass] / nvote) / cutoff[j];
                if (crit > cmax) {
                    jet[n] = j + 1;
                    cmax = crit;
                    ntie = 1;
                }
                /* Break ties at random: */
                if (crit == cmax) {
                    // if (unif_rand() < 1.0 / ntie) {
                    if (oRandom.nextDouble() < 1.0 / (double) ntie) {
                        jet[n] = j + 1;
                        cmax = crit;
                    }
                    ntie++;
                }
            }
        }
        if (labelts > 0) {
            rfutils.zeroDouble((ScootArray<Double>) errts.copy(), nclass + 1);
            for (n = 0; n < ntest; ++n) {
                if (jet[n] != clts[n]) {
                    errts.set(0, errts.get(0) + 1.0);
                    errts.set(clts[n], clts[n] + 1.0);
                }
            }
            errts.set(0, errts.get(0) / ntest);
            for (n = 1; n <= nclass; ++n)
                errts.set(n, errts.get(n) / nclts[n - 1]);
        }
    }

    /// <summary>
    /// C# wrapper for random forests
    /// </summary>
    /// <param name="x">"matrix of predictors (transposed!)" - variables are
    /// rows, observations are columns</param>
    /// <param name="dimx">two integers: number of variables and number of
    /// cases</param>
    /// <param name="cl">class labels of the data</param>
    /// <param name="ncl">number of classes in the response</param>
    /// <param name="cat">integer vector of number of classes in the predictor:
    /// 1=continuous</param>
    /// <param name="maxcat">maximum of cat</param>
    /// <param name="sampsize">Sizes of samples to draw. For classification, if
    /// sampsize is a vector of the length
    /// of the number of strata, then sampling is stratified by strata, and teh
    /// elements of sampsize indicate the numbers
    /// to be drawn from the strata</param>
    /// <param name="strata">A variable that is used for stratified
    /// sampling</param>
    /// <param name="Options">
    /// <list type="num">
    /// <listheader>7 integers: (0=no, 1=yes)</listheader>
    /// <item>add a second class (for unsupervised RF)?</item>
    /// <list type="bull">
    /// <item>1: sampling from product of marginals</item>
    /// <item>2: sampling from product of uniforms</item>
    /// </list>
    /// <item>assess variable importance?</item>
    /// <item>calculate proximity?</item>
    /// <item>calculate proximity based on OOB predictions?</item>
    /// <item>calculate outlying measure?</item>
    /// <item>how often to print output?</item>
    /// <item>keep the forest for future prediction?</item>
    /// </list>
    /// </param>
    /// <param name="ntree">number of trees</param>
    /// <param name="nvar">number of predictors to use for each split</param>
    /// <param name="ipi">0=use class proportion as prob.; 1=use supplied
    /// priors</param>
    /// <param name="classwt">Priors of the classes. Need not add up to one.
    /// Ignored for regression</param>
    /// <param name="cut">(Classification Only) A vector of length equal to the
    /// number of classes. The 'winning'
    /// class for an observation is the one with teh maximum ratio of proportion
    /// of votes to cutoff. Default is
    /// 1/k where k is the number of classes (i.e., majority vote wins)</param>
    /// <param name="nodesize">minimum node size: no node with fewer than ndsize
    /// cases will be split</param>
    /// <param name="outcl">Output: class predicted by RF</param>
    /// <param name="counttr">Output: matrix of votes (transposed!)</param>
    /// <param name="prox">Output: matrix of proximity (if iprox=1)</param>
    /// <param name="imprt">Output: matrix of variable importance measures. A
    /// matrix with <code>nclass+2</code>
    /// (for classification) or two (for regression) columns. For
    /// classification, the first <code>nclass</code>
    /// columns are the class-specific measures computed as mean decrease in
    /// accuracy. The <code>nclass+1st</code>
    /// columns is the mean decrease in accuracy over all classes. Teh last
    /// column is the mean decrease in Gini index.
    /// For regression, the first column is the mean decrease in accuracy and
    /// the second column is the mean decrease in
    /// MSE. if <code>importance=false</code>, the last measure is still
    /// returned as a vector</param>
    /// <param name="impsd">Output: The "standard errors" of the
    /// permutation-based importance measure. For classification,
    /// a <code>p x (nclass + 1)</code> matrix corresponding to the first
    /// <code>nclass+1</code> columns of the importance
    /// matrix. For regression, a lenght <code>p</code> vector</param>
    /// <param name="impmat">Output: matrix of local variable importance
    /// measures. A <code>p x n</code> matrix containing
    /// the casewise importance measures, the <code>[i,j]</code> element of
    /// which is the importance of the <code>i-th</code>
    /// variable on the <code>j-th</code> case. NULL if
    /// <code>localImp=false</code></param>
    /// <param name="nrnodes"></param>
    /// <param name="ndbigtree"></param>
    /// <param name="nodestatus"></param>
    /// <param name="bestvar"></param>
    /// <param name="treemap"></param>
    /// <param name="nodeclass"></param>
    /// <param name="xbestsplit"></param>
    /// <param name="errtr"></param>
    /// <param name="testdat"></param>
    /// <param name="xts">Test data set</param>
    /// <param name="clts">Test response data (optional)</param>
    /// <param name="nts"></param>
    /// <param name="countts"></param>
    /// <param name="outclts">Predicted output class for the test data</param>
    /// <param name="labelts"></param>
    /// <param name="proxts">proximity for the test data</param>
    /// <param name="errts">error rate for the test data</param>
    /// <param name="inbag"></param>
    /// <param name="oUniformRandom">The random number generator used by the
    /// random forest</param>
    public static void classRF(double[] x, int[] dimx, ScootArray<Integer> cl,
            int ncl,
            ScootArray<Integer> cat, int maxcat,
            int[] sampsize, int[] strata, int[] Options, int ntree, int nvar,
            int ipi, double[] classwt, double[] cut, int nodesize,
            int[] outcl, int[] counttr, double[] prox,
            double[] imprt, double[] impsd, double[] impmat, int nrnodes,
            ScootArray<Integer> ndbigtree, ScootArray<Integer> nodestatus,
            ScootArray<Integer> bestvar,
            ScootArray<Integer> treemap,
            ScootArray<Integer> nodeclass, ScootArray<Double> xbestsplit,
            ScootArray<Double> errtr,
            int testdat, double[] xts, int[] clts, int nts, double[] countts,
            int[] outclts, int labelts, double[] proxts,
            ScootArray<Double> errts,
            int[] inbag, Random oUniformRandom) throws Exception {
        int nsample0, mdim, nclass, addClass, mtry, ntest, nsample, ndsize, 
                near, nuse, noutall, nrightall, nrightimpall,
                keepInbag, nstrata;
//        int mimp, nimp;
        int jb, j, n, m, k, idxByNnode, idxByNsample, imp, localImp, iprox,
                oobprox, keepf, replace, stratify, trace;
        int[] nright, nrightimp, nout, nclts;
        int Ntree;

        nstrata = 0;
        nuse = 0;

        int[] _out, jin, jerr, classFreq, at, nind, oobpair;
        oobpair = null;
        ArrayList<int[]> strata_idx = null;
        int[] strata_size = null;
        int last, ktmp, nEmpty, ntry;

        double av = 0.0;

        double[] tx;
//        double[] tp;

//        int[] a_mem, b_mem, bestsplit_mem, bestsplitnext_mem, nodepop_mem,
//                nodestart_mem, ta_mem, idmove_mem, ncase_mem, varUsed_mem,
//                mind_mem, nodex_mem, nodexts_mem, jts_mem, jtr_mem, jvr_mem;
        ScootArray<Integer> a, b, bestsplit, bestsplitnext, nodepop, nodestart,
                ta, idmove,
                ncase, varUsed, mind, nodex, nodexts, jts, jtr, jvr;
//        double[] tgini_mem, classpop_mem, tclasspop_mem, tclasscat_mem, win_mem,
//                wr_mem, wl_mem;
        ScootArray<Double> tgini, classpop, tclasspop, tclasscat, win, wr, wl;

        if (Options.length != 10)
            throw new Exception("error message 1");
        addClass = Options[0];
        imp = Options[1];
        localImp = Options[2];
        iprox = Options[3];
        oobprox = Options[4];
        trace = Options[5];
        keepf = Options[6];
        replace = Options[7];
        stratify = Options[8];
        keepInbag = Options[9];
        if (dimx.length != 2)
            throw new Exception("error message 2");
        mdim = dimx[0];
        nsample0 = dimx[1];
        nclass = (ncl == 1) ? 2 : ncl;
        ndsize = nodesize;
        Ntree = ntree;
        mtry = nvar;
        ntest = nts;
        nsample = (addClass > 0) ? (nsample0 + nsample0) : nsample0;
//        mimp = (imp > 0) ? mdim : 1;
//        nimp = (imp > 0) ? nsample : 1;
        near = (iprox > 0) ? nsample0 : 1;
        if (trace == 0)
            trace = Ntree + 1;

//        tgini_mem = new double[mdim];// (double *] S_alloc(mdim,
//                                     // sizeof(double]];
//        wl_mem = new double[nclass]; // (double *] S_alloc(nclass,
//                                     // sizeof(double]];
//        wr_mem = new double[nclass]; // (double *] S_alloc(nclass,
//                                     // sizeof(double]];
//        classpop_mem = new double[nclass * nrnodes];// (double *]
//                                                    // S_alloc(nclass* *nrnodes,
//                                                    // sizeof(double]];
//        tclasscat_mem = new double[nclass * 32];// (double *] S_alloc(nclass*32,
//                                                // sizeof(double]];
//        tclasspop_mem = new double[nclass];// (double *] S_alloc(nclass,
//                                           // sizeof(double]];
//        win_mem = new double[nsample];// (double *] S_alloc(nsample,
                                      // sizeof(double]];
        tgini = new ScootArray<Double>();// (double *) S_alloc(mdim,
                                              // sizeof(double));
        wl = new ScootArray<Double>(); // (double *) S_alloc(nclass,
                                         // sizeof(double));
        wr = new ScootArray<Double>(); // (double *) S_alloc(nclass,
                                         // sizeof(double));
        classpop = new ScootArray<Double>();// (double *)
                                                    // S_alloc(nclass*
                                                    // *nrnodes,
                                                    // sizeof(double));
        tclasscat = new ScootArray<Double>();// (double *)
                                                      // S_alloc(nclass*32,
                                                      // sizeof(double));
        tclasspop = new ScootArray<Double>();// (double *)
                                                      // S_alloc(nclass,
                                                      // sizeof(double));
        win = new ScootArray<Double>();// (double *) S_alloc(nsample,
                                          // sizeof(double));
        tx = new double[nsample];// (double *) S_alloc(nsample, sizeof(double));
//        tp = new double[nsample];// (double *) S_alloc(nsample, sizeof(double));

//        bestsplitnext_mem = new int[nrnodes];// (int *] S_alloc(*nrnodes,
//                                             // sizeof(int]];
//        bestsplit_mem = new int[nrnodes];// (int *] S_alloc(*nrnodes,
//                                         // sizeof(int]];
//        nodepop_mem = new int[nrnodes];// (int *] S_alloc(*nrnodes,
//                                       // sizeof(int]];
//        nodestart_mem = new int[nrnodes];// (int *] S_alloc(*nrnodes,
//                                         // sizeof(int]];
//        nodex_mem = new int[nsample];// (int *] S_alloc(nsample, sizeof(int]];
//        nodexts_mem = new int[ntest];// (int *] S_alloc(ntest, sizeof(int]];
//        ta_mem = new int[nsample];// (int *] S_alloc(nsample, sizeof(int]];
//        ncase_mem = new int[nsample];// (int *] S_alloc(nsample, sizeof(int]];
//        varUsed_mem = new int[mdim];// (int *] S_alloc(mdim, sizeof(int]];
//        jtr_mem = new int[nsample];// (int *] S_alloc(nsample, sizeof(int]];
//        jvr_mem = new int[nsample];// (int *] S_alloc(nsample, sizeof(int]];
//        jts_mem = new int[ntest];// (int *] S_alloc(ntest, sizeof(int]];
//        idmove_mem = new int[nsample];// (int *] S_alloc(nsample, sizeof(int]];
//        a_mem = new int[mdim * nsample];// (int *] S_alloc(mdim*nsample,
//                                        // sizeof(int]];
//        b_mem = new int[mdim * nsample];// (int *] S_alloc(mdim*nsample,
//                                        // sizeof(int]];
//        mind_mem = new int[mdim];// (int *] S_alloc(mdim, sizeof(int]];

        bestsplitnext = new ScootArray<Integer>();// (int *)
                                                               // S_alloc(*nrnodes,
                                                               // sizeof(int));
        bestsplit = new ScootArray<Integer>();// (int *)
                                                       // S_alloc(*nrnodes,
                                                       // sizeof(int));
        nodepop = new ScootArray<Integer>();// (int *) S_alloc(*nrnodes,
                                                   // sizeof(int));
        nodestart = new ScootArray<Integer>();// (int *)
                                                       // S_alloc(*nrnodes,
                                                       // sizeof(int));
        nodex = new ScootArray<Integer>();// (int *) S_alloc(nsample,
                                               // sizeof(int));
        nodexts = new ScootArray<Integer>();// (int *) S_alloc(ntest,
                                                   // sizeof(int));
        ta = new ScootArray<Integer>();// (int *) S_alloc(nsample,
                                         // sizeof(int));
        ncase = new ScootArray<Integer>();// (int *) S_alloc(nsample,
                                               // sizeof(int));
        varUsed = new ScootArray<Integer>();// (int *) S_alloc(mdim,
                                                   // sizeof(int));
        jtr = new ScootArray<Integer>();// (int *) S_alloc(nsample,
                                           // sizeof(int));
        jvr = new ScootArray<Integer>();// (int *) S_alloc(nsample,
                                           // sizeof(int));
        jts = new ScootArray<Integer>();// (int *) S_alloc(ntest,
                                           // sizeof(int));
        idmove = new ScootArray<Integer>();// (int *) S_alloc(nsample,
                                                 // sizeof(int));
        a = new ScootArray<Integer>();// (int *) S_alloc(mdim*nsample,
                                       // sizeof(int));
        b = new ScootArray<Integer>();// (int *) S_alloc(mdim*nsample,
                                       // sizeof(int));
        mind = new ScootArray<Integer>();// (int *) S_alloc(mdim,
                                             // sizeof(int));

        jin = new int[nsample];// (int *) S_alloc(nsample, sizeof(int));
        _out = new int[nsample];// (int *) S_alloc(nsample, sizeof(int));
        jerr = new int[nsample];// (int *) S_alloc(nsample, sizeof(int));
        classFreq = new int[nclass];// (int *) S_alloc(nclass, sizeof(int));
        at = new int[mdim * nsample];// (int *) S_alloc(mdim*nsample,
                                     // sizeof(int));
        nright = new int[nclass];// (int *) S_alloc(nclass, sizeof(int));
        nrightimp = new int[nclass];// (int *) S_alloc(nclass, sizeof(int));
        nout = new int[nclass];// (int *) S_alloc(nclass, sizeof(int));

        if (oobprox > 0) {
            oobpair = new int[near * near]; // (int *) S_alloc(near*near,
                                            // sizeof(int));
        }

        /* Count number of cases in each class. */
        rfutils.zeroInt(classFreq, nclass);
        for (n = 0; n < nsample; ++n)
            classFreq[cl.get(n) - 1]++;
        /* Normalize class weights. */
        rfutils.normClassWt((ScootArray<Integer>) cl.copy(), nsample, nclass,
                ipi, classwt, classFreq);

        nind = null;
        if (stratify > 0) {
            /* Count number of strata and frequency of each stratum. */
            nstrata = 0;
            for (n = 0; n < nsample0; ++n)
                if (strata[n] > nstrata)
                    nstrata = strata[n];
            /*
             * Create the array of pointers, each pointing to a vector of
             * indices of where data of each stratum is.
             */
            strata_size = new int[nstrata]; // (int *) S_alloc(nstrata,
                                            // sizeof(int));
            for (n = 0; n < nsample0; ++n) {
                strata_size[strata[n] - 1]++;
            }
            strata_idx = new ArrayList<int[]>(nstrata); // (int **)
                                                        // S_alloc(nstrata,
            // sizeof(int *));
            for (n = 0; n < nstrata; ++n) {
                // strata_idx[n] = (int *) S_alloc(strata_size[n], sizeof(int));
                strata_idx.add(new int[strata_size[n]]);
            }
            rfutils.zeroInt(strata_size, nstrata);
            for (n = 0; n < nsample0; ++n) {
                strata_size[strata[n] - 1]++;
                strata_idx.get(strata[n] - 1)[strata_size[strata[n] - 1]
                        - 1] = n;
            }
        } else

        {
            nind = (replace > 0) ? null : new int[nsample];// (int *)
                                                           // S_alloc(nsample,
                                                           // sizeof(int));
        }

        /* INITIALIZE FOR RUN */
        if (testdat > 0)
            rfutils.zeroDouble(countts, ntest * nclass);
        rfutils.zeroInt(counttr, nclass * nsample);
        rfutils.zeroInt(_out, nsample);
        rfutils.zeroDouble((ScootArray<Double>) tgini.copy(), mdim);
        rfutils.zeroDouble((ScootArray<Double>) errtr.copy(),
                (nclass + 1) * Ntree);

        if (labelts > 0) {
            nclts = new int[nclass]; // (int *) S_alloc(nclass, sizeof(int));
            for (n = 0; n < ntest; ++n)
                nclts[clts[n] - 1]++;
            rfutils.zeroDouble((ScootArray<Double>) errts.copy(),
                    (nclass + 1) * Ntree);
        }
        // #ifndef USER
        else {
            // if xtest != null => testdat = true
            // if ytest == null => labelts = false
            // in that situation, nclts was not being initialized.
            nclts = new int[1];// (int *) S_alloc(1, sizeof(int));
            nclts[0] = 0;
        }
        // #endif

        if (imp > 0) {
            rfutils.zeroDouble(imprt, (nclass + 2) * mdim);
            rfutils.zeroDouble(impsd, (nclass + 1) * mdim);
            if (localImp > 0)
                rfutils.zeroDouble(impmat, nsample * mdim);
        }
        if (iprox > 0) {
            rfutils.zeroDouble(prox, nsample0 * nsample0);
            if (testdat > 0)
                rfutils.zeroDouble(proxts, ntest * (ntest + nsample0));
        }
        rfutils.makeA(x, mdim, nsample, (ScootArray<Integer>) cat.copy(),
                at, (ScootArray<Integer>) b.copy());
        assert (b.pointer == 0);

        // R_CheckUserInterrupt();

        /* Starting the main loop over number of trees. */
        // GetRNGstate();
        if (trace <= Ntree) {
            /* Print header for running output. */
            // Rprintf("ntree OOB");
            // for (n = 1; n <= nclass; ++n) Rprintf("%7i", n);
            // if (*labelts) {
            // Rprintf("| Test");
            // for (n = 1; n <= nclass; ++n) Rprintf("%7i", n);
            // }
            // Rprintf("\n");
        }
        idxByNnode = 0;
        idxByNsample = 0;
        for (jb = 0; jb < Ntree; jb++) {
            /* Do we need to simulate data for the second class? */
            if (addClass > 0)
                rfutils.createClass(x, nsample0, nsample, mdim, oUniformRandom);
            do {
                rfutils.zeroInt(nodestatus.scootRight(idxByNnode), nrnodes);
                rfutils.zeroInt(treemap.scootRight(2 * idxByNnode),
                        2 * nrnodes);
                rfutils.zeroDouble(xbestsplit.scootRight(nuse), nrnodes);
                rfutils.zeroInt(nodeclass.scootRight(idxByNnode), nrnodes);
                rfutils.zeroInt(varUsed.copy(), mdim);
                /* TODO: Put all sampling code into a function. */
                /* drawSample(sampsize, nsample, ); */
                if (stratify > 0) { /* stratified sampling */
                    rfutils.zeroInt(jin, nsample);
                    rfutils.zeroDouble((ScootArray<Double>) tclasspop.copy(),
                            nclass);
                    rfutils.zeroDouble((ScootArray<Double>) win.copy(),
                            nsample);
                    if (replace > 0) { /* with replacement */
                        for (n = 0; n < nstrata; ++n) {
                            for (j = 0; j < sampsize[n]; ++j) {
                                // ktmp = (int) (unif_rand() * strata_size[n]);
                                ktmp = (int) (Math
                                        .floor(oUniformRandom.nextDouble()
                                                * (double) (strata_size[n])));
                                k = strata_idx.get(n)[ktmp];
                                tclasspop.set(cl.get(k) - 1,
                                        tclasspop.get(cl.get(k) - 1)
                                                + classwt[cl.get(k) - 1]);
                                win.set(k, win.get(k) + classwt[cl.get(k) - 1]);
                                jin[k] = 1;
                            }
                        }
                    } else { /* stratified sampling w/o replacement */
                        /* re-initialize the index array */
                        rfutils.zeroInt(strata_size, nstrata);
                        for (j = 0; j < nsample; ++j) {
                            strata_size[strata[j] - 1]++;
                            strata_idx.get(
                                    strata[j] - 1)[strata_size[strata[j] - 1]
                                            - 1] = j;
                        }
                        /* sampling without replacement */
                        for (n = 0; n < nstrata; ++n) {
                            last = strata_size[n] - 1;
                            for (j = 0; j < sampsize[n]; ++j) {
                                // ktmp = (int) (unif_rand() * (last+1));
                                ktmp = (int) (Math
                                        .floor(oUniformRandom.nextDouble()
                                                * (double) (last + 1)));
                                k = strata_idx.get(n)[ktmp];
                                swapInt(strata_idx.get(n)[last],
                                        strata_idx.get(n)[ktmp]);
                                last--;
                                tclasspop.set(cl.get(k) - 1,
                                        tclasspop.get(cl.get(k) - 1)
                                                + classwt[cl.get(k) - 1]);
                                win.set(k, win.get(k) + classwt[cl.get(k) - 1]);
                                jin[k] = 1;
                            }
                        }
                    }
                } else { /* unstratified sampling */
                    ntry = 0;
                    do {
                        nEmpty = 0;
                        rfutils.zeroInt(jin, nsample);
                        rfutils.zeroDouble(
                                (ScootArray<Double>) tclasspop.copy(), nclass);
                        rfutils.zeroDouble((ScootArray<Double>) win.copy(),
                                nsample);
                        if (replace > 0) {
                            for (n = 0; n < sampsize[0]; ++n) {
                                /*
                                 * #ifdef USER k = unif_rand() * nsample; #else
                                 * k = (int) (unif_rand() * (double) nsample);
                                 * #endif
                                 */
                                k = (int) (Math
                                        .floor(oUniformRandom.nextDouble()
                                                * (double) (nsample)));
                                tclasspop.set(cl.get(k) - 1,
                                        tclasspop.get(cl.get(k) - 1)
                                                + classwt[cl.get(k) - 1]);
                                win.set(k, win.get(k) + classwt[cl.get(k) - 1]);
                                jin[k] = 1;
                            }
                        } else {
                            for (n = 0; n < nsample; ++n)
                                nind[n] = n;
                            last = nsample - 1;
                            for (n = 0; n < sampsize[0]; ++n) {
                                // ktmp = (int) (unif_rand() * (last+1));
                                ktmp = (int) (Math
                                        .floor(oUniformRandom.nextDouble()
                                                * (double) (last + 1)));
                                k = nind[ktmp];
                                swapInt(nind[ktmp], nind[last]);
                                last--;
                                tclasspop.set(cl.get(k) - 1,
                                        tclasspop.get(cl.get(k) - 1)
                                                + classwt[cl.get(k) - 1]);
                                win.set(k, win.get(k) + classwt[cl.get(k) - 1]);
                                jin[k] = 1;
                            }
                        }
                        /* check if any class is missing in the sample */
                        for (n = 0; n < nclass; ++n) {
                            if (tclasspop.get(n) == 0.0)
                                nEmpty++;
                        }
                        ntry++;
                    } while (nclass - nEmpty < 2 && ntry <= 30);
                    /*
                     * If there are still fewer than two classes in the data,
                     * throw an error.
                     */
                    if (nclass - nEmpty < 2)
                        // error("Still have fewer than two classes in the
                        // in-bag sample after 30 attempts.");
                        throw new Exception(
                                "Still have fewer than two classes in the in-bag sample after 30 attempts.");
                }

                /* If need to keep indices of inbag data, do that here. */
                if (keepInbag > 0) {
                    for (n = 0; n < nsample0; ++n) {
                        inbag[n + idxByNsample] = jin[n];
                    }
                }

                /* Copy the original a matrix back. */
                // memcpy(a, at, sizeof(int) * mdim * nsample); // destination,
                // source, size
                assert (a.size() == at.length);
                for (int tempi = 0; tempi < a.size(); tempi++)
                    a.set(tempi, at[tempi]);
                rfutils.modA(a.copy(), nuse, nsample, mdim, cat.copy(),
                        maxcat, ncase.copy(), jin);
                assert (a.pointer == 0);

                rfsub.buildtree(a.copy(), b.copy(), cl.copy(), cat.copy(),
                        maxcat, mdim, nsample,
                        nclass,
                        treemap.scootRight(2 * idxByNnode),
                        bestvar.scootRight(idxByNnode),
                        bestsplit.copy(), bestsplitnext.copy(), tgini.copy(),
                        nodestatus.scootRight(idxByNnode), nodepop.copy(),
                        nodestart.copy(), classpop.copy(), tclasspop.copy(),
                        tclasscat.copy(),
                        ta.copy(), nrnodes, idmove.copy(), ndsize,
                        ncase.copy(),
                        mtry, varUsed.copy(), nodeclass.scootRight(idxByNnode),
                        ndbigtree.scootRight(jb), win.copy(), wr.copy(),
                        wl.copy(),
                        mdim,
                        nuse, mind.copy(), oUniformRandom); // fails here on jb
                                                            // >= 4
                assert (a.pointer == 0 && b.pointer == 0
                        && bestsplit.pointer == 0);
                assert (
                        win.pointer == 0 && wl.pointer == 0 && wl.pointer == 0);

                /* if the "tree" has only the root node, start over */
            } while (ndbigtree.get(jb) == 1);

            rfutils.Xtranslate(x, mdim, nrnodes, nsample,
                    bestvar.scootRight(idxByNnode),
                    bestsplit.copy(), bestsplitnext.copy(),
                    xbestsplit.scootRight(idxByNnode),
                    nodestatus.scootRight(idxByNnode), cat.copy(),
                    ndbigtree.get(jb));
            assert (bestsplit.pointer == 0);

            /* Get test set error */
            if (testdat > 0) {
                Tree.predictClassTree(xts, ntest, mdim,
                        treemap.scootRight(2 * idxByNnode),
                        nodestatus.scootRight(idxByNnode),
                        xbestsplit.scootRight(idxByNnode),
                        bestvar.scootRight(idxByNnode),
                        nodeclass.scootRight(idxByNnode), ndbigtree.get(jb),
                        cat.copy(), nclass, jts.copy(), nodexts.copy(),
                        maxcat);
                TestSetError(countts, jts.copy(), clts, outclts, ntest, nclass,
                        jb + 1,
                        errts.scootRight(jb * (nclass + 1)), labelts, nclts,
                        cut,
                        oUniformRandom);
            }

            /* Get out-of-bag predictions and errors. */
            Tree.predictClassTree(x, nsample, mdim,
                    treemap.scootRight(2 * idxByNnode),
                    nodestatus.scootRight(idxByNnode),
                    xbestsplit.scootRight(idxByNnode),
                    bestvar.scootRight(idxByNnode),
                    nodeclass.scootRight(idxByNnode), ndbigtree.get(jb),
                    cat.copy(), nclass, jtr.copy(), nodex.copy(), maxcat);

            rfutils.zeroInt(nout, nclass);
            noutall = 0;
            for (n = 0; n < nsample; ++n) {
                if (jin[n] == 0) {
                    /* increment the OOB votes */
                    counttr[n * nclass + jtr.get(n) - 1]++;
                    /* count number of times a case is OOB */
                    _out[n]++;
                    /*
                     * count number of OOB cases in the current iteration.
                     * nout[n] is the number of OOB cases for the n-th class.
                     * noutall is the number of OOB cases overall.
                     */
                    nout[cl.get(n) - 1]++;
                    noutall++;
                }
            }

            /* Compute out-of-bag error rate. */
            oob(nsample, nclass, jin, cl.copy(), jtr.copy(), jerr, counttr,
                    _out,
                    errtr.scootRight(jb * (nclass + 1)), outcl, cut,
                    oUniformRandom);
            assert (jtr.pointer == 0);

            if ((jb + 1) % trace == 0) {
                /*
                 * Rprintf("%5i: %6.2f%%", jb+1, 100.0*errtr[jb * (nclass+1)]);
                 * for (n = 1; n <= nclass; ++n) { Rprintf("%6.2f%%", 100.0 *
                 * errtr[n + jb * (nclass+1)]); } if (*labelts) { Rprintf("| ");
                 * for (n = 0; n <= nclass; ++n) { Rprintf("%6.2f%%", 100.0 *
                 * errts[n + jb * (nclass+1)]); } } Rprintf("\n"); #ifdef WIN32
                 * R_FlushConsole(); R_ProcessEvents(); #endif
                 * R_CheckUserInterrupt();
                 */
            }

            /* DO PROXIMITIES */
            if (iprox > 0) {
                rfutils.computeProximity(prox, oobprox,
                        (ScootArray<Integer>) nodex.copy(), jin, oobpair, near);
                /* proximity for test data */
                if (testdat > 0) {
                    rfutils.computeProximity(proxts, 0,
                            (ScootArray<Integer>) nodexts.copy(), jin, oobpair,
                            ntest);
                    /* Compute proximity between testset and training set. */
                    for (n = 0; n < ntest; ++n) {
                        for (k = 0; k < near; ++k) {
                            if (nodexts.get(n) == nodex.get(k))
                                proxts[n + ntest * (k + ntest)] += 1.0;
                        }
                    }
                }
            }

            /* DO VARIABLE IMPORTANCE */
            if (imp > 0) {
                nrightall = 0;
                /*
                 * Count the number of correct prediction by the current tree
                 * among the OOB samples, by class.
                 */
                rfutils.zeroInt(nright, nclass);
                for (n = 0; n < nsample; ++n) {
                    /* out-of-bag and predicted correctly: */
                    if (jin[n] == 0 && jtr.get(n) == cl.get(n)) {
                        nright[cl.get(n) - 1]++;
                        nrightall++;
                    }
                }
                for (m = 0; m < mdim; ++m) {
                    if (varUsed.get(m) > 0) {
                        nrightimpall = 0;
                        rfutils.zeroInt(nrightimp, nclass);
                        for (n = 0; n < nsample; ++n)
                            tx[n] = x[m + n * mdim];
                        /* Permute the m-th variable. */
                        rfutils.permuteOOB(m, x, jin, nsample, mdim,
                                oUniformRandom);
                        /* Predict the modified data using the current tree. */
                        Tree.predictClassTree(x, nsample, mdim,
                                treemap.scootRight(2 * idxByNnode),
                                nodestatus.scootRight(idxByNnode),
                                xbestsplit.scootRight(idxByNnode),
                                bestvar.scootRight(idxByNnode),
                                nodeclass.scootRight(idxByNnode),
                                ndbigtree.get(jb),
                                cat.copy(), nclass, jvr.copy(), nodex.copy(),
                                maxcat);

                        /*
                         * Count how often correct predictions are made with the
                         * modified data.
                         */
                        for (n = 0; n < nsample; n++) {
                            /* Restore the original data for that variable. */
                            x[m + n * mdim] = tx[n];
                            if (jin[n] == 0) {
                                if (jvr.get(n) == cl.get(n)) {
                                    nrightimp[cl.get(n) - 1]++;
                                    nrightimpall++;
                                }
                                if (localImp > 0 && jvr.get(n) != jtr.get(n)) {
                                    if (cl.get(n) == jvr.get(n)) {
                                        impmat[m + n * mdim] -= 1.0;
                                    } else {
                                        impmat[m + n * mdim] += 1.0;
                                    }
                                }
                            }
                        }
                        /*
                         * Accumulate decrease in proportions of correct
                         * predictions.
                         */
                        for (n = 0; n < nclass; ++n) {
                            if (nout[n] > 0) {
                                imprt[m + n * mdim] += ((double) (nright[n]
                                        - nrightimp[n])) /
                                        nout[n];
                                impsd[m + n * mdim] += ((double) (nright[n]
                                        - nrightimp[n]) *
                                        (nright[n] - nrightimp[n])) / nout[n];
                            }
                        }
                        if (noutall > 0) {
                            imprt[m + nclass * mdim] += ((double) (nrightall
                                    - nrightimpall)) / noutall;
                            impsd[m + nclass * mdim] += ((double) (nrightall
                                    - nrightimpall) *
                                    (nrightall - nrightimpall)) / noutall;
                        }
                    }
                }
            }

            /*
             * R_CheckUserInterrupt(); #ifdef WIN32 R_ProcessEvents(); #endif
             */
            if (keepf > 0)
                idxByNnode += nrnodes;
            if (keepInbag > 0)
                idxByNsample += nsample0;
        }
        // PutRNGstate();

        /* Final processing of variable importance. */
        for (m = 0; m < mdim; m++)
            tgini.set(m, tgini.get(m) / Ntree);
        if (imp > 0) {
            for (m = 0; m < mdim; ++m) {
                if (localImp > 0) { /* casewise measures */
                    for (n = 0; n < nsample; ++n)
                        impmat[m + n * mdim] /= _out[n];
                }
                /* class-specific measures */
                for (k = 0; k < nclass; ++k) {
                    av = imprt[m + k * mdim] / Ntree;
                    impsd[m + k * mdim] = Math
                            .sqrt(((impsd[m + k * mdim] / Ntree) - av * av)
                                    / (double) (Ntree));
                    imprt[m + k * mdim] = av;
                    /*
                     * imprt[m + k*mdim] = (se <= 0.0) ? -1000.0 - av : av / se;
                     */
                }
                /* overall measures */
                av = imprt[m + nclass * mdim] / Ntree;
                impsd[m + nclass * mdim] = Math
                        .sqrt(((impsd[m + nclass * mdim] / Ntree) - av * av)
                                / (double) (Ntree));
                imprt[m + nclass * mdim] = av;
                imprt[m + (nclass + 1) * mdim] = tgini.get(m);
            }
        } else {
            for (m = 0; m < mdim; ++m)
                imprt[m] = tgini.get(m);
        }

        /* PROXIMITY DATA ++++++++++++++++++++++++++++++++ */
        if (iprox > 0) {
            for (n = 0; n < near; ++n) {
                for (k = n + 1; k < near; ++k) {
                    double denom = (oobprox > 0)
                            ? ((oobpair[near * k + n] > 0)
                                    ? oobpair[near * k + n]
                                    : 1)
                            : Ntree;
                    prox[near * k + n] /= denom;
                    prox[near * n + k] = prox[near * k + n];
                }
                prox[near * n + n] = 1.0;
            }
            if (testdat > 0) {
                for (n = 0; n < ntest; ++n) {
                    for (k = 0; k < ntest + nsample; ++k)
                        proxts[ntest * k + n] /= Ntree;
                    proxts[ntest * n + n] = 1.0;
                }
            }
        }
    }

    @SuppressWarnings("unused")
    private static void classForest(int mdim, int ntest, int nclass, int maxcat,
            int nrnodes, int ntree, double[] x, ScootArray<Double> xbestsplit,
            double[] pid, double[] cutoff, double[] countts,
            ScootArray<Integer> treemap,
            ScootArray<Integer> nodestatus, ScootArray<Integer> cat,
            ScootArray<Integer> nodeclass,
            ScootArray<Integer> jts,
            int[] jet, ScootArray<Integer> bestvar, ScootArray<Integer> node,
            int[] treeSize,
            int keepPred, int[] prox, double[] proxMat, int[] nodes,
            Random oRandom) {
        int j, n, n1, n2, idxNodes, offset1, offset2, ntie;
        int[] junk;
        double crit, cmax;

        rfutils.zeroDouble(countts, nclass * ntest);
        idxNodes = 0;
        offset1 = 0;
        offset2 = 0;
        junk = null;

        for (j = 0; j < ntree; ++j) {
            /* predict by the j-th tree */
            Tree.predictClassTree(x, ntest, mdim,
                    treemap.scootRight(2 * idxNodes),
                    nodestatus.scootRight(idxNodes),
                    xbestsplit.scootRight(idxNodes),
                    bestvar.scootRight(idxNodes),
                    nodeclass.scootRight(idxNodes),
                    treeSize[j], cat.copy(), nclass,
                    jts.scootRight(offset1), node.scootRight(offset2), maxcat);
            /* accumulate votes: */
            for (n = 0; n < ntest; ++n) {
                countts[jts.get(n + offset1) - 1 + n * nclass] += 1.0;
            }

            /* if desired, do proximities for this round */
            if (prox[0] > 0)
                rfutils.computeProximity(proxMat, 0, node.scootRight(offset2),
                        junk, junk,
                        ntest);
            idxNodes += nrnodes;
            if (keepPred > 0)
                offset1 += ntest;
            if (nodes[0] > 0)
                offset2 += ntest;
        }

        /* Aggregated prediction is the class with the maximum votes/cutoff */
        for (n = 0; n < ntest; ++n)

        {
            cmax = 0.0;
            ntie = 1;
            for (j = 0; j < nclass; ++j) {
                crit = (countts[j + n * nclass] / ntree) / cutoff[j];
                if (crit > cmax) {
                    jet[n] = j + 1;
                    cmax = crit;
                    ntie = 1;
                }
                /* Break ties at random: */
                if (crit == cmax) {
                    // if (unif_rand() < 1.0 / ntie)
                    if (oRandom.nextDouble() < 1.0 / (double) (ntie))
                        jet[n] = j + 1;
                    ntie++;
                }
            }
        }

        /*
         * if proximities requested, do the final adjustment (division by number
         * of trees)
         */
        if (prox[0] > 0) {
            for (n1 = 0; n1 < ntest; ++n1) {
                for (n2 = n1 + 1; n2 < ntest; ++n2) {
                    proxMat[n1 + n2 * ntest] /= ntree;
                    proxMat[n2 + n1 * ntest] = proxMat[n1 + n2 * ntest];
                }
                proxMat[n1 + n1 * ntest] = 1.0;
            }
        }
    }

    /*
     * Modified by A. Liaw 1/10/2003 (Deal with cutoff) Re-written in C by A.
     * Liaw 3/08/2004
     */
    private static void oob(int nsample, int nclass, int[] jin,
            ScootArray<Integer> cl, ScootArray<Integer> jtr, int[] jerr,
            int[] counttr, int[] _out, ScootArray<Double> errtr, int[] jest,
            double[] cutoff, Random oRandom) {

        int j, n, noob, ntie;
        int[] noobcl;
        double qq, smax, smaxtr;

        noobcl = new int[nclass]; // (int *) S_alloc(nclass, sizeof(int));
        rfutils.zeroInt(jerr, nsample);
        rfutils.zeroDouble((ScootArray<Double>) errtr.copy(), nclass + 1);

        noob = 0;
        for (n = 0; n < nsample; ++n) {
            if (_out[n] > 0) {
                noob++;
                noobcl[cl.get(n) - 1]++;
                smax = 0.0;
                smaxtr = 0.0;
                ntie = 1;
                for (j = 0; j < nclass; ++j) {
                    qq = (((double) counttr[j + n * nclass]) / _out[n])
                            / cutoff[j];
                    if (j + 1 != cl.get(n))
                        smax = (qq > smax) ? qq : smax;
                    /*
                     * if vote / cutoff is larger than current max, re-set max
                     * and change predicted class to the current class
                     */
                    if (qq > smaxtr) {
                        smaxtr = qq;
                        jest[n] = j + 1;
                        ntie = 1;
                    }
                    /* break tie at random */
                    if (qq == smaxtr) {
                        // if (unif_rand() < 1.0 / ntie) {
                        if (oRandom.nextDouble() < 1.0
                                / (double) (ntie)) {
                            smaxtr = qq;
                            jest[n] = j + 1;
                        }
                        ntie++;
                    }
                }
                if (jest[n] != cl.get(n)) {
                    errtr.set(cl.get(n), cl.get(n) + 1.0);
                    errtr.set(0, errtr.get(0) + 1.0);
                    jerr[n] = 1;
                }
            }
        }
        errtr.set(0, errtr.get(0) / noob);
        for (n = 1; n <= nclass; ++n)
            errtr.set(n, errtr.get(n) / noobcl[n - 1]);
    }

}