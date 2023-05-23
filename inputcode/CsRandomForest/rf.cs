/*
 *  Copyright (C) 2012 Rob Carnell
 *  
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *  
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA,
 * or visit http://www.gnu.org/licenses/gpl-2.0.html
 * 
 * rf.cs is based on rf.c from the R package randomForest v4.6-6 licensed under the GPL v2 or later
 * modifications by Robert Carnell, August 2012
 */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Diagnostics;

namespace CsRandomForest
{
    public class rf
    {
        private static void swapInt(ref int a, ref int b)
        {
            int oldb = b;
            b = a;
            a = oldb;
        }

        private static void TestSetError(double[] countts, CArray<int> jts, int[] clts, 
                  int[] jet, int ntest,
		          int nclass, int nvote, CArray<double> errts,
		          int labelts, int[] nclts, double[] cutoff, Random oRandom) {
            int j, n, ntie;
            double cmax, crit;

            for (n = 0; n < ntest; ++n) 
                countts[jts[n]-1 + n*nclass] += 1.0;

            /*  Prediction is the class with the maximum votes */
            for (n = 0; n < ntest; ++n) {
		        cmax=0.0;
		        ntie = 1;
		        for (j = 0; j < nclass; ++j) {
			        crit = (countts[j + n*nclass] / nvote) / cutoff[j];
			        if (crit > cmax) {
				        jet[n] = j+1;
				        cmax = crit;
				        ntie = 1;
			        }
			        /*  Break ties at random: */
			        if (crit == cmax) {
				        //if (unif_rand() < 1.0 / ntie) {
                        if (oRandom.NextDouble() < 1.0 / Convert.ToDouble(ntie)) {
					        jet[n] = j+1;
					        cmax = crit;
				        }
				        ntie++;
			        }
		        }
            }
            if (labelts > 0) {
                rfutils.zeroDouble(errts.copy(), nclass + 1);
		        for (n = 0; n < ntest; ++n) {
			        if (jet[n] != clts[n]) {
				        errts[0] += 1.0;
				        errts[clts[n]] += 1.0;
			        }
		        }
		        errts[0] /= ntest;
		        for (n = 1; n <= nclass; ++n) errts[n] /= nclts[n-1];
            }
        }

        /// <summary>
        /// C# wrapper for random forests
        /// </summary>
        /// <param name="x">"matrix of predictors (transposed!)" - variables are rows, observations are columns</param>
        /// <param name="dimx">two integers: number of variables and number of cases</param>
        /// <param name="cl">class labels of the data</param>
        /// <param name="ncl">number of classes in the response</param>
        /// <param name="cat">integer vector of number of classes in the predictor: 1=continuous</param>
        /// <param name="maxcat">maximum of cat</param>
        /// <param name="sampsize">Sizes of samples to draw.  For classification, if sampsize is a vector of the length 
        /// of the number of strata, then sampling is stratified by strata, and teh elements of sampsize indicate the numbers 
        /// to be drawn from the strata</param>
        /// <param name="strata">A variable that is used for stratified sampling</param>
        /// <param name="Options">
        ///     <list type="num">
        ///         <listheader>7 integers: (0=no, 1=yes)</listheader>
        ///         <item>add a second class (for unsupervised RF)?</item>
        ///         <list type="bull">
        ///             <item>1: sampling from product of marginals</item>
        ///             <item>2: sampling from product of uniforms</item>
        ///         </list>
        ///         <item>assess variable importance?</item>
        ///         <item>calculate proximity?</item>
        ///         <item>calculate proximity based on OOB predictions?</item>
        ///         <item>calculate outlying measure?</item>
        ///         <item>how often to print output?</item>
        ///         <item>keep the forest for future prediction?</item>
        ///     </list>
        /// </param>
        /// <param name="ntree">number of trees</param>
        /// <param name="nvar">number of predictors to use for each split</param>
        /// <param name="ipi">0=use class proportion as prob.; 1=use supplied priors</param>
        /// <param name="classwt">Priors of the classes.  Need not add up to one.  Ignored for regression</param>
        /// <param name="cut">(Classification Only) A vector of length equal to the number of classes.  The 'winning' 
        /// class for an observation is the one with teh maximum ratio of proportion of votes to cutoff.  Default is 
        /// 1/k where k is the number of classes (i.e., majority vote wins)</param>
        /// <param name="nodesize">minimum node size: no node with fewer than ndsize cases will be split</param>
        /// <param name="outcl">Output:  class predicted by RF</param>
        /// <param name="counttr">Output:  matrix of votes (transposed!)</param>
        /// <param name="prox">Output:  matrix of proximity (if iprox=1)</param>
        /// <param name="imprt">Output:  matrix of variable importance measures.  A matrix with <code>nclass+2</code> 
        /// (for classification) or two (for regression) columns.  For classification, the first <code>nclass</code> 
        /// columns are the class-specific measures computed as mean decrease in accuracy.  The <code>nclass+1st</code> 
        /// columns is the mean decrease in accuracy over all classes.  Teh last column is the mean decrease in Gini index.  
        /// For regression, the first column is the mean decrease in accuracy and the second column is the mean decrease in 
        /// MSE.  if <code>importance=false</code>, the last measure is still returned as a vector</param>
        /// <param name="impsd">Output:  The "standard errors" of the permutation-based importance measure.  For classification, 
        /// a <code>p x (nclass + 1)</code> matrix corresponding to the first <code>nclass+1</code> columns of the importance 
        /// matrix.  For regression, a lenght <code>p</code> vector</param>
        /// <param name="impmat">Output:  matrix of local variable importance measures.  A <code>p x n</code> matrix containing
        /// the casewise importance measures, the <code>[i,j]</code> element of which is the importance of the <code>i-th</code>
        /// variable on the <code>j-th</code> case.  NULL if <code>localImp=false</code></param>
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
        /// <param name="oUniformRandom">The random number generator used by the random forest</param>
        public static void classRF(double[] x, int[] dimx, CArray<int> cl, ref int ncl, 
                 CArray<int> cat, ref int maxcat,
	             int[] sampsize, int[] strata, int[] Options, ref int ntree, ref int nvar,
	             ref int ipi, double[] classwt, double[] cut, ref int nodesize,
	             int[] outcl, int[] counttr, double[] prox,
	             double[] imprt, double[] impsd, double[] impmat, ref int nrnodes,
	             CArray<int> ndbigtree, CArray<int> nodestatus, CArray<int> bestvar, 
                 CArray<int> treemap,
	             CArray<int> nodeclass, CArray<double> xbestsplit, CArray<double> errtr,
	             ref int testdat, double[] xts, int[] clts, ref int nts, double[] countts,
	             int[] outclts, ref int labelts, double[] proxts, CArray<double> errts,
                 int[] inbag, Random oUniformRandom) 
        {
            int nsample0, mdim, nclass, addClass, mtry, ntest, nsample, ndsize,
                mimp, nimp, near, nuse, noutall, nrightall, nrightimpall,
	        keepInbag, nstrata;
            int jb, j, n, m, k, idxByNnode, idxByNsample, imp, localImp, iprox,
	        oobprox, keepf, replace, stratify, trace;
            int[] nright, nrightimp, nout, nclts;
            int Ntree;

            nstrata = 0;
            nuse = 0;

            int[] _out, jin, jerr, classFreq, at, nind, oobpair;
            oobpair = null;
            List<int[]> strata_idx = null;
            int[] strata_size = null;
            int last, ktmp, nEmpty, ntry;

            double av=0.0;

            double[] tx, tp;

            int[] a_mem, b_mem, bestsplit_mem, bestsplitnext_mem, nodepop_mem,
                nodestart_mem, ta_mem, idmove_mem, ncase_mem, varUsed_mem,
                mind_mem, nodex_mem, nodexts_mem, jts_mem, jtr_mem, jvr_mem;
            CArray<int> a, b, bestsplit, bestsplitnext, nodepop, nodestart, ta, idmove,
                ncase, varUsed, mind, nodex, nodexts, jts, jtr, jvr;
            double[] tgini_mem, classpop_mem, tclasspop_mem, tclasscat_mem, win_mem, wr_mem, wl_mem;
            CArray<double> tgini, classpop, tclasspop, tclasscat, win, wr, wl;

            if (Options.Length != 10)
                throw new ArgumentException("error message 1");
            addClass = Options[0];
            imp      = Options[1];
            localImp = Options[2];
            iprox    = Options[3];
            oobprox  = Options[4];
            trace    = Options[5];
            keepf    = Options[6];
            replace  = Options[7];
            stratify = Options[8];
            keepInbag = Options[9];
            if (dimx.Length != 2)
                throw new ArgumentException("error message 2");
            mdim     = dimx[0];
            nsample0 = dimx[1];
            nclass   = (ncl == 1) ? 2 : ncl;
            ndsize   = nodesize;
            Ntree    = ntree;
            mtry     = nvar;
            ntest    = nts;
            nsample = (addClass > 0) ? (nsample0 + nsample0) : nsample0;
            mimp = (imp > 0) ? mdim : 1;
            nimp = (imp > 0) ? nsample : 1;
            near = (iprox > 0) ? nsample0 : 1;
            if (trace == 0) trace = Ntree + 1;

            tgini_mem = new double[mdim];//     (double *] S_alloc(mdim, sizeof(double]];
            wl_mem = new double[nclass]; //(double *] S_alloc(nclass, sizeof(double]];
            wr_mem = new double[nclass]; //(double *] S_alloc(nclass, sizeof(double]];
            classpop_mem = new double[nclass * nrnodes];//(double *] S_alloc(nclass* *nrnodes, sizeof(double]];
            tclasscat_mem = new double[nclass * 32];//(double *] S_alloc(nclass*32, sizeof(double]];
            tclasspop_mem = new double[nclass];//(double *] S_alloc(nclass, sizeof(double]];
            win_mem = new double[nsample];//(double *] S_alloc(nsample, sizeof(double]];
            tgini = new CArray<double>(tgini_mem);//     (double *) S_alloc(mdim, sizeof(double));
            wl = new CArray<double>(wl_mem); //(double *) S_alloc(nclass, sizeof(double));
            wr = new CArray<double>(wr_mem); //(double *) S_alloc(nclass, sizeof(double));
            classpop = new CArray<double>(classpop_mem);//(double *) S_alloc(nclass* *nrnodes, sizeof(double));
            tclasscat = new CArray<double>(tclasscat_mem);//(double *) S_alloc(nclass*32, sizeof(double));
            tclasspop = new CArray<double>(tclasspop_mem);//(double *) S_alloc(nclass, sizeof(double));
            win = new CArray<double>(win_mem);//(double *) S_alloc(nsample, sizeof(double));
            tx = new double[nsample];//(double *) S_alloc(nsample, sizeof(double));
            tp = new double[nsample];//(double *) S_alloc(nsample, sizeof(double));

            bestsplitnext_mem = new int[nrnodes];//(int *] S_alloc(*nrnodes, sizeof(int]];
            bestsplit_mem = new int[nrnodes];//(int *] S_alloc(*nrnodes, sizeof(int]];
            nodepop_mem = new int[nrnodes];//(int *] S_alloc(*nrnodes, sizeof(int]];
            nodestart_mem = new int[nrnodes];//(int *] S_alloc(*nrnodes, sizeof(int]];
            nodex_mem = new int[nsample];//(int *] S_alloc(nsample, sizeof(int]];
            nodexts_mem = new int[ntest];//(int *] S_alloc(ntest, sizeof(int]];
            ta_mem = new int[nsample];//(int *] S_alloc(nsample, sizeof(int]];
            ncase_mem = new int[nsample];//(int *] S_alloc(nsample, sizeof(int]];
            varUsed_mem = new int[mdim];//(int *] S_alloc(mdim, sizeof(int]];
            jtr_mem = new int[nsample];//(int *] S_alloc(nsample, sizeof(int]];
            jvr_mem = new int[nsample];//(int *] S_alloc(nsample, sizeof(int]];
            jts_mem = new int[ntest];//(int *] S_alloc(ntest, sizeof(int]];
            idmove_mem = new int[nsample];//(int *] S_alloc(nsample, sizeof(int]];
            a_mem = new int[mdim * nsample];//(int *] S_alloc(mdim*nsample, sizeof(int]];
            b_mem = new int[mdim * nsample];//(int *] S_alloc(mdim*nsample, sizeof(int]];
            mind_mem = new int[mdim];//(int *] S_alloc(mdim, sizeof(int]];

            bestsplitnext = new CArray<int>(bestsplitnext_mem);//(int *) S_alloc(*nrnodes, sizeof(int));
            bestsplit = new CArray<int>(bestsplit_mem);//(int *) S_alloc(*nrnodes, sizeof(int));
            nodepop = new CArray<int>(nodepop_mem);//(int *) S_alloc(*nrnodes, sizeof(int));
            nodestart = new CArray<int>(nodestart_mem);//(int *) S_alloc(*nrnodes, sizeof(int));
            nodex = new CArray<int>(nodex_mem);//(int *) S_alloc(nsample, sizeof(int));
            nodexts = new CArray<int>(nodexts_mem);//(int *) S_alloc(ntest, sizeof(int));
            ta = new CArray<int>(ta_mem);//(int *) S_alloc(nsample, sizeof(int));
            ncase = new CArray<int>(ncase_mem);//(int *) S_alloc(nsample, sizeof(int));
            varUsed = new CArray<int>(varUsed_mem);//(int *) S_alloc(mdim, sizeof(int));
            jtr = new CArray<int>(jtr_mem);//(int *) S_alloc(nsample, sizeof(int));
            jvr = new CArray<int>(jvr_mem);//(int *) S_alloc(nsample, sizeof(int));
            jts = new CArray<int>(jts_mem);//(int *) S_alloc(ntest, sizeof(int));
            idmove = new CArray<int>(idmove_mem);//(int *) S_alloc(nsample, sizeof(int));
            a = new CArray<int>(a_mem);//(int *) S_alloc(mdim*nsample, sizeof(int));
            b = new CArray<int>(b_mem);//(int *) S_alloc(mdim*nsample, sizeof(int));
            mind = new CArray<int>(mind_mem);//(int *) S_alloc(mdim, sizeof(int));
            
            jin = new int[nsample];//(int *) S_alloc(nsample, sizeof(int));
            _out = new int[nsample];//(int *) S_alloc(nsample, sizeof(int));
            jerr = new int[nsample];//(int *) S_alloc(nsample, sizeof(int));
            classFreq = new int[nclass];//(int *) S_alloc(nclass, sizeof(int));
            at = new int[mdim * nsample];//(int *) S_alloc(mdim*nsample, sizeof(int));
            nright = new int[nclass];//(int *) S_alloc(nclass, sizeof(int));
            nrightimp = new int[nclass];//(int *) S_alloc(nclass, sizeof(int));
            nout = new int[nclass];//(int *) S_alloc(nclass, sizeof(int));

            if (oobprox > 0) {
	            oobpair = new int[near*near]; //(int *) S_alloc(near*near, sizeof(int));
            }

            /* Count number of cases in each class. */
            rfutils.zeroInt(ref classFreq, nclass);
            for (n = 0; n < nsample; ++n) classFreq[cl[n] - 1] ++;
            /* Normalize class weights. */
            rfutils.normClassWt(cl.copy(), nsample, nclass, ipi, ref classwt, ref classFreq);

            nind = null;
            if (stratify > 0) {
	            /* Count number of strata and frequency of each stratum. */
	            nstrata = 0;
	            for (n = 0; n < nsample0; ++n)
	                if (strata[n] > nstrata) nstrata = strata[n];
                    /* Create the array of pointers, each pointing to a vector
	                of indices of where data of each stratum is. */
                        strata_size = new int[nstrata]; //(int  *) S_alloc(nstrata, sizeof(int));
	            for (n = 0; n < nsample0; ++n) {
	                strata_size[strata[n] - 1] ++;
	            }
	            strata_idx =  new List<int[]>(nstrata); //(int **) S_alloc(nstrata, sizeof(int *));
	            for (n = 0; n < nstrata; ++n) {
	                //strata_idx[n] = (int *) S_alloc(strata_size[n], sizeof(int));
                    strata_idx.Add(new int[strata_size[n]]);
	            }
	            rfutils.zeroInt(ref strata_size, nstrata);
	            for (n = 0; n < nsample0; ++n) {
	                strata_size[strata[n] - 1] ++;
	                strata_idx[strata[n] - 1][strata_size[strata[n] - 1] - 1] = n;
	            }
            } else {
	            nind = (replace > 0) ? null : new int[nsample];//(int *) S_alloc(nsample, sizeof(int));
            }

            /*    INITIALIZE FOR RUN */
            if (testdat > 0) rfutils.zeroDouble(ref countts, ntest * nclass);
            rfutils.zeroInt(ref counttr, nclass * nsample);
            rfutils.zeroInt(ref _out, nsample);
            rfutils.zeroDouble(tgini.copy(), mdim);
            rfutils.zeroDouble(errtr.copy(), (nclass + 1) * Ntree);

            if (labelts > 0) {
	            nclts  = new int[nclass]; //(int *) S_alloc(nclass, sizeof(int));
	            for (n = 0; n < ntest; ++n) nclts[clts[n]-1]++;
	                rfutils.zeroDouble(errts.copy(), (nclass + 1) * Ntree);
            }
        //#ifndef USER
	        else
	        {
		        // if xtest != null => testdat = true
		        // if ytest == null => labelts = false
		        // in that situation, nclts was not being initialized.
		        nclts = new int[1];//(int *) S_alloc(1, sizeof(int));
		        nclts[0] = 0;
	        }
        //#endif

            if (imp > 0) {
                rfutils.zeroDouble(ref imprt, (nclass+2) * mdim);
                rfutils.zeroDouble(ref impsd, (nclass+1) * mdim);
	            if (localImp > 0) 
                    rfutils.zeroDouble(ref impmat, nsample * mdim);
            }
            if (iprox > 0) {
                rfutils.zeroDouble(ref prox, nsample0 * nsample0);
                if (testdat > 0) 
                    rfutils.zeroDouble(ref proxts, ntest * (ntest + nsample0));
            }
            rfutils.makeA(ref x, mdim, nsample, cat.copy(), ref at, b.copy());
            Debug.Assert(b.Offset == 0);

            //R_CheckUserInterrupt();

            /* Starting the main loop over number of trees. */
            //GetRNGstate();
            if (trace <= Ntree) {
	            /* Print header for running output. */
	            //Rprintf("ntree      OOB");
	            //for (n = 1; n <= nclass; ++n) Rprintf("%7i", n);
	            //if (*labelts) {
	                //Rprintf("|    Test");
	                //for (n = 1; n <= nclass; ++n) Rprintf("%7i", n);
	            //}
	            //Rprintf("\n");
            }
            idxByNnode = 0;
            idxByNsample = 0;
            for (jb = 0; jb < Ntree; jb++) {
                /* Do we need to simulate data for the second class? */
                if (addClass > 0) 
                    rfutils.createClass(ref x, nsample0, nsample, mdim, oUniformRandom);
		        do {
			        rfutils.zeroInt(nodestatus + idxByNnode, nrnodes);
			        rfutils.zeroInt(treemap + 2*idxByNnode, 2 * nrnodes);
			        rfutils.zeroDouble(xbestsplit + idxByNnode, nrnodes);
			        rfutils.zeroInt(nodeclass + idxByNnode, nrnodes);
                    rfutils.zeroInt(varUsed.copy(), mdim);
                    /* TODO: Put all sampling code into a function. */
                    /* drawSample(sampsize, nsample, ); */
			        if (stratify > 0) {  /* stratified sampling */
				        rfutils.zeroInt(ref jin, nsample);
				        rfutils.zeroDouble(tclasspop.copy(), nclass);
				        rfutils.zeroDouble(win.copy(), nsample);
				        if (replace > 0) {  /* with replacement */
					        for (n = 0; n < nstrata; ++n) {
						        for (j = 0; j < sampsize[n]; ++j) {
							        //ktmp = (int) (unif_rand() * strata_size[n]);
                                    ktmp = Convert.ToInt32(Math.Floor(oUniformRandom.NextDouble() * Convert.ToDouble(strata_size[n])));
							        k = strata_idx[n][ktmp];
							        tclasspop[cl[k] - 1] += classwt[cl[k] - 1];
							        win[k] += classwt[cl[k] - 1];
							        jin[k] = 1;
						        }
					        }
				        } else { /* stratified sampling w/o replacement */
					        /* re-initialize the index array */
					        rfutils.zeroInt(ref strata_size, nstrata);
					        for (j = 0; j < nsample; ++j) {
						        strata_size[strata[j] - 1] ++;
						        strata_idx[strata[j] - 1][strata_size[strata[j] - 1] - 1] = j;
					        }
					        /* sampling without replacement */
					        for (n = 0; n < nstrata; ++n) {
						        last = strata_size[n] - 1;
						        for (j = 0; j < sampsize[n]; ++j) {
							        //ktmp = (int) (unif_rand() * (last+1));
                                    ktmp = Convert.ToInt32(Math.Floor(oUniformRandom.NextDouble() * Convert.ToDouble(last+1)));
							        k = strata_idx[n][ktmp];
                                    swapInt(ref strata_idx[n][last], ref strata_idx[n][ktmp]);
							        last--;
							        tclasspop[cl[k] - 1] += classwt[cl[k]-1];
							        win[k] += classwt[cl[k]-1];
							        jin[k] = 1;
						        }
					        }
				        }
			        } else {  /* unstratified sampling */
				        ntry = 0;
				        do {
					        nEmpty = 0;
					        rfutils.zeroInt(ref jin, nsample);
					        rfutils.zeroDouble(tclasspop.copy(), nclass);
					        rfutils.zeroDouble(win.copy(), nsample);
					        if (replace > 0) {
						        for (n = 0; n < sampsize[0]; ++n) {
        /*#ifdef USER
							        k = unif_rand() * nsample;
        #else
							        k = (int) (unif_rand() * (double) nsample);
        #endif*/
                                    k = Convert.ToInt32(Math.Floor(oUniformRandom.NextDouble() * Convert.ToDouble(nsample)));
							        tclasspop[cl[k] - 1] += classwt[cl[k] - 1];
							        win[k] += classwt[cl[k]-1];
							        jin[k] = 1;
						        }
					        } else {
						        for (n = 0; n < nsample; ++n) nind[n] = n;
						        last = nsample - 1;
						        for (n = 0; n < sampsize[0]; ++n) {
							        //ktmp = (int) (unif_rand() * (last+1));
                                    ktmp = Convert.ToInt32(Math.Floor(oUniformRandom.NextDouble() * Convert.ToDouble(last + 1)));
							        k = nind[ktmp];
                                    swapInt(ref nind[ktmp], ref nind[last]);
							        last--;
							        tclasspop[cl[k] - 1] += classwt[cl[k]-1];
							        win[k] += classwt[cl[k]-1];
							        jin[k] = 1;
						        }
					        }
					        /* check if any class is missing in the sample */
					        for (n = 0; n < nclass; ++n) {
						        if (tclasspop[n] == 0.0) nEmpty++;
					        }
					        ntry++;
				        } while (nclass - nEmpty < 2 && ntry <= 30);
				        /* If there are still fewer than two classes in the data, throw an error. */
				        if (nclass - nEmpty < 2) 
                            //error("Still have fewer than two classes in the in-bag sample after 30 attempts.");
                            throw new Exception("Still have fewer than two classes in the in-bag sample after 30 attempts.");
			        }

                    /* If need to keep indices of inbag data, do that here. */
                    if (keepInbag > 0) {
                        for (n = 0; n < nsample0; ++n) {
                            inbag[n + idxByNsample] = jin[n];
                        }
                    }

			        /* Copy the original a matrix back. */
			        //memcpy(a, at, sizeof(int) * mdim * nsample); // destination, source, size
                    Debug.Assert(a.Length == at.Length);
                    for (int tempi = 0; tempi < a.Length; tempi++)
                        a[tempi] = at[tempi];
      	            rfutils.modA(a.copy(), ref nuse, nsample, mdim, cat.copy(), ref maxcat, ncase.copy(), ref jin);
                    Debug.Assert(a.Offset == 0);

			        rfsub.buildtree(a.copy(), b.copy(), cl.copy(), cat.copy(), ref maxcat, ref mdim, ref nsample,
								        ref nclass,
								        treemap + 2*idxByNnode, bestvar + idxByNnode,
								        bestsplit.copy(), bestsplitnext.copy(), tgini.copy(),
								        nodestatus + idxByNnode, nodepop.copy(),
								        nodestart.copy(), classpop.copy(), tclasspop.copy(), tclasscat.copy(),
								        ta.copy(), ref nrnodes, idmove.copy(), ref ndsize, ncase.copy(),
								        ref mtry, varUsed.copy(), nodeclass + idxByNnode,
								        ndbigtree + jb, win.copy(), wr.copy(), wl.copy(), ref mdim,
								        ref nuse, mind.copy(), oUniformRandom); // fails here on jb >= 4
                    Debug.Assert(a.Offset == 0 && b.Offset == 0 && bestsplit.Offset == 0);
                    Debug.Assert(win.Offset == 0 && wl.Offset == 0 && wl.Offset == 0);

			        /* if the "tree" has only the root node, start over */
		        } while (ndbigtree[jb] == 1);

		        rfutils.Xtranslate(x, mdim, nrnodes, nsample, bestvar + idxByNnode,
				   bestsplit.copy(), bestsplitnext.copy(), xbestsplit + idxByNnode,
				   nodestatus + idxByNnode, cat.copy(), ndbigtree[jb]);
                Debug.Assert(bestsplit.Offset == 0);

		        /*  Get test set error */
		        if (testdat > 0) {
                    classTree.predictClassTree(xts, ntest, mdim, treemap + 2*idxByNnode,
                                     nodestatus + idxByNnode, xbestsplit + idxByNnode,
                                     bestvar + idxByNnode,
                                     nodeclass + idxByNnode, ndbigtree[jb],
                                     cat.copy(), nclass, jts.copy(), nodexts.copy(), maxcat);
			        TestSetError(countts, jts.copy(), clts, outclts, ntest, nclass, jb+1,
						         errts + jb*(nclass+1), labelts, nclts, cut, oUniformRandom);
                }

		        /*  Get out-of-bag predictions and errors. */
                classTree.predictClassTree(x, nsample, mdim, treemap + 2*idxByNnode,
                                 nodestatus + idxByNnode, xbestsplit + idxByNnode,
                                 bestvar + idxByNnode,
                                 nodeclass + idxByNnode, ndbigtree[jb],
                                 cat.copy(), nclass, jtr.copy(), nodex.copy(), maxcat);
                
		        rfutils.zeroInt(ref nout, nclass);
		        noutall = 0;
		        for (n = 0; n < nsample; ++n) {
			        if (jin[n] == 0) {
				        /* increment the OOB votes */
				        counttr[n*nclass + jtr[n] - 1] ++;
				        /* count number of times a case is OOB */
				        _out[n]++;
				        /* count number of OOB cases in the current iteration.
				           nout[n] is the number of OOB cases for the n-th class.
				           noutall is the number of OOB cases overall. */
				        nout[cl[n] - 1]++;
				        noutall++;
			        }
		        }

                /* Compute out-of-bag error rate. */
		        oob(nsample, nclass, ref jin, cl.copy(), jtr.copy(), ref jerr, ref counttr, ref _out,
			        errtr + jb*(nclass+1), ref outcl, ref cut, oUniformRandom);
                Debug.Assert(jtr.Offset == 0);

		        if ((jb+1) % trace == 0) {
			        /*Rprintf("%5i: %6.2f%%", jb+1, 100.0*errtr[jb * (nclass+1)]);
			        for (n = 1; n <= nclass; ++n) {
				        Rprintf("%6.2f%%", 100.0 * errtr[n + jb * (nclass+1)]);
			        }
			        if (*labelts) {
				        Rprintf("| ");
				        for (n = 0; n <= nclass; ++n) {
					        Rprintf("%6.2f%%", 100.0 * errts[n + jb * (nclass+1)]);
				        }
			        }
			        Rprintf("\n");
        #ifdef WIN32
			        R_FlushConsole();
			        R_ProcessEvents();
        #endif
			        R_CheckUserInterrupt();*/
		        }

		        /*  DO PROXIMITIES */
		        if (iprox > 0) {
                    rfutils.computeProximity(ref prox, oobprox, nodex.copy(), ref jin, ref oobpair, near);
			        /* proximity for test data */
			        if (testdat > 0) {
                        rfutils.computeProximity(ref proxts, 0, nodexts.copy(), ref jin, ref oobpair, ntest);
                        /* Compute proximity between testset and training set. */
				        for (n = 0; n < ntest; ++n) {
					        for (k = 0; k < near; ++k) {
						        if (nodexts[n] == nodex[k])
							        proxts[n + ntest * (k+ntest)] += 1.0;
					        }
				        }
			        }
		        }

		        /*  DO VARIABLE IMPORTANCE  */
		        if (imp > 0) {
			        nrightall = 0;
			        /* Count the number of correct prediction by the current tree
			           among the OOB samples, by class. */
			        rfutils.zeroInt(ref nright, nclass);
			        for (n = 0; n < nsample; ++n) {
       	                /* out-of-bag and predicted correctly: */
				        if (jin[n] == 0 && jtr[n] == cl[n]) {
					        nright[cl[n] - 1]++;
					        nrightall++;
				        }
			        }
			        for (m = 0; m < mdim; ++m) {
				        if (varUsed[m] > 0) {
					        nrightimpall = 0;
					        rfutils.zeroInt(ref nrightimp, nclass);
					        for (n = 0; n < nsample; ++n) tx[n] = x[m + n*mdim];
					        /* Permute the m-th variable. */
                            rfutils.permuteOOB(m, ref x, ref jin, nsample, mdim, oUniformRandom);
					        /* Predict the modified data using the current tree. */
                            classTree.predictClassTree(x, nsample, mdim, treemap + 2*idxByNnode,
                                             nodestatus + idxByNnode,
                                             xbestsplit + idxByNnode,
                                             bestvar + idxByNnode,
                                             nodeclass + idxByNnode, ndbigtree[jb],
                                             cat.copy(), nclass, jvr.copy(), nodex.copy(), maxcat);

                            /* Count how often correct predictions are made with
					           the modified data. */
					        for (n = 0; n < nsample; n++) {
						        /* Restore the original data for that variable. */
						        x[m + n*mdim] = tx[n];
						        if (jin[n] == 0) {
							        if (jvr[n] == cl[n]) {
								        nrightimp[cl[n] - 1]++;
								        nrightimpall++;
							        }
							        if (localImp > 0 && jvr[n] != jtr[n]) {
								        if (cl[n] == jvr[n]) {
									        impmat[m + n*mdim] -= 1.0;
								        } else {
									        impmat[m + n*mdim] += 1.0;
								        }
							        }
						        }
					        }
					        /* Accumulate decrease in proportions of correct
					           predictions. */
					        for (n = 0; n < nclass; ++n) {
						        if (nout[n] > 0) {
							        imprt[m + n*mdim] +=
								        ((double) (nright[n] - nrightimp[n])) /
								        nout[n];
							        impsd[m + n*mdim] +=
								        ((double) (nright[n] - nrightimp[n]) *
								         (nright[n] - nrightimp[n])) / nout[n];
						        }
					        }
					        if (noutall > 0) {
						        imprt[m + nclass*mdim] +=
							        ((double)(nrightall - nrightimpall)) / noutall;
						        impsd[m + nclass*mdim] +=
							        ((double) (nrightall - nrightimpall) *
							         (nrightall - nrightimpall)) / noutall;
					        }
				        }
			        }
		        }

		        /*R_CheckUserInterrupt();
        #ifdef WIN32
		        R_ProcessEvents();
        #endif*/
                if (keepf > 0) idxByNnode += nrnodes;
                if (keepInbag > 0) idxByNsample += nsample0;
            }
            //PutRNGstate();

            /*  Final processing of variable importance. */
            for (m = 0; m < mdim; m++) tgini[m] /= Ntree;
            if (imp > 0) {
		        for (m = 0; m < mdim; ++m) {
			        if (localImp > 0) { /* casewise measures */
				        for (n = 0; n < nsample; ++n) impmat[m + n*mdim] /= _out[n];
			        }
			        /* class-specific measures */
			        for (k = 0; k < nclass; ++k) {
				        av = imprt[m + k*mdim] / Ntree;
				        impsd[m + k*mdim] =
                            Math.Sqrt(((impsd[m + k*mdim] / Ntree) - av*av) / Convert.ToDouble(Ntree));
				        imprt[m + k*mdim] = av;
				        /* imprt[m + k*mdim] = (se <= 0.0) ? -1000.0 - av : av / se; */
			        }
			        /* overall measures */
			        av = imprt[m + nclass*mdim] / Ntree;
			        impsd[m + nclass*mdim] =
                        Math.Sqrt(((impsd[m + nclass*mdim] / Ntree) - av*av) / Convert.ToDouble(Ntree));
			        imprt[m + nclass*mdim] = av;
			        imprt[m + (nclass+1)*mdim] = tgini[m];
		        }
            } else {
		        for (m = 0; m < mdim; ++m) imprt[m] = tgini[m];
            }

            /*  PROXIMITY DATA ++++++++++++++++++++++++++++++++*/
            if (iprox > 0) {
		        for (n = 0; n < near; ++n) {
			        for (k = n + 1; k < near; ++k) {
                        double denom = (oobprox > 0) ? ((oobpair[near * k + n] > 0) ? oobpair[near * k + n] : 1) : Ntree;
                        prox[near * k + n] /= denom;
                        prox[near * n + k] = prox[near * k + n];
			        }
			        prox[near*n + n] = 1.0;
		        }
		        if (testdat > 0) {
			        for (n = 0; n < ntest; ++n) {
				        for (k = 0; k < ntest + nsample; ++k)
					        proxts[ntest*k + n] /= Ntree;
				        proxts[ntest * n + n] = 1.0;
			        }
		        }
            }
        }

        private static void classForest(ref int mdim, ref int ntest, ref int nclass, ref int maxcat,
                         ref int nrnodes, ref int ntree, ref double[] x, CArray<double> xbestsplit,
                         ref double[] pid, ref double[] cutoff, ref double[] countts, CArray<int> treemap,
                         CArray<int> nodestatus, CArray<int> cat, CArray<int> nodeclass, 
                         CArray<int> jts,
                         ref int[] jet, CArray<int> bestvar, CArray<int> node, ref int[] treeSize,
                         ref int keepPred, ref int[] prox, ref double[] proxMat, int[] nodes,
                         Random oRandom) {
            int j, n, n1, n2, idxNodes, offset1, offset2, ntie;
            int[] junk;
            double crit, cmax;

            rfutils.zeroDouble(ref countts, nclass * ntest);
            idxNodes = 0;
            offset1 = 0;
            offset2 = 0;
            junk = null;

            for (j = 0; j < ntree; ++j) {
		        /* predict by the j-th tree */
                classTree.predictClassTree(x, ntest, mdim, treemap + 2*idxNodes,
			         nodestatus + idxNodes, xbestsplit + idxNodes,
			         bestvar + idxNodes, nodeclass + idxNodes,
			         treeSize[j], cat.copy(), nclass,
			         jts + offset1, node + offset2, maxcat);
		        /* accumulate votes: */
		        for (n = 0; n < ntest; ++n) {
			        countts[jts[n + offset1] - 1 + n * nclass] += 1.0;
		        }

		        /* if desired, do proximities for this round */
		        if (prox[0] > 0) 
                    rfutils.computeProximity(ref proxMat, 0, node + offset2, ref junk, ref junk, ntest);
		        idxNodes += nrnodes;
		        if (keepPred > 0) offset1 += ntest;
		        if (nodes[0] > 0)    offset2 += ntest;
            }

            /* Aggregated prediction is the class with the maximum votes/cutoff */
            for (n = 0; n < ntest; ++n) {
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
				        //if (unif_rand() < 1.0 / ntie) 
                        if (oRandom.NextDouble() < 1.0 / Convert.ToDouble(ntie))
                            jet[n] = j + 1;
				        ntie++;
			        }
		        }
            }

            /* if proximities requested, do the final adjustment
               (division by number of trees) */
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
          Modified by A. Liaw 1/10/2003 (Deal with cutoff)
          Re-written in C by A. Liaw 3/08/2004
        */
        private static void oob(int nsample, int nclass, ref int[] jin, 
            CArray<int> cl, CArray<int> jtr, ref int[] jerr,
	        ref int[] counttr, ref int[] _out, CArray<double> errtr, ref int[] jest,
	        ref double[] cutoff, Random oRandom) {

            int j, n, noob, ntie;
            int[] noobcl;
            double qq, smax, smaxtr;

            noobcl  = new int[nclass]; //(int *) S_alloc(nclass, sizeof(int));
            rfutils.zeroInt(ref jerr, nsample);
            rfutils.zeroDouble(errtr.copy(), nclass+1);

            noob = 0;
            for (n = 0; n < nsample; ++n) {
	        if (_out[n] > 0) {
	            noob++;
	            noobcl[cl[n]-1]++;
	            smax = 0.0;
	            smaxtr = 0.0;
	            ntie = 1;
	            for (j = 0; j < nclass; ++j) {
	    	        qq = (((double) counttr[j + n*nclass]) / _out[n]) / cutoff[j];
	    	        if (j+1 != cl[n]) smax = (qq > smax) ? qq : smax;
	    	        /* if vote / cutoff is larger than current max, re-set max and
		   	           change predicted class to the current class */
	    	        if (qq > smaxtr) {
	    		        smaxtr = qq;
	    		        jest[n] = j+1;
	    		        ntie = 1;
	    	        }
	    	        /* break tie at random */
	    	        if (qq == smaxtr) {
	    		        //if (unif_rand() < 1.0 / ntie) {
                        if (oRandom.NextDouble() < 1.0 / Convert.ToDouble(ntie)) {
	    			        smaxtr = qq;
	    			        jest[n] = j+1;
	    		        }
	    		        ntie++;
	    	        }
	            }
	            if (jest[n] != cl[n]) {
		        errtr[cl[n]] += 1.0;
		        errtr[0] += 1.0;
		        jerr[n] = 1;
	            }
	        }
            }
            errtr[0] /= noob;
            for (n = 1; n <= nclass; ++n) errtr[n] /= noobcl[n-1];
        }

    }
}
