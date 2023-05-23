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
 * classTree.cs is based on classTree.c from the R package randomForest v4.6-6 licensed under the GPL v2 or later
 * modifications by Robert Carnell, August 2012
 */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace CsRandomForest
{
    /// <summary>
    /// 
    /// </summary>
    public class classTree
    {
        private static int NODE_TERMINAL = -1;
        //private static int NODE_TOSPLIT = -2; // not used
        //private static int NODE_INTERIOR = -3; // not used

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
        public static void predictClassTree(double[] x, int n, int mdim, CArray<int> treemap,
		              CArray<int> nodestatus, CArray<double> xbestsplit,
		              CArray<int> bestvar, CArray<int> nodeclass,
		              int treeSize, CArray<int> cat, int nclass,
		              CArray<int> jts, CArray<int> nodex, int maxcat) {
            int m, i, j, k;
            int[] cbestsplit = null;
	        uint npack;

            /* decode the categorical splits */
            if (maxcat > 1) {
                cbestsplit = new int[maxcat * treeSize]; //(int *) Calloc(maxcat * treeSize, int);
                rfutils.zeroInt(ref cbestsplit, maxcat * treeSize);
                for (i = 0; i < treeSize; ++i) {
                    if (nodestatus[i] != NODE_TERMINAL) {
                        if (cat[bestvar[i] - 1] > 1) {
                            npack = Convert.ToUInt32(xbestsplit[i]);
                            /* unpack `npack' into bits */
                            for (j = 0; npack > 0; npack >>= 1, ++j) {
                                cbestsplit[j + i*maxcat] = Convert.ToInt32(npack & 01);
                            }
                        }
                    }
                }
            }
            for (i = 0; i < n; ++i) {
		        k = 0;
		        while (nodestatus[k] != NODE_TERMINAL) {
                    m = bestvar[k] - 1;
                    if (cat[m] == 1) {
				        /* Split by a numerical predictor */
				        k = (x[m + i * mdim] <= xbestsplit[k]) ?
					        treemap[k * 2] - 1 : treemap[1 + k * 2] - 1;
			        } else {
				        /* Split by a categorical predictor */
				        k = (cbestsplit[(int) x[m + i * mdim] - 1 + k * maxcat] > 0)?
					        treemap[k * 2] - 1 : treemap[1 + k * 2] - 1;
			        }
		        }
		        /* Terminal node: assign class label */
		        jts[i] = nodeclass[k];
		        nodex[i] = k + 1;
            }
            //if (maxcat > 1) 
            //    Free(cbestsplit);
        }

        /// <summary>
        /// This finds the best split of a categorical variable with <code>lcat</code> categories 
        /// and <code>nclass</code> classes The method uses an exhaustive search over all partitions 
        /// of the category values if the number of categories is 10 or fewer.  Otherwise ncsplit randomly 
        /// selected splits are tested and best used.
        /// </summary>
        /// <param name="parentDen"></param>
        /// <param name="tclasscat"><code>tclasscat(j,k) is the number of cases in class <code>j</code> with category value <code>k</code></param>
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
        public static void catmax(ref double parentDen, CArray<double> tclasscat,
                              CArray<double> tclasspop, ref int nclass, ref int lcat,
                              ref int ncatsp, ref double critmax, ref int nhit,
                              ref int maxcat, ref int ncmax, ref int ncsplit,
                              Random oRandom) {
            int j, k, n, nsplit;
            int[] icat = new int[32];
            double leftNum, leftDen, rightNum, decGini;
            double[] leftCatClassCount = null;

            leftCatClassCount = new double[nclass]; //(double *) Calloc(*nclass, double);
            nhit = 0;
            nsplit = (lcat > ncmax) ?
                ncsplit : Convert.ToInt32(Math.Pow(2.0, Convert.ToDouble(lcat - 1)) - 1.0);

            for (n = 0; n < nsplit; ++n) {
                rfutils.zeroInt(ref icat, 32);
                if (lcat > ncmax) {
                    /* Generate random split.
                       TODO: consider changing to generating random bits with more
                       efficient algorithm */
                    for (j = 0; j < lcat; ++j) 
                        //icat[j] = unif_rand() > 0.5 ? 1 : 0;
                        icat[j] = oRandom.NextDouble() > 0.5 ? 1 : 0;
                } else {
                    rfutils.unpack(lcat, Convert.ToUInt32(n + 1), ref icat);
                }
                for (j = 0; j < nclass; ++j) {
                    leftCatClassCount[j] = 0;
                    for (k = 0; k < lcat; ++k) {
                        if (icat[k] > 0) {
                            leftCatClassCount[j] += tclasscat[j + k * nclass];
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
                    leftCatClassCount[j] = tclasspop[j] - leftCatClassCount[j];
                    rightNum += leftCatClassCount[j] * leftCatClassCount[j];
                }
                decGini = (leftNum / leftDen) + (rightNum / (parentDen - leftDen));
                if (decGini > critmax) {
                    critmax = decGini;
                    nhit = 1;
                    ncatsp = (lcat > ncmax) ? Convert.ToInt32(rfutils.pack(lcat, ref icat)) : n + 1;
                }
            }
            //Free(leftCatClassCount);
        }



        /// <summary>
        /// Find the best split of with categorical variable when there are two classes
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
        public static void catmaxb(ref double totalWt, CArray<double> tclasscat, 
            CArray<double> classCount, 
                               ref int nclass, ref int nCat, ref int nbest, ref double critmax,
                               ref int nhit, ref double[] catCount) {

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
                    tclasscat[i * nclass] / catCount[i] : 0.0;
                kcat[i] = i + 1;
            }
            //R_qsort_I(catProportion, kcat, 1, nCat);
            Array.Sort(catProportion, kcat); // sort the items (kcat) by the keys (catProportion) and the keys (catProportion)
            for (i = 0; i < nclass; ++i)
            {
                cp[i] = 0;
                cm[i] = classCount[i];
            }
            rightDen = totalWt;
            leftDen = 0.0;
            for (i = 0; i < nCat - 1; ++i) {
                leftDen += catCount[kcat[i]-1];
                rightDen -= catCount[kcat[i]-1];
                leftNum = 0.0;
                rightNum = 0.0;
                for (j = 0; j < nclass; ++j) {
                    cp[j] += tclasscat[j + (kcat[i]-1) * nclass];
                    cm[j] -= tclasscat[j + (kcat[i]-1) * nclass];
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
                rfutils.zeroInt(ref kcat, nCat);
                for (i = 0; i < nCat; ++i) {
                    catProportion[i] = (catCount[i] > 0) ?
                        tclasscat[i * nclass] / catCount[i] : 0.0;
                    kcat[i] = catProportion[i] < bestsplit ? 1 : 0;
			        /* Rprintf("%i ", kcat[i]); */
                }
                nbest = Convert.ToInt32(rfutils.pack(nCat, ref kcat));
		        /* Rprintf("\nnbest=%u\nnbest=%i\n", *nbest, *nbest); */
            }
        }

    }
}
