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
 * rfsub.cs is based on rfsub.f from the R package randomForest v4.6-6 licensed under the GPL v2 or later
 * modifications by Robert Carnell, August 2012
 */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace CsRandomForest
{
    /// <summary>
    /// Random Forest subroutines
    /// </summary>
    public class rfsub
    {
        /* Table of constant values */

        private static int c__32 = 32;

        private static int zerv_(CArray<int> ix, ref int m1)
        {
            /* System generated locals */
            int i__1;

            /* Local variables */
            int n;

            /* Parameter adjustments */
            --ix;

            /* Function Body */
            i__1 = m1;
            for (n = 1; n <= i__1; ++n) {
	            ix[n] = 0;
            }
            return 0;
        } /* zerv_ */

        /* Subroutine */ 
        private static int zermr_(CArray<double> rx, ref int m1, ref int m2)
        {
            /* System generated locals */
            int rx_dim1, rx_offset, i__1, i__2;

            /* Local variables */
            int i__, j;

            /* Parameter adjustments */
            rx_dim1 = m1;
            rx_offset = 1 + rx_dim1;
            rx -= rx_offset;

            /* Function Body */
            i__1 = m1;
            for (i__ = 1; i__ <= i__1; ++i__) {
	        i__2 = m2;
	        for (j = 1; j <= i__2; ++j) {
	            rx[i__ + j * rx_dim1] = 0.0;
	        }
            }
            return 0;
        } /* zermr_ */

        private static int zervr_(CArray<double> rx, ref int m1)
        {
            /* System generated locals */
            int i__1;

            /* Local variables */
            int n;

            /* Parameter adjustments */
            --rx;

            /* Function Body */
            i__1 = m1;
            for (n = 1; n <= i__1; ++n) {
	        rx[n] = 0.0;
            }
            return 0;
        } /* zervr_ */

        private static int zerm_(CArray<int> mx, ref int m1, ref int m2)
        {
            /* System generated locals */
            int mx_dim1, mx_offset, i__1, i__2;

            /* Local variables */
            int i__, j;

            /* Parameter adjustments */
            mx_dim1 = m1;
            mx_offset = 1 + mx_dim1;
            mx -= mx_offset;

            /* Function Body */
            i__1 = m1;
            for (i__ = 1; i__ <= i__1; ++i__) {
	        i__2 = m2;
	        for (j = 1; j <= i__2; ++j) {
	            mx[i__ + j * mx_dim1] = 0;
	        }
            }
            return 0;
        } /* zerm_ */

        private static int zermd_(CArray<double> rx, ref int m1, ref int m2)
        {
            /* System generated locals */
            int rx_dim1, rx_offset, i__1, i__2;

            /* Local variables */
            int i__, j;

            /* Parameter adjustments */
            rx_dim1 = m1;
            rx_offset = 1 + rx_dim1;
            rx -= rx_offset;

            /* Function Body */
            i__1 = m1;
            for (i__ = 1; i__ <= i__1; ++i__) {
	        i__2 = m2;
	        for (j = 1; j <= i__2; ++j) {
	            rx[i__ + j * rx_dim1] = 0.0;
	        }
            }
            return 0;
        } /* zermd_ */

        /// <summary>
        /// Buildtree consists of repeated calls to two subroutines, Findbestsplit and Movedata.  Findbestsplit does 
        /// just that--it finds the best split of the current node.  Movedata moves the data in the split node 
        /// right and left so that the data corresponding to each child node is contiguous. The buildtree bookkeeping 
        /// is different from that in Friedman's original CART program.  <code>ncur</code> is the total number of nodes to date. 
        /// <code>nodestatus(k)=1</code> if the <code>kth</code> node has been split.  <code>nodestatus(k)=2</code> if the node exists but has not yet been 
        /// split, and <code>=-1</code> of the node is terminal. A node is terminal if its size is below a threshold value, or 
        /// if it is all one class, or if all the x-values are equal.  If the current node <code>k</code> is split, then its 
        /// children are numbered <code>ncur+1</code> (left), and <code>ncur+2</code> (right), <code>ncur</code> increases to <code>ncur+2</code> and the next node 
        /// to be split is numbered <code>k+1</code>.  When no more nodes can be split, buildtree returns to the main program.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="cl"></param>
        /// <param name="cat"></param>
        /// <param name="maxcat"></param>
        /// <param name="mdim"></param>
        /// <param name="nsample"></param>
        /// <param name="nclass"></param>
        /// <param name="treemap"></param>
        /// <param name="bestvar"></param>
        /// <param name="bestsplit"></param>
        /// <param name="bestsplitnext"></param>
        /// <param name="tgini"></param>
        /// <param name="nodestatus"></param>
        /// <param name="nodepop"></param>
        /// <param name="nodestart"></param>
        /// <param name="classpop"></param>
        /// <param name="tclasspop"></param>
        /// <param name="tclasscat"></param>
        /// <param name="ta"></param>
        /// <param name="nrnodes"></param>
        /// <param name="idmove"></param>
        /// <param name="ndsize"></param>
        /// <param name="ncase"></param>
        /// <param name="mtry"></param>
        /// <param name="iv"></param>
        /// <param name="nodeclass"></param>
        /// <param name="ndbigtree"></param>
        /// <param name="win"></param>
        /// <param name="wr"></param>
        /// <param name="wl"></param>
        /// <param name="mred"></param>
        /// <param name="nuse"></param>
        /// <param name="mind"></param>
        /// <param name="oRandom"></param>
        /// <returns></returns>
        public static int buildtree(CArray<int> a, CArray<int> b, CArray<int> cl, 
            CArray<int> cat, ref int maxcat, ref int mdim, ref int nsample, 
            ref int nclass, 
            CArray<int> treemap, CArray<int> bestvar, 
            CArray<int> bestsplit, CArray<int> bestsplitnext, CArray<double> tgini, 
            CArray<int> nodestatus, CArray<int> nodepop, 
            CArray<int> nodestart, CArray<double> classpop, CArray<double> tclasspop, 
            CArray<double> tclasscat, 
            CArray<int> ta, ref int nrnodes, CArray<int> idmove, ref int ndsize, CArray<int> ncase, 
            ref int mtry, CArray<int> iv, CArray<int> nodeclass, 
            CArray<int> ndbigtree, CArray<double> win, CArray<double> wr, CArray<double> wl, 
            ref int mred, 
            ref int nuse, CArray<int> mind,
            Random oRandom)
        {
            /* System generated locals */
            int a_dim1, a_offset, b_dim1, b_offset, classpop_dim1, 
	            classpop_offset, tclasscat_dim1, tclasscat_offset, i__1, i__2;

            /* Local variables */
            double decsplit = 0.0;
            int j, k, n, nc, kn;
            double pp;
            int ntie, ncur;
            double popt1, popt2;
            int ndend;
            int nbest = 0;
            double xrand;
            int jstat;
            int ndendl = 0;
            int kbuild, msplit, ndstart;

            /* Parameter adjustments */
            --tgini;
            --cat;
            --win;
            --ncase;
            --idmove;
            --ta;
            --cl;
            b_dim1 = mdim;
            b_offset = 1 + b_dim1;
            b -= b_offset;
            a_dim1 = mdim;
            a_offset = 1 + a_dim1;
            a -= a_offset;
            --wl;
            --wr;
            tclasscat_dim1 = nclass;
            tclasscat_offset = 1 + tclasscat_dim1;
            tclasscat -= tclasscat_offset;
            --tclasspop;
            --nodeclass;
            classpop_dim1 = nclass;
            classpop_offset = 1 + classpop_dim1;
            classpop -= classpop_offset;
            --nodestart;
            --nodepop;
            --nodestatus;
            --bestsplitnext;
            --bestsplit;
            --bestvar;
            treemap -= 3;
            --mind;
            --iv;

            /* Function Body */
            msplit = 0;
            //zerv_(&nodestatus[1], nrnodes);
            zerv_(nodestatus + 1, ref nrnodes);
            zerv_(nodestart + 1, ref nrnodes);
            zerv_(nodepop + 1, ref nrnodes);
            zermr_(classpop + classpop_offset, ref nclass, ref nrnodes);

            i__1 = nclass;
            for (j = 1; j <= i__1; ++j) {
		        classpop[j + classpop_dim1] = tclasspop[j];
            }
            ncur = 1;
            nodestart[1] = 1;
            nodepop[1] = nuse;
            nodestatus[1] = 2;
        /*     start main loop */
            i__1 = nrnodes;
            for (kbuild = 1; kbuild <= i__1; ++kbuild) {
        /*         call intpr("kbuild", 6, kbuild, 1) */
        /*         call intpr("ncur", 4, ncur, 1) */
		        if (kbuild > ncur) {
			        goto L50;
		        }
		        if (nodestatus[kbuild] != 2) {
			        goto L30;
		        }
	        /*     initialize for next call to findbestsplit */
		        ndstart = nodestart[kbuild];
		        ndend = ndstart + nodepop[kbuild] - 1;
		        i__2 = nclass;
		        for (j = 1; j <= i__2; ++j) {
			        tclasspop[j] = classpop[j + kbuild * classpop_dim1];
		        }
		        jstat = 0;
		        findbestsplit(a + a_offset, b + b_offset, cl + 1, ref mdim, ref nsample, 
			        ref nclass, cat + 1, ref maxcat, ref ndstart, ref ndend, tclasspop + 1, 
                    tclasscat + tclasscat_offset, ref msplit, ref decsplit, ref nbest, 
                    ncase + 1, ref jstat, ref mtry, win + 1, wr + 1, wl + 1, ref mred, mind + 1,
                    oRandom);
	        /*         call intpr("jstat", 5, jstat, 1) */
	        /*         call intpr("msplit", 6, msplit, 1) */
	        /*     If the node is terminal, move on.  Otherwise, split. */
		        if (jstat == -1) {
			        nodestatus[kbuild] = -1;
			        goto L30;
		        } else {
			        bestvar[kbuild] = msplit;
			        iv[msplit] = 1;
			        if (decsplit < 0.0) {
			        decsplit = 0.0;
			        }
			        tgini[msplit] += decsplit;
			        if (cat[msplit] == 1) {
			        bestsplit[kbuild] = a[msplit + nbest * a_dim1];
			        bestsplitnext[kbuild] = a[msplit + (nbest + 1) * a_dim1];
			        } else {
			        bestsplit[kbuild] = nbest;
			        bestsplitnext[kbuild] = 0;
			        }
		        }
                movedata(a + a_offset, ta + 1, ref mdim, ref nsample, ref ndstart, ref ndend,
                    idmove + 1, ncase + 1, ref msplit, cat + 1, ref nbest, ref ndendl);

	        /*         call intpr("ndend", 5, ndend, 1) */
	        /*         call intpr("ndendl", 6, ndendl, 1) */
	        /*     leftnode no.= ncur+1, rightnode no. = ncur+2. */
		        nodepop[ncur + 1] = ndendl - ndstart + 1;
		        nodepop[ncur + 2] = ndend - ndendl;
		        nodestart[ncur + 1] = ndstart;
		        nodestart[ncur + 2] = ndendl + 1;
	        /*     find class populations in both nodes */
		        i__2 = ndendl;
		        for (n = ndstart; n <= i__2; ++n) {
			        nc = ncase[n];
			        j = cl[nc];
			        classpop[j + (ncur + 1) * classpop_dim1] += win[nc];
		        }
		        i__2 = ndend;
		        for (n = ndendl + 1; n <= i__2; ++n) {
			        nc = ncase[n];
			        j = cl[nc];
			        classpop[j + (ncur + 2) * classpop_dim1] += win[nc];
		        }
	        /*         call intpr("nL", 2, nodepop(ncur+1), 1) */
	        /*         call intpr("nR", 2, nodepop(ncur+2), 1) */
	        /*     check on nodestatus */
		        nodestatus[ncur + 1] = 2;
		        nodestatus[ncur + 2] = 2;
		        if (nodepop[ncur + 1] <= ndsize) {
			        nodestatus[ncur + 1] = -1;
		        }
		        if (nodepop[ncur + 2] <= ndsize) {
			        nodestatus[ncur + 2] = -1;
		        }
		        popt1 = 0.0;
		        popt2 = 0.0;
		        i__2 = nclass;
		        for (j = 1; j <= i__2; ++j) {
			        popt1 += classpop[j + (ncur + 1) * classpop_dim1];
			        popt2 += classpop[j + (ncur + 2) * classpop_dim1];
		        }
		        i__2 = nclass;
		        for (j = 1; j <= i__2; ++j) {
			        if (classpop[j + (ncur + 1) * classpop_dim1] == popt1) {
			        nodestatus[ncur + 1] = -1;
			        }
			        if (classpop[j + (ncur + 2) * classpop_dim1] == popt2) {
			        nodestatus[ncur + 2] = -1;
			        }
		        }
		        treemap[(kbuild << 1) + 1] = ncur + 1;
		        treemap[(kbuild << 1) + 2] = ncur + 2;
		        nodestatus[kbuild] = 1;
		        ncur += 2;
		        if (ncur >= nrnodes) {
			        goto L50;
		        }
	        L30:
		        ;
	        }
	        L50:
		        ndbigtree[0] = nrnodes;
		        for (k = nrnodes; k >= 1; --k) {
		        if (nodestatus[k] == 0) {
			        --(ndbigtree[0]);
		        }
		        if (nodestatus[k] == 2) {
			        nodestatus[k] = -1;
		        }
	        }
        /*     form prediction in terminal nodes */
            i__1 = ndbigtree[0];
            for (kn = 1; kn <= i__1; ++kn) {
		        if (nodestatus[kn] == -1) {
			        pp = 0.0;
			        ntie = 1;
			        i__2 = nclass;
			        for (j = 1; j <= i__2; ++j) {
				        if (classpop[j + kn * classpop_dim1] > pp) {
					        nodeclass[kn] = j;
					        pp = classpop[j + kn * classpop_dim1];
					        ntie = 1;
				        }
		        /*     Break ties at random: */
				        if (classpop[j + kn * classpop_dim1] == pp) {
                            xrand = oRandom.NextDouble();
					        //rrand_(&xrand);
					        if (xrand < 1.0 / Convert.ToDouble(ntie)) {
						        nodeclass[kn] = j;
						        pp = classpop[j + kn * classpop_dim1];
					        }
					        ++ntie;
				        }
			        }
		        }
        /*         call intpr("node", 4, kn, 1) */
        /*         call intpr("status", 6, nodestatus(kn), 1) */
        /*         call intpr("pred", 4, nodeclass(kn), 1) */
        /*         call dblepr("pop1", 4, classpop(1, kn), 1) */
        /*         call dblepr("pop2", 4, classpop(2, kn), 1) */
            }
            return 0;
        } /* buildtree_ */

        /*     SUBROUTINE FINDBESTSPLIT */
        /*     For the best split, msplit is the variable split on. decsplit is the */
        /*     dec. in impurity.  If msplit is numerical, nsplit is the case number */
        /*     of value of msplit split on, and nsplitnext is the case number of the */
        /*     next larger value of msplit.  If msplit is categorical, then nsplit is */
        /*     the coding into an integer of the categories going left. */
        private static int findbestsplit(CArray<int> a, CArray<int> b, CArray<int> cl, 
	        ref int mdim, ref int nsample, ref int nclass, CArray<int> cat, 
	        ref int maxcat, ref int ndstart, ref int ndend, CArray<double> tclasspop, 
            CArray<double> tclasscat, ref int msplit, ref double decsplit, 
            ref int nbest, CArray<int> ncase, ref int jstat, ref int mtry, 
            CArray<double> win, CArray<double> wr, CArray<double> wl, ref int mred, 
	        CArray<int> mind, Random oRandom)
        {
            /* System generated locals */
            int a_dim1, a_offset, b_dim1, b_offset, tclasscat_dim1, 
	            tclasscat_offset, i__1, i__2, i__3;

            /* Local variables */
            int i__, j, k, l;
            double u;
            int nc;
            double[] dn = new double[32];
            int nn, mt;
            double su, rld, pdo, rrd, rln, pno;
            int nsp;
            double rrn;
            int nnz, lcat, ntie;
            double crit;
            int nhit, mvar;
            double crit0;
            int ncmax;
            double xrand;
            double critmax;
            int ncsplit;

            /* Parameter adjustments */
            --cat;
            --win;
            --ncase;
            --cl;
            b_dim1 = mdim;
            b_offset = 1 + b_dim1;
            b -= b_offset;
            a_dim1 = mdim;
            a_offset = 1 + a_dim1;
            a -= a_offset;
            --wl;
            --wr;
            tclasscat_dim1 = nclass;
            tclasscat_offset = 1 + tclasscat_dim1;
            tclasscat -= tclasscat_offset;
            --tclasspop;
            --mind;

            /* Function Body */
            ncmax = 10;
            ncsplit = 512;
        /*     compute initial values of numerator and denominator of Gini */
            pno = 0.0;
            pdo = 0.0;
            i__1 = nclass;
            for (j = 1; j <= i__1; ++j) {
	        pno += tclasspop[j] * tclasspop[j];
	        pdo += tclasspop[j];
            }
            crit0 = pno / pdo;
            jstat = 0;
        /*     start main loop through variables to find best split */
            critmax = -1e25f;
            i__1 = mred;
            for (k = 1; k <= i__1; ++k) {
	        mind[k] = k;
            }
            nn = mred;
        /*     sampling mtry variables w/o replacement. */
            i__1 = mtry;
            for (mt = 1; mt <= i__1; ++mt) {
                xrand = oRandom.NextDouble();
	        //rrand_(&xrand);
	        j = (int) (nn * xrand) + 1;
	        mvar = mind[j];
	        mind[j] = mind[nn];
	        mind[nn] = mvar;
	        --nn;
	        lcat = cat[mvar]; // failed here on mvar
	        if (lcat == 1) {
        /*     Split on a numerical predictor. */
	            rrn = pno;
	            rrd = pdo;
	            rln = 0.0;
	            rld = 0.0;
	            zervr_(wl+1, ref nclass);
	            i__2 = nclass;
	            for (j = 1; j <= i__2; ++j) {
		        wr[j] = tclasspop[j];
	            }
	            ntie = 1;
	            i__2 = ndend - 1;
	            for (nsp = ndstart; nsp <= i__2; ++nsp) {
		        nc = a[mvar + nsp * a_dim1];
		        u = win[nc];
		        k = cl[nc];
		        rln += u * (wl[k] * 2 + u);
		        rrn += u * (wr[k] * -2 + u);
		        rld += u;
		        rrd -= u;
		        wl[k] += u;
		        wr[k] -= u;
		        if (b[mvar + nc * b_dim1] < b[mvar + a[mvar + (nsp + 1) * 
			        a_dim1] * b_dim1]) {
        /*     If neither nodes is empty, check the split. */
		            if (Math.Min(rrd, rld) > 1e-5f) {
			        crit = rln / rld + rrn / rrd;
			        if (crit > critmax) {
			            nbest = nsp;
			            critmax = crit;
			            msplit = mvar;
			            ntie = 1;
			        }
        /*     Break ties at random: */
			        if (crit == critmax) {
                        xrand = oRandom.NextDouble();
			            //rrand_(&xrand);
			            if (xrand < 1.0 / Convert.ToDouble(ntie)) {
				        nbest = nsp;
				        critmax = crit;
				        msplit = mvar;
			            }
			            ++ntie;
			        }
		            }
		        }
	            }
	        } else {
        /*     Split on a categorical predictor.  Compute the decrease in impurity. */
	            zermr_(tclasscat + tclasscat_offset, ref nclass, ref c__32);
	            i__2 = ndend;
	            for (nsp = ndstart; nsp <= i__2; ++nsp) {
		        nc = ncase[nsp];
		        l = a[mvar + ncase[nsp] * a_dim1];
		        tclasscat[cl[nc] + l * tclasscat_dim1] += win[nc];
	            }
	            nnz = 0;
	            i__2 = lcat;
	            for (i__ = 1; i__ <= i__2; ++i__) {
		        su = 0.0;
		        i__3 = nclass;
		        for (j = 1; j <= i__3; ++j) {
		            su += tclasscat[j + i__ * tclasscat_dim1];
		        }
		        dn[i__ - 1] = su;
		        if (su > 0.0) {
		            ++nnz;
		        }
	            }
	            nhit = 0;
	            if (nnz > 1) {
		        if (nclass == 2 && lcat > ncmax) {
		            classTree.catmaxb(ref pdo, tclasscat + tclasscat_offset, ++tclasspop, 
                        ref nclass, ref lcat, ref nbest, ref critmax, ref nhit, ref dn);
		        } else {
		            classTree.catmax(ref pdo, tclasscat + tclasscat_offset, ++tclasspop,
			             ref nclass, ref lcat, ref nbest, ref critmax, ref nhit, ref maxcat, 
                         ref ncmax, ref ncsplit, oRandom);
                }
		        if (nhit == 1) {
		            msplit = mvar;
		        }
        /*            else */
        /*               critmax = -1.0e25 */
	            }
	        }
            }
            if (critmax < -1e10f || msplit == 0) {
	        jstat = -1;
            }
            decsplit = critmax - crit0;
            return 0;
        } /* findbestsplit_ */

        /*     SUBROUTINE MOVEDATA */
        /*     This subroutine is the heart of the buildtree construction. Based on the */
        /*     best split the data in the part of the a matrix corresponding to the */
        /*     current node is moved to the left if it belongs to the left child and */
        /*     right if it belongs to the right child. */
        private static int movedata(CArray<int> a, CArray<int> ta, ref int mdim, 
	        ref int nsample, ref int ndstart, ref int ndend, CArray<int> idmove, 
	        CArray<int> ncase, ref int msplit, CArray<int> cat, ref int nbest, 
	        ref int ndendl)
        {
            /* System generated locals */
            int a_dim1, a_offset, i__1, i__2;

            /* Local variables */
            int k, l, n, nc, ih, ndo, msh, nsp;
            int[] icat = new int[32];

        /*     compute idmove=indicator of case nos. going left */
            /* Parameter adjustments */
            --cat;
            --ncase;
            --idmove;
            --ta;
            a_dim1 = mdim;
            a_offset = 1 + a_dim1;
            a -= a_offset;

            /* Function Body */
            if (cat[msplit] == 1) {
	        i__1 = nbest;
	        for (nsp = ndstart; nsp <= i__1; ++nsp) {
	            nc = a[msplit + nsp * a_dim1];
	            idmove[nc] = 1;
	        }
	        i__1 = ndend;
	        for (nsp = nbest + 1; nsp <= i__1; ++nsp) {
	            nc = a[msplit + nsp * a_dim1];
	            idmove[nc] = 0;
	        }
	        ndendl = nbest;
            } else {
	        ndendl = ndstart - 1;
	        l = cat[msplit];
	        rfutils.unpack(l, Convert.ToUInt32(nbest), ref icat);
	        i__1 = ndend;
	        for (nsp = ndstart; nsp <= i__1; ++nsp) {
	            nc = ncase[nsp];
	            if (icat[a[msplit + nc * a_dim1] - 1] == 1) {
		        idmove[nc] = 1;
		        ++(ndendl);
	            } else {
		        idmove[nc] = 0;
	            }
	        }
            }
        /*     shift case. nos. right and left for numerical variables. */
            i__1 = mdim;
            for (msh = 1; msh <= i__1; ++msh) {
	        if (cat[msh] == 1) {
	            k = ndstart - 1;
	            i__2 = ndend;
	            for (n = ndstart; n <= i__2; ++n) {
		        ih = a[msh + n * a_dim1];
		        if (idmove[ih] == 1) {
		            ++k;
		            ta[k] = a[msh + n * a_dim1];
		        }
        /* L50: */
	            }
	            i__2 = ndend;
	            for (n = ndstart; n <= i__2; ++n) {
		        ih = a[msh + n * a_dim1];
		        if (idmove[ih] == 0) {
		            ++k;
		            ta[k] = a[msh + n * a_dim1];
		        }
        /* L60: */
	            }
	            i__2 = ndend;
	            for (k = ndstart; k <= i__2; ++k) {
		        a[msh + k * a_dim1] = ta[k];
        /* L70: */
	            }
	        }
        /* L40: */
            }
            ndo = 0;
            if (ndo == 1) {
	        i__1 = mdim;
	        for (msh = 1; msh <= i__1; ++msh) {
	            if (cat[msh] > 1) {
		        k = ndstart - 1;
		        i__2 = ndend;
		        for (n = ndstart; n <= i__2; ++n) {
		            ih = ncase[n];
		            if (idmove[ih] == 1) {
			        ++k;
			        ta[k] = a[msh + ih * a_dim1];
		            }
        /* L150: */
		        }
		        i__2 = ndend;
		        for (n = ndstart; n <= i__2; ++n) {
		            ih = ncase[n];
		            if (idmove[ih] == 0) {
			        ++k;
			        ta[k] = a[msh + ih * a_dim1];
		            }
        /* L160: */
		        }
		        i__2 = ndend;
		        for (k = ndstart; k <= i__2; ++k) {
		            a[msh + k * a_dim1] = ta[k];
        /* L170: */
		        }
	            }
        /* L140: */
	        }
            }
        /*     compute case nos. for right and left nodes. */
            if (cat[msplit] == 1) {
	        i__1 = ndend;
	        for (n = ndstart; n <= i__1; ++n) {
	            ncase[n] = a[msplit + n * a_dim1];
        /* L80: */
	        }
            } else {
	        k = ndstart - 1;
	        i__1 = ndend;
	        for (n = ndstart; n <= i__1; ++n) {
	            if (idmove[ncase[n]] == 1) {
		        ++k;
		        ta[k] = ncase[n];
	            }
        /* L90: */
	        }
	        i__1 = ndend;
	        for (n = ndstart; n <= i__1; ++n) {
	            if (idmove[ncase[n]] == 0) {
		        ++k;
		        ta[k] = ncase[n];
	            }
        /* L100: */
	        }
	        i__1 = ndend;
	        for (k = ndstart; k <= i__1; ++k) {
	            ncase[k] = ta[k];
        /* L110: */
	        }
            }
            return 0;
        } /* movedata_ */

        /*      subroutine myunpack(l,npack,icat) */

        /*     npack is a long integer.  The sub. returns icat, an integer of zeroes and */
        /*     ones corresponding to the coefficients in the binary expansion of npack. */

        /*      integer icat(32),npack */
        /*      do j=1,32 */
        /*         icat(j)=0 */
        /*      end do */
        /*      n=npack */
        /*      icat(1)=mod(n,2) */
        /*      do k=2,l */
        /*         n=(n-icat(k-1))/2 */
        /*         icat(k)=mod(n,2) */
        /*      end do */
        /*      end */
    }
}
