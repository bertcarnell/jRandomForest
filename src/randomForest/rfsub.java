package randomForest;

import java.util.Random;

/// <summary>
/// Random Forest subroutines
/// </summary>
public class rfsub {
    /* Table of constant values */

    private static int c__32 = 32;

    private static int zerv_(ScootArray<Integer> ix, int m1) {
        /* System generated locals */
        int i__1;

        /* Local variables */
        int n;

        /* Parameter adjustments */
        ix.scootLeft();
        ;

        /* Function Body */
        i__1 = m1;
        for (n = 1; n <= i__1; ++n) {
            ix.set(n, 0);
        }
        return 0;
    } /* zerv_ */

    /* Subroutine */
    private static int zermr_(ScootArray<Double> rx, int m1, int m2) {
        /* System generated locals */
        int rx_dim1, rx_offset, i__1, i__2;

        /* Local variables */
        int i__, j;

        /* Parameter adjustments */
        rx_dim1 = m1;
        rx_offset = 1 + rx_dim1;
        rx.scootLeft(rx_offset);

        /* Function Body */
        i__1 = m1;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = m2;
            for (j = 1; j <= i__2; ++j) {
                rx.set(i__ + j * rx_dim1, 0.0);
            }
        }
        return 0;
    } /* zermr_ */

    private static int zervr_(ScootArray<Double> rx, int m1) {
        /* System generated locals */
        int i__1;

        /* Local variables */
        int n;

        /* Parameter adjustments */
        rx.scootLeft();

        /* Function Body */
        i__1 = m1;
        for (n = 1; n <= i__1; ++n) {
            rx.set(n, 0.0);
        }
        return 0;
    } /* zervr_ */

    @SuppressWarnings("unused")
    private static int zerm_(ScootArray<Integer> mx, int m1, int m2) {
        /* System generated locals */
        int mx_dim1, mx_offset, i__1, i__2;

        /* Local variables */
        int i__, j;

        /* Parameter adjustments */
        mx_dim1 = m1;
        mx_offset = 1 + mx_dim1;
        mx.scootLeft(mx_offset);

        /* Function Body */
        i__1 = m1;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = m2;
            for (j = 1; j <= i__2; ++j) {
                mx.set(i__ + j * mx_dim1, 0);
            }
        }
        return 0;
    } /* zerm_ */

    @SuppressWarnings("unused")
    private static int zermd_(ScootArray<Double> rx, int m1, int m2) {
        /* System generated locals */
        int rx_dim1, rx_offset, i__1, i__2;

        /* Local variables */
        int i__, j;

        /* Parameter adjustments */
        rx_dim1 = m1;
        rx_offset = 1 + rx_dim1;
        rx.scootLeft(rx_offset);

        /* Function Body */
        i__1 = m1;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = m2;
            for (j = 1; j <= i__2; ++j) {
                rx.set(i__ + j * rx_dim1, 0.0);
            }
        }
        return 0;
    } /* zermd_ */

    /// <summary>
    /// Buildtree consists of repeated calls to two subroutines, Findbestsplit
    /// and Movedata. Findbestsplit does
    /// just that--it finds the best split of the current node. Movedata moves
    /// the data in the split node
    /// right and left so that the data corresponding to each child node is
    /// contiguous. The buildtree bookkeeping
    /// is different from that in Friedman's original CART program.
    /// <code>ncur</code> is the total number of nodes to date.
    /// <code>nodestatus(k)=1</code> if the <code>kth</code> node has been
    /// split. <code>nodestatus(k)=2</code> if the node exists but has not yet
    /// been
    /// split, and <code>=-1</code> of the node is terminal. A node is terminal
    /// if its size is below a threshold value, or
    /// if it is all one class, or if all the x-values are equal. If the current
    /// node <code>k</code> is split, then its
    /// children are numbered <code>ncur+1</code> (left), and
    /// <code>ncur+2</code> (right), <code>ncur</code> increases to
    /// <code>ncur+2</code> and the next node
    /// to be split is numbered <code>k+1</code>. When no more nodes can be
    /// split, buildtree returns to the main program.
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
    public static int buildtree(ScootArray<Integer> a, ScootArray<Integer> b,
            ScootArray<Integer> cl,
            ScootArray<Integer> cat, int maxcat, int mdim, int nsample,
            int nclass,
            ScootArray<Integer> treemap, ScootArray<Integer> bestvar,
            ScootArray<Integer> bestsplit, ScootArray<Integer> bestsplitnext,
            ScootArray<Double> tgini,
            ScootArray<Integer> nodestatus, ScootArray<Integer> nodepop,
            ScootArray<Integer> nodestart, ScootArray<Double> classpop,
            ScootArray<Double> tclasspop,
            ScootArray<Double> tclasscat,
            ScootArray<Integer> ta, int nrnodes, ScootArray<Integer> idmove,
            int ndsize, ScootArray<Integer> ncase,
            int mtry, ScootArray<Integer> iv, ScootArray<Integer> nodeclass,
            ScootArray<Integer> ndbigtree, ScootArray<Double> win,
            ScootArray<Double> wr, ScootArray<Double> wl,
            int mred,
            int nuse, ScootArray<Integer> mind,
            Random oRandom) {
        /* System generated locals */
        int a_dim1, a_offset, b_dim1, b_offset, classpop_dim1,
                classpop_offset, tclasscat_dim1, tclasscat_offset, i__1, i__2;

        /* Local variables */
        double decsplit = 0.0;
        int j, n, nc;
//        int k, kn, ntie;
//        double pp;
        int ncur;
        double popt1, popt2;
        int ndend;
        int nbest = 0;
//        double xrand;
        int jstat;
        int ndendl = 0;
        int kbuild, msplit, ndstart;

        /* Parameter adjustments */
        tgini.scootLeft();
        cat.scootLeft();
        win.scootLeft();
        ncase.scootLeft();
        idmove.scootLeft();
        ta.scootLeft();
        cl.scootLeft();
        b_dim1 = mdim;
        b_offset = 1 + b_dim1;
        b.scootLeft(b_offset);
        a_dim1 = mdim;
        a_offset = 1 + a_dim1;
        a.scootLeft(a_offset);
        wl.scootLeft();
        wr.scootLeft();
        tclasscat_dim1 = nclass;
        tclasscat_offset = 1 + tclasscat_dim1;
        tclasscat.scootLeft(tclasscat_offset);
        tclasspop.scootLeft();
        nodeclass.scootLeft();
        classpop_dim1 = nclass;
        classpop_offset = 1 + classpop_dim1;
        classpop.scootLeft(classpop_offset);
        nodestart.scootLeft();
        nodepop.scootLeft();
        nodestatus.scootLeft();
        bestsplitnext.scootLeft();
        bestsplit.scootLeft();
        bestvar.scootLeft();
        treemap.scootLeft(3);
        mind.scootLeft();
        iv.scootLeft();

        /* Function Body */
        msplit = 0;
        // zerv_(&nodestatus[1], nrnodes);
        zerv_(nodestatus.scootRight(), nrnodes);
        zerv_(nodestart.scootRight(), nrnodes);
        zerv_(nodepop.scootRight(), nrnodes);
        zermr_(classpop.scootRight(classpop_offset), nclass, nrnodes);

        i__1 = nclass;
        for (j = 1; j <= i__1; ++j) {
            classpop.set(j + classpop_dim1, tclasspop.get(j));
        }
        ncur = 1;
        nodestart.set(1, 1);
        nodepop.set(1, nuse);
        nodestatus.set(1, 2);
        /* start main loop */
        i__1 = nrnodes;
        for (kbuild = 1; kbuild <= i__1; ++kbuild) {
            /* call intpr("kbuild", 6, kbuild, 1) */
            /* call intpr("ncur", 4, ncur, 1) */
            if (kbuild > ncur) {
                ndbigtree.set(0, nrnodes);
                L50(nclass, nodestatus, classpop, nrnodes, nodeclass, ndbigtree,
                        oRandom, classpop_dim1, kbuild);
                return 0;
            }
            if (nodestatus.get(kbuild) != 2) {
                continue;
            }
            /* initialize for next call to findbestsplit */
            ndstart = nodestart.get(kbuild);
            ndend = ndstart + nodepop.get(kbuild) - 1;
            i__2 = nclass;
            for (j = 1; j <= i__2; ++j) {
                tclasspop.set(j, classpop.get(j + kbuild * classpop_dim1));
            }
            jstat = 0;
            findbestsplit(a.scootRight(a_offset), b.scootRight(b_offset),
                    cl.scootRight(), mdim, nsample,
                    nclass, cat.scootRight(), maxcat, ndstart, ndend,
                    tclasspop.scootRight(),
                    tclasscat.scootRight(tclasscat_offset), msplit, decsplit,
                    nbest,
                    ncase.scootRight(), jstat, mtry, win.scootRight(),
                    wr.scootRight(), wl.scootRight(), mred, mind.scootRight(),
                    oRandom);
            /* call intpr("jstat", 5, jstat, 1) */
            /* call intpr("msplit", 6, msplit, 1) */
            /* If the node is terminal, move on. Otherwise, split. */
            if (jstat == -1) {
                nodestatus.set(kbuild, -1);
                continue;
            } else {
                bestvar.set(kbuild, msplit);
                iv.set(msplit, 1);
                if (decsplit < 0.0) {
                    decsplit = 0.0;
                }
                tgini.set(msplit, tgini.get(msplit) + decsplit);
                if (cat.get(msplit) == 1) {
                    bestsplit.set(kbuild, a.get(msplit + nbest * a_dim1));
                    bestsplitnext.set(kbuild,
                            a.get(msplit + (nbest + 1) * a_dim1));
                } else {
                    bestsplit.set(kbuild, nbest);
                    bestsplitnext.set(kbuild, 0);
                }
            }
            movedata(a.scootRight(a_offset), ta.scootRight(), mdim, nsample,
                    ndstart, ndend,
                    idmove.scootRight(), ncase.scootRight(), msplit,
                    cat.scootRight(), nbest, ndendl);

            /* call intpr("ndend", 5, ndend, 1) */
            /* call intpr("ndendl", 6, ndendl, 1) */
            /* leftnode no.= ncur+1, rightnode no. = ncur+2. */
            nodepop.set(ncur + 1, ndendl - ndstart + 1);
            nodepop.set(ncur + 2, ndend - ndendl);
            nodestart.set(ncur + 1, ndstart);
            nodestart.set(ncur + 2, ndendl + 1);
            /* find class populations in both nodes */
            i__2 = ndendl;
            for (n = ndstart; n <= i__2; ++n) {
                nc = ncase.get(n);
                j = cl.get(nc);
                classpop.set(j + (ncur + 1) * classpop_dim1,
                        classpop.get(j + (ncur + 1) * classpop_dim1)
                                + win.get(nc));
            }
            i__2 = ndend;
            for (n = ndendl + 1; n <= i__2; ++n) {
                nc = ncase.get(n);
                j = cl.get(nc);
                classpop.set(j + (ncur + 2) * classpop_dim1,
                        classpop.get(j + (ncur + 2) * classpop_dim1)
                                + win.get(nc));
            }
            /* call intpr("nL", 2, nodepop(ncur+1), 1) */
            /* call intpr("nR", 2, nodepop(ncur+2), 1) */
            /* check on nodestatus */
            nodestatus.set(ncur + 1, 2);
            nodestatus.set(ncur + 2, 2);
            if (nodepop.get(ncur + 1) <= ndsize) {
                nodestatus.set(ncur + 1, -1);
            }
            if (nodepop.get(ncur + 2) <= ndsize) {
                nodestatus.set(ncur + 2, -1);
            }
            popt1 = 0.0;
            popt2 = 0.0;
            i__2 = nclass;
            for (j = 1; j <= i__2; ++j) {
                popt1 += classpop.get(j + (ncur + 1) * classpop_dim1);
                popt2 += classpop.get(j + (ncur + 2) * classpop_dim1);
            }
            i__2 = nclass;
            for (j = 1; j <= i__2; ++j) {
                if (classpop.get(j + (ncur + 1) * classpop_dim1) == popt1) {
                    nodestatus.set(ncur + 1, -1);
                }
                if (classpop.get(j + (ncur + 2) * classpop_dim1) == popt2) {
                    nodestatus.set(ncur + 2, -1);
                }
            }
            treemap.set((kbuild << 1) + 1, ncur + 1);
            treemap.set((kbuild << 1) + 2, ncur + 2);
            nodestatus.set(kbuild, 1);
            ncur += 2;
            if (ncur >= nrnodes) {
                ndbigtree.set(0, nrnodes);
                L50(nclass, nodestatus, classpop, nrnodes, nodeclass, ndbigtree,
                        oRandom, classpop_dim1, kbuild);
                return 0;
            }
        }
        ndbigtree.set(0, nrnodes);
        L50(nclass, nodestatus, classpop, nrnodes, nodeclass, ndbigtree,
                oRandom, classpop_dim1, kbuild);
        return 0;
    } /* buildtree_ */

    private static void L50(int nclass, ScootArray<Integer> nodestatus,
            ScootArray<Double> classpop, int nrnodes,
            ScootArray<Integer> nodeclass, ScootArray<Integer> ndbigtree,
            Random oRandom, int classpop_dim1, int kbuild) {
        int i__1;
        int i__2;
        int j;
        int k;
        int kn;
        double pp;
        int ntie;
        double xrand;
        for (k = nrnodes; k >= 1; --k) {
            if (nodestatus.get(k) == 0) {
                ndbigtree.set(0, ndbigtree.get(0) - 1);
            }
            if (nodestatus.get(k) == 2) {
                nodestatus.set(kbuild, -1);
            }
        }
        /* form prediction in terminal nodes */
        i__1 = ndbigtree.get(0);
        for (kn = 1; kn <= i__1; ++kn) {
            if (nodestatus.get(kn) == -1) {
                pp = 0.0;
                ntie = 1;
                i__2 = nclass;
                for (j = 1; j <= i__2; ++j) {
                    if (classpop.get(j + kn * classpop_dim1) > pp) {
                        nodeclass.set(kn, j);
                        pp = classpop.get(j + kn * classpop_dim1);
                        ntie = 1;
                    }
                    /* Break ties at random: */
                    if (classpop.get(j + kn * classpop_dim1) == pp) {
                        xrand = oRandom.nextDouble();
                        // rrand_(&xrand);
                        if (xrand < 1.0 / (double) (ntie)) {
                            nodeclass.set(kn, j);
                            pp = classpop.get(j + kn * classpop_dim1);
                        }
                        ++ntie;
                    }
                }
            }
            /* call intpr("node", 4, kn, 1) */
            /* call intpr("status", 6, nodestatus(kn), 1) */
            /* call intpr("pred", 4, nodeclass(kn), 1) */
            /* call dblepr("pop1", 4, classpop(1, kn), 1) */
            /* call dblepr("pop2", 4, classpop(2, kn), 1) */
        }
    }

    /* SUBROUTINE FINDBESTSPLIT */
    /* For the best split, msplit is the variable split on. decsplit is the */
    /* dec. in impurity. If msplit is numerical, nsplit is the case number */
    /* of value of msplit split on, and nsplitnext is the case number of the */
    /* next larger value of msplit. If msplit is categorical, then nsplit is */
    /* the coding into an integer of the categories going left. */
    private static int findbestsplit(ScootArray<Integer> a,
            ScootArray<Integer> b,
            ScootArray<Integer> cl,
            int mdim, int nsample, int nclass, ScootArray<Integer> cat,
            int maxcat, int ndstart, int ndend, ScootArray<Double> tclasspop,
            ScootArray<Double> tclasscat, int msplit, double decsplit,
            int nbest, ScootArray<Integer> ncase, int jstat, int mtry,
            ScootArray<Double> win, ScootArray<Double> wr,
            ScootArray<Double> wl,
            int mred,
            ScootArray<Integer> mind, Random oRandom) {
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
        cat.scootLeft();
        win.scootLeft();
        ncase.scootLeft();
        cl.scootLeft();
        b_dim1 = mdim;
        b_offset = 1 + b_dim1;
        b.scootLeft(b_offset);
        a_dim1 = mdim;
        a_offset = 1 + a_dim1;
        a.scootLeft(a_offset);
        wl.scootLeft();
        wr.scootLeft();
        tclasscat_dim1 = nclass;
        tclasscat_offset = 1 + tclasscat_dim1;
        tclasscat.scootLeft(tclasscat_offset);
        tclasspop.scootLeft();
        mind.scootLeft();

        /* Function Body */
        ncmax = 10;
        ncsplit = 512;
        /* compute initial values of numerator and denominator of Gini */
        pno = 0.0;
        pdo = 0.0;
        i__1 = nclass;
        for (j = 1; j <= i__1; ++j) {
            pno += tclasspop.get(j) * tclasspop.get(j);
            pdo += tclasspop.get(j);
        }
        crit0 = pno / pdo;
        jstat = 0;
        /* start main loop through variables to find best split */
        critmax = -1e25f;
        i__1 = mred;
        for (k = 1; k <= i__1; ++k) {
            mind.set(k, k);
        }
        nn = mred;
        /* sampling mtry variables w/o replacement. */
        i__1 = mtry;
        for (mt = 1; mt <= i__1; ++mt) {
            xrand = oRandom.nextDouble();
            // rrand_(&xrand);
            j = (int) (nn * xrand) + 1;
            mvar = mind.get(j);
            mind.set(j, mind.get(nn));
            mind.set(nn, mvar);
            --nn;
            lcat = cat.get(mvar); // failed here on mvar
            if (lcat == 1) {
                /* Split on a numerical predictor. */
                rrn = pno;
                rrd = pdo;
                rln = 0.0;
                rld = 0.0;
                zervr_(wl.scootRight(), nclass);
                i__2 = nclass;
                for (j = 1; j <= i__2; ++j) {
                    wr.set(j, tclasspop.get(j));
                }
                ntie = 1;
                i__2 = ndend - 1;
                for (nsp = ndstart; nsp <= i__2; ++nsp) {
                    nc = a.get(mvar + nsp * a_dim1);
                    u = win.get(nc);
                    k = cl.get(nc);
                    rln += u * (wl.get(k) * 2 + u);
                    rrn += u * (wr.get(k) * -2 + u);
                    rld += u;
                    rrd -= u;
                    wl.set(k, wl.get(k) + u);
                    wr.set(k, wr.get(k) - u);
                    if (b.get(mvar + nc * b_dim1) < b
                            .get(mvar + a.get(mvar + (nsp + 1) *
                                    a_dim1) * b_dim1)) {
                        /* If neither nodes is empty, check the split. */
                        if (Math.min(rrd, rld) > 1e-5f) {
                            crit = rln / rld + rrn / rrd;
                            if (crit > critmax) {
                                nbest = nsp;
                                critmax = crit;
                                msplit = mvar;
                                ntie = 1;
                            }
                            /* Break ties at random: */
                            if (crit == critmax) {
                                xrand = oRandom.nextDouble();
                                // rrand_(&xrand);
                                if (xrand < 1.0 / (double) (ntie)) {
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
                /*
                 * Split on a categorical predictor. Compute the decrease in
                 * impurity.
                 */
                zermr_(tclasscat.scootRight(tclasscat_offset), nclass, c__32);
                i__2 = ndend;
                for (nsp = ndstart; nsp <= i__2; ++nsp) {
                    nc = ncase.get(nsp);
                    l = a.get(mvar + ncase.get(nsp) * a_dim1);
                    tclasscat.set(cl.get(nc) + l * tclasscat_dim1,
                            tclasscat.get(cl.get(nc) + l * tclasscat_dim1)
                                    + win.get(nc));
                }
                nnz = 0;
                i__2 = lcat;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    su = 0.0;
                    i__3 = nclass;
                    for (j = 1; j <= i__3; ++j) {
                        su += tclasscat.get(j + i__ * tclasscat_dim1);
                    }
                    dn[i__ - 1] = su;
                    if (su > 0.0) {
                        ++nnz;
                    }
                }
                nhit = 0;
                if (nnz > 1) {
                    if (nclass == 2 && lcat > ncmax) {
                        Tree.catmaxb(pdo,
                                tclasscat.scootRight(tclasscat_offset),
                                tclasspop.scootRight(),
                                nclass, lcat, nbest, critmax, nhit, dn);
                    } else {
                        Tree.catmax(pdo, tclasscat.scootRight(tclasscat_offset),
                                tclasspop.scootRight(),
                                nclass, lcat, nbest, critmax, nhit, maxcat,
                                ncmax, ncsplit, oRandom);
                    }
                    if (nhit == 1) {
                        msplit = mvar;
                    }
                    /* else */
                    /* critmax = -1.0e25 */
                }
            }
        }
        if (critmax < -1e10f || msplit == 0) {
            jstat = -1;
        }
        decsplit = critmax - crit0;
        return 0;
    } /* findbestsplit_ */

    /* SUBROUTINE MOVEDATA */
    /*
     * This subroutine is the heart of the buildtree construction. Based on the
     */
    /* best split the data in the part of the a matrix corresponding to the */
    /* current node is moved to the left if it belongs to the left child and */
    /* right if it belongs to the right child. */
    private static int movedata(ScootArray<Integer> a, ScootArray<Integer> ta,
            int mdim,
            int nsample, int ndstart, int ndend, ScootArray<Integer> idmove,
            ScootArray<Integer> ncase, int msplit, ScootArray<Integer> cat,
            int nbest,
            int ndendl) {
        /* System generated locals */
        int a_dim1, a_offset, i__1, i__2;

        /* Local variables */
        int k, l, n, nc, ih, ndo, msh, nsp;
        int[] icat = new int[32];

        /* compute idmove=indicator of case nos. going left */
        /* Parameter adjustments */
        cat.scootLeft();
        ncase.scootLeft();
        idmove.scootLeft();
        ta.scootLeft();
        a_dim1 = mdim;
        a_offset = 1 + a_dim1;
        a.scootLeft(a_offset);

        /* Function Body */
        if (cat.get(msplit) == 1) {
            i__1 = nbest;
            for (nsp = ndstart; nsp <= i__1; ++nsp) {
                nc = a.get(msplit + nsp * a_dim1);
                idmove.set(nc, 1);
            }
            i__1 = ndend;
            for (nsp = nbest + 1; nsp <= i__1; ++nsp) {
                nc = a.get(msplit + nsp * a_dim1);
                idmove.set(nc, 0);
            }
            ndendl = nbest;
        } else {
            ndendl = ndstart - 1;
            l = cat.get(msplit);
            rfutils.unpack(l, (int) (nbest), icat);
            i__1 = ndend;
            for (nsp = ndstart; nsp <= i__1; ++nsp) {
                nc = ncase.get(nsp);
                if (icat[a.get(msplit + nc * a_dim1) - 1] == 1) {
                    idmove.set(nc, 1);
                    ++(ndendl);
                } else {
                    idmove.set(nc, 0);
                }
            }
        }
        /* shift case. nos. right and left for numerical variables. */
        i__1 = mdim;
        for (msh = 1; msh <= i__1; ++msh) {
            if (cat.get(msh) == 1) {
                k = ndstart - 1;
                i__2 = ndend;
                for (n = ndstart; n <= i__2; ++n) {
                    ih = a.get(msh + n * a_dim1);
                    if (idmove.get(ih) == 1) {
                        ++k;
                        ta.set(k, a.get(msh + n * a_dim1));
                    }
                    /* L50: */
                }
                i__2 = ndend;
                for (n = ndstart; n <= i__2; ++n) {
                    ih = a.get(msh + n * a_dim1);
                    if (idmove.get(ih) == 0) {
                        ++k;
                        ta.set(k, a.get(msh + n * a_dim1));
                    }
                    /* L60: */
                }
                i__2 = ndend;
                for (k = ndstart; k <= i__2; ++k) {
                    a.set(msh + k * a_dim1, ta.get(k));
                    /* L70: */
                }
            }
            /* L40: */
        }
        ndo = 0;
        if (ndo == 1) {
            i__1 = mdim;
            for (msh = 1; msh <= i__1; ++msh) {
                if (cat.get(msh) > 1) {
                    k = ndstart - 1;
                    i__2 = ndend;
                    for (n = ndstart; n <= i__2; ++n) {
                        ih = ncase.get(n);
                        if (idmove.get(ih) == 1) {
                            ++k;
                            ta.set(k, a.get(msh + ih * a_dim1));
                        }
                        /* L150: */
                    }
                    i__2 = ndend;
                    for (n = ndstart; n <= i__2; ++n) {
                        ih = ncase.get(n);
                        if (idmove.get(ih) == 0) {
                            ++k;
                            ta.set(k, a.get(msh + ih * a_dim1));
                        }
                        /* L160: */
                    }
                    i__2 = ndend;
                    for (k = ndstart; k <= i__2; ++k) {
                        a.set(msh + k * a_dim1, ta.get(k));
                        /* L170: */
                    }
                }
                /* L140: */
            }
        }
        /* compute case nos. for right and left nodes. */
        if (cat.get(msplit) == 1) {
            i__1 = ndend;
            for (n = ndstart; n <= i__1; ++n) {
                ncase.set(n, a.get(msplit + n * a_dim1));
                /* L80: */
            }
        } else {
            k = ndstart - 1;
            i__1 = ndend;
            for (n = ndstart; n <= i__1; ++n) {
                if (idmove.get(ncase.get(n)) == 1) {
                    ++k;
                    ta.set(k, ncase.get(n));
                }
                /* L90: */
            }
            i__1 = ndend;
            for (n = ndstart; n <= i__1; ++n) {
                if (idmove.get(ncase.get(n)) == 0) {
                    ++k;
                    ta.set(k, ncase.get(n));
                }
                /* L100: */
            }
            i__1 = ndend;
            for (k = ndstart; k <= i__1; ++k) {
                ncase.set(k, ta.get(k));
                /* L110: */
            }
        }
        return 0;
    } /* movedata_ */

    /* subroutine myunpack(l,npack,icat) */

    /*
     * npack is a long integer. The sub. returns icat, an integer of zeroes and
     */
    /*
     * ones corresponding to the coefficients in the binary expansion of npack.
     */

    /* integer icat(32),npack */
    /* do j=1,32 */
    /* icat(j)=0 */
    /* end do */
    /* n=npack */
    /* icat(1)=mod(n,2) */
    /* do k=2,l */
    /* n=(n-icat(k-1))/2 */
    /* icat(k)=mod(n,2) */
    /* end do */
    /* end */
}
