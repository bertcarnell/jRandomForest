package randomForestTest;
import static org.junit.Assert.*;

import java.util.Random;

import org.junit.Test;

import randomForest.ScootArray;
import randomForest.rf;

/// <summary>
///This is a test class for rfTest and is intended
///to contain all rfTest Unit Tests
///</summary>
public class rfTest
{


//    private TestContext testContextInstance;

    /// <summary>
    ///Gets or sets the test context which provides
    ///information about and functionality for the current test run.
    ///</summary>
//    public TestContext TestContext
//    {
//        get
//        {
//            return testContextInstance;
//        }
//        set
//        {
//            testContextInstance = value;
//        }
//    }
//
//    #region Additional test attributes
    // 
    //You can use the following additional attributes as you write your tests:
    //
    //Use ClassInitialize to run code before running the first test in the class
    //[ClassInitialize()]
    //public static void MyClassInitialize(TestContext testContext)
    //{
    //}
    //
    //Use ClassCleanup to run code after all tests in a class have run
    //[ClassCleanup()]
    //public static void MyClassCleanup()
    //{
    //}
    //
    //Use TestInitialize to run code before running each test
    //[TestInitialize()]
    //public void MyTestInitialize()
    //{
    //}
    //
    //Use TestCleanup to run code after each test has run
    //[TestCleanup()]
    //public void MyTestCleanup()
    //{
    //}
    //
//    #endregion


    /// <summary>
    ///A test for classRF
    ///</summary>
//    [TestMethod()]
//    [DeploymentItem("CsRandomForest.dll")]
    @Test
    public void classRFTest()
    {
        /*     [,1] [,2] [,3] [,4]
       [1,]    1    2    3    4
       [2,]    5    6    7    8
       [3,]    9   10   11   12*/
        double[] x = { 1.0, 5.0, 9.0, 2.0, 6.0, 10.0, 3.0, 7.0, 11.0, 4.0, 8.0, 12.0 };
        int[] dimx = { 3, 4 };
        Integer[] cl_mem = { 1, 1, 2, 2 };
        ScootArray<Integer> cl = new ScootArray<Integer>(cl_mem);
        int ncl = 2;
        Integer[] cat_mem = { 1, 1, 1 }; // all are continuous
        ScootArray<Integer> cat = new ScootArray<Integer>(cat_mem);
        int maxcat = 1;
        int[] sampsize = { 4 };
        int[] strata = { 0 };
        int[] Options = { 0, 0, 0, 0, 0, 0, 1, 1, 0, 0 };
        int ntree = 15;
        int nvar = 1;
        int ipi = 0;
        double[] classwt = { 1, 1 };
        double[] cut = { 0.5, 0.5 };
        int nodesize = 1;
        int[] outcl = { 0, 0, 0, 0 };
        int[] counttr = { 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] prox = { 0.0 };
        double[] imprt = { 0.0, 0.0, 0.0 };
        double[] impsd = { 0.0 };
        double[] impmat = { 0.0 };
        int nrnodes = 9;
        Integer[] ndbigtree_mem = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        ScootArray<Integer> ndbigtree = new ScootArray<Integer>(ndbigtree_mem);
        Integer[] nodestatus_mem = new Integer[15 * 9];
        ScootArray<Integer> nodestatus = new ScootArray<Integer>(nodestatus_mem);
        Integer[] bestvar_mem = new Integer[15 * 9];
        ScootArray<Integer> bestvar = new ScootArray<Integer>(bestvar_mem);
        Integer[] treemap_mem = new Integer[15 * 2 * 9];
        ScootArray<Integer> treemap = new ScootArray<Integer>(treemap_mem);
        Integer[] nodeclass_mem = new Integer[15 * 9];
        ScootArray<Integer> nodeclass = new ScootArray<Integer>(nodeclass_mem);
        Double[] xbestsplit_mem = new Double[15 * 9];
        ScootArray<Double> xbestsplit = new ScootArray<Double>(xbestsplit_mem);
        Double[] errtr_mem = new Double[(2 + 1) * 15];
        ScootArray<Double> errtr = new ScootArray<Double>(errtr_mem);
        int testdat = 0;
        double[] xts = { 0.0 };
        int[] clts = { 0 };
        int nts = 1;
        double[] countts = { 0.0, 0.0 };
        int[] outclts = { 0 };
        int labelts = 0;
        double[] proxts = { 0.0 };
        Double[] errts_mem = { 0.0 };
        ScootArray<Double> errts = new ScootArray<Double>(errts_mem);
        int[] inbag = { 0, 0, 0 };

        Random oRandom = new Random(1976);
        try {
            rf.classRF(x, dimx, cl, ncl, cat, maxcat, sampsize, strata, Options, 
                    ntree, nvar, ipi, classwt, cut, nodesize, outcl, counttr, prox, imprt, 
                    impsd, impmat, nrnodes, ndbigtree, nodestatus, bestvar, treemap, 
                    nodeclass, xbestsplit, errtr, testdat, xts, clts, nts, countts, outclts, 
                    labelts, proxts, errts, inbag, oRandom);
        } catch (Exception e) {
            e.printStackTrace();
        }
        

        int[] outclExpected = { 1, 1, 2, 2 };

        for (int i = 0; i < outcl.length; i++)
            assertEquals(outclExpected[i], outcl[i]);

        oRandom = new Random(1976);
        try {
            rf.classRF(x, dimx, cl, ncl, cat, maxcat, sampsize, strata, Options,
                    ntree, nvar, ipi, classwt, cut, nodesize, outcl, counttr, prox, imprt,
                    impsd, impmat, nrnodes, ndbigtree, nodestatus, bestvar, treemap,
                    nodeclass, xbestsplit, errtr, testdat, xts, clts, nts, countts, outclts,
                    labelts, proxts, errts, inbag, oRandom);
        } catch (Exception e) {
            e.printStackTrace();
        }

        for (int i = 0; i < outcl.length; i++)
            assertEquals(outclExpected[i], outcl[i]);

        oRandom = new Random(2000);
        try {
            rf.classRF(x, dimx, cl, ncl, cat, maxcat, sampsize, strata, Options,
                    ntree, nvar, ipi, classwt, cut, nodesize, outcl, counttr, prox, imprt,
                    impsd, impmat, nrnodes, ndbigtree, nodestatus, bestvar, treemap,
                    nodeclass, xbestsplit, errtr, testdat, xts, clts, nts, countts, outclts,
                    labelts, proxts, errts, inbag, oRandom);
        } catch (Exception e) {
            e.printStackTrace();
        }

        int[] outclExpected2 = { 1, 1, 1, 2 };

        for (int i = 0; i < outcl.length; i++)
            assertEquals(outclExpected2[i], outcl[i]);

        oRandom = new Random(2000);
        try {
            rf.classRF(x, dimx, cl, ncl, cat, maxcat, sampsize, strata, Options,
                    ntree, nvar, ipi, classwt, cut, nodesize, outcl, counttr, prox, imprt,
                    impsd, impmat, nrnodes, ndbigtree, nodestatus, bestvar, treemap,
                    nodeclass, xbestsplit, errtr, testdat, xts, clts, nts, countts, outclts,
                    labelts, proxts, errts, inbag, oRandom);
        } catch (Exception e) {
            e.printStackTrace();
        }

        for (int i = 0; i < outcl.length; i++)
            assertEquals(outclExpected2[i], outcl[i]);
    }

 }
