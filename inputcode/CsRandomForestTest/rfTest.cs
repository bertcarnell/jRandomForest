using CsRandomForest;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Runtime.Serialization.Formatters.Binary;
using System.Linq;

namespace CsRandomForestTest
{
    
    
    /// <summary>
    ///This is a test class for rfTest and is intended
    ///to contain all rfTest Unit Tests
    ///</summary>
    [TestClass()]
    public class rfTest
    {


        private TestContext testContextInstance;

        /// <summary>
        ///Gets or sets the test context which provides
        ///information about and functionality for the current test run.
        ///</summary>
        public TestContext TestContext
        {
            get
            {
                return testContextInstance;
            }
            set
            {
                testContextInstance = value;
            }
        }

        #region Additional test attributes
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
        #endregion


        /// <summary>
        ///A test for classRF
        ///</summary>
        [TestMethod()]
        [DeploymentItem("CsRandomForest.dll")]
        public void classRFTest()
        {
            /*     [,1] [,2] [,3] [,4]
           [1,]    1    2    3    4
           [2,]    5    6    7    8
           [3,]    9   10   11   12*/
            double[] x = { 1.0, 5.0, 9.0, 2.0, 6.0, 10.0, 3.0, 7.0, 11.0, 4.0, 8.0, 12.0 };
            int[] dimx = { 3, 4 };
            int[] cl_mem = { 1, 1, 2, 2 };
            CArray<int> cl = new CArray<int>(cl_mem);
            int ncl = 2;
            int[] cat_mem = { 1, 1, 1 }; // all are continuous
            CArray<int> cat = new CArray<int>(cat_mem);
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
            int[] ndbigtree_mem = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
            CArray<int> ndbigtree = new CArray<int>(ndbigtree_mem);
            int[] nodestatus_mem = new int[15 * 9];
            CArray<int> nodestatus = new CArray<int>(nodestatus_mem);
            int[] bestvar_mem = new int[15 * 9];
            CArray<int> bestvar = new CArray<int>(bestvar_mem);
            int[] treemap_mem = new int[15 * 2 * 9];
            CArray<int> treemap = new CArray<int>(treemap_mem);
            int[] nodeclass_mem = new int[15 * 9];
            CArray<int> nodeclass = new CArray<int>(nodeclass_mem);
            double[] xbestsplit_mem = new double[15 * 9];
            CArray<double> xbestsplit = new CArray<double>(xbestsplit_mem);
            double[] errtr_mem = new double[(2 + 1) * 15];
            CArray<double> errtr = new CArray<double>(errtr_mem);
            int testdat = 0;
            double[] xts = { 0.0 };
            int[] clts = { 0 };
            int nts = 1;
            double[] countts = { 0.0, 0.0 };
            int[] outclts = { 0 };
            int labelts = 0;
            double[] proxts = { 0.0 };
            double[] errts_mem = { 0.0 };
            CArray<double> errts = new CArray<double>(errts_mem);
            int[] inbag = { 0, 0, 0 };

            Random oRandom = new Random(1976);
            rf.classRF(x, dimx, cl, ref ncl, cat, ref maxcat, sampsize, strata, Options, 
                ref ntree, ref nvar, ref ipi, classwt, cut, ref nodesize, outcl, counttr, prox, imprt, 
                impsd, impmat, ref nrnodes, ndbigtree, nodestatus, bestvar, treemap, 
                nodeclass, xbestsplit, errtr, ref testdat, xts, clts, ref nts, countts, outclts, 
                ref labelts, proxts, errts, inbag, oRandom);

            int[] outclExpected = { 1, 1, 2, 2 };

            for (int i = 0; i < outcl.Length; i++)
                Assert.AreEqual(outclExpected[i], outcl[i]);

            oRandom = new Random(1976);
            rf.classRF(x, dimx, cl, ref ncl, cat, ref maxcat, sampsize, strata, Options,
                ref ntree, ref nvar, ref ipi, classwt, cut, ref nodesize, outcl, counttr, prox, imprt,
                impsd, impmat, ref nrnodes, ndbigtree, nodestatus, bestvar, treemap,
                nodeclass, xbestsplit, errtr, ref testdat, xts, clts, ref nts, countts, outclts,
                ref labelts, proxts, errts, inbag, oRandom);

            for (int i = 0; i < outcl.Length; i++)
                Assert.AreEqual(outclExpected[i], outcl[i]);

            oRandom = new Random(2000);
            rf.classRF(x, dimx, cl, ref ncl, cat, ref maxcat, sampsize, strata, Options,
                ref ntree, ref nvar, ref ipi, classwt, cut, ref nodesize, outcl, counttr, prox, imprt,
                impsd, impmat, ref nrnodes, ndbigtree, nodestatus, bestvar, treemap,
                nodeclass, xbestsplit, errtr, ref testdat, xts, clts, ref nts, countts, outclts,
                ref labelts, proxts, errts, inbag, oRandom);

            int[] outclExpected2 = { 1, 1, 1, 2 };

            for (int i = 0; i < outcl.Length; i++)
                Assert.AreEqual(outclExpected2[i], outcl[i]);

            oRandom = new Random(2000);
            rf.classRF(x, dimx, cl, ref ncl, cat, ref maxcat, sampsize, strata, Options,
                ref ntree, ref nvar, ref ipi, classwt, cut, ref nodesize, outcl, counttr, prox, imprt,
                impsd, impmat, ref nrnodes, ndbigtree, nodestatus, bestvar, treemap,
                nodeclass, xbestsplit, errtr, ref testdat, xts, clts, ref nts, countts, outclts,
                ref labelts, proxts, errts, inbag, oRandom);

            for (int i = 0; i < outcl.Length; i++)
                Assert.AreEqual(outclExpected2[i], outcl[i]);
        }

     }
}
