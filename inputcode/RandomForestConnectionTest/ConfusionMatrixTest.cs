/*
 * Copyright 2013 Rob Carnell
 * 
 */

using RandomForestConnection;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections;

namespace RandomForestConnectionTest
{
    
    
    /// <summary>
    ///This is a test class for ConfusionMatrixTest and is intended
    ///to contain all ConfusionMatrixTest Unit Tests
    ///</summary>
    [TestClass()]
    public class ConfusionMatrixTest
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
        ///A test for ConfusionMatrix Constructor
        ///</summary>
        [TestMethod()]
        public void ConfusionMatrixConstructorTest()
        {
            ConfusionMatrix target = new ConfusionMatrix();
            Assert.IsNotNull(target);
        }

        /// <summary>
        ///A test for CreateConfusion
        ///</summary>
        [TestMethod()]
        public void CreateConfusionTest()
        {
            int[] truthClass = {1,1,2,2,3,3,4,4,5,5};
            int[] predClass =  {1,2,2,4,3,1,4,4,1,2};
            int[] expectedClasses = {1,2,3,4,5};
            int[,] expectedConfusion = {{1,1,0,0,0},
                                        {0,1,0,1,0},
                                        {1,0,1,0,0},
                                        {0,0,0,2,0},
                                        {1,1,0,0,0}};
            double[] expectedErrorRate = {0.5, 0.5, 0.5, 0.0, 1.0};

            ConfusionMatrix actual = ConfusionMatrix.CreateConfusion(truthClass, predClass);
            // check classes
            Assert.AreEqual(expectedClasses.Length, actual.classes.Length);
            for (int i = 0; i < expectedClasses.Length; i++)
                Assert.AreEqual(expectedClasses[i], actual.classes[i]);
            // check confusion matrix
            IEnumerator ece = expectedConfusion.GetEnumerator();
            IEnumerator ae = actual.confusionMatrix.GetEnumerator();
            while (ece.MoveNext() && ae.MoveNext())
                Assert.AreEqual(ece.Current, ae.Current);
            
            // check error rate
            Assert.AreEqual(expectedErrorRate.Length, actual.error_rate.Length);
            for (int i = 0; i < expectedErrorRate.Length; i++)
                Assert.AreEqual(expectedErrorRate[i], actual.error_rate[i]);

            // check unequal lengths
            int[] truthClassBad = { 1, 1, 2, 2, 3, 3, 4, 4, 5, 5 };
            int[] predClassBad = { 1, 2};
            ExceptionAssert.Throws<ArgumentException>(() => ConfusionMatrix.CreateConfusion(truthClassBad, predClassBad));

            // check predClass not in truthClass
            int[] truthClassTemp = { 1, 2 };
            int[] predClassTemp = { 1, 3 };
            ConfusionMatrix actual2 = ConfusionMatrix.CreateConfusion(truthClassTemp, predClassTemp);
            Assert.AreEqual(actual2.classes.Length, 3);
            Assert.AreEqual(actual2.error_rate[0], 0.0);
            Assert.AreEqual(actual2.error_rate[1], 1.0);
            Assert.AreEqual(actual2.error_rate[2], 0.0);
        }
    }
}
