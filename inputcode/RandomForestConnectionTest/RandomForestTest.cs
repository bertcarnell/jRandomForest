/*
 * Copyright 2013 Rob Carnell
 * 
 */

using RandomForestConnection;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Runtime.Serialization.Formatters.Binary;
using System.IO;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace RandomForestConnectionTest
{
    /// <summary>
    ///This is a test class for RandomForestTest and is intended
    ///to contain all RandomForestTest Unit Tests
    ///</summary>
    [TestClass()]
    public class RandomForestTest
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
        /// Calculates the lenght in bytes of an object 
        /// and returns the size 
        /// </summary>
        /// <param name="TestObject"></param>
        /// <returns></returns>
        private long GetObjectSize(object TestObject)
        {
            BinaryFormatter bf = new BinaryFormatter();
            MemoryStream ms = new MemoryStream();
            byte[] Array;
            bf.Serialize(ms, TestObject);
            Array = ms.ToArray();
            return Array.Length;
        }

        /// <summary>
        ///A test for Classification
        ///</summary>
        [TestMethod()]
        [DeploymentItem("RandomForestConnection.dll")]
        public void ClassificationTest()
        {
            double[] temp = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0 };
            double[,] x = new MatrixVector<double>(temp, 4, 3).getMatrix();
            int[] y = { 1, 1, 2, 2 };
            int ntree = 15;
            int mtry = 1;
            bool bReturnImportance = true;
            bool bReturnProximity = true;
            int[] outclExpected = { 1, 1, 1, 2 };

            RandomForest target = new RandomForest();
            target.Classification(x, y, ntree, mtry, bReturnImportance, bReturnProximity);

            for (int i = 0; i < outclExpected.Length; i++)
                Assert.AreEqual(outclExpected[i], target.PredictedClass[i]);

            double[,] xtest = new double[1,3]{{1.5, 5.5, 12.1}};
            int[] expected = {1};
            RandomForest target2 = new RandomForest();
            target2.Classification(x, y, ntree, mtry, bReturnImportance, bReturnProximity, xtest);

            Assert.AreEqual(expected.Length, target2.PredictedClassTest.Length);
            Assert.AreEqual(expected[0], target2.PredictedClassTest[0]);
            Assert.AreEqual(2, target2.VotesTest.Length);
            Assert.AreEqual(1.0, target2.VotesTest[0] + target2.VotesTest[1], 1.0E-12);

            // operation tests, TODO:  add output tests
            int[] ytest = { 2 };
            int[] iCatVar = {0};
            double[] classwt = { 0.5, 0.5 };
            int[] strat = { 1, 1, 2, 2 };
            int[] samps = { 2, 1 };
            target2.Classification(x, y, ntree, mtry, bReturnImportance, bReturnProximity, xtest, ytest);
            target2.Classification(x, y, ntree, mtry, bReturnImportance, bReturnProximity, xtest, null, iCatVar);
            target2.Classification(x, y, ntree, mtry, bReturnImportance, bReturnProximity, xtest, null, null, classwt);
            target2.Classification(x, y, ntree, mtry, bReturnImportance, bReturnProximity, xtest, null, iCatVar, classwt);
            target2.Classification(x, y, ntree, mtry, bReturnImportance, bReturnProximity, xtest, null, null, null, strat, samps);
            target2.Classification(x, y, ntree, mtry, bReturnImportance, bReturnProximity, xtest, null, iCatVar, classwt, strat, samps);
        }

        /// <summary>
        /// 
        /// </summary>
        [TestMethod()]
        [DeploymentItem("RandomForestConnection.dll")]
        public void MemorySizeTestSmall()
        {
            Console.WriteLine("Small Dataset");
            double[] temp = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0 };
            double[,] x = new MatrixVector<double>(temp, 4, 3).getMatrix();
            int[] y = { 1, 1, 2, 2 };
            int ntree = 15;
            int mtry = 1;
            bool bReturnImportance = true;
            bool bReturnProximity = true;
            int[] outclExpected = { 1, 1, 1, 2 };

            long time1 = GC.GetTotalMemory(true);
            RandomForest target = new RandomForest();
            long time2 = GC.GetTotalMemory(true);
            long size1 = GetObjectSize((object)target);

            Console.WriteLine(string.Format("Estimate of initial object size: {0} MB or {1} MB", (double)size1 / 1024.0 / 1024.0, (double)(time2 - time1) / 1024.0 / 1024.0));

            long time3 = GC.GetTotalMemory(true);
            target.Classification(x, y, ntree, mtry, bReturnImportance, bReturnProximity);
            long time4 = GC.GetTotalMemory(true);
            long size2 = GetObjectSize((object)target);

            Console.WriteLine(string.Format("Estimate of final object size: {0} MB or {1} MB", (double)size2 / 1024.0 / 1024.0, (double)(time4 - time3) / 1024.0 / 1024.0));
        }

        /// <summary>
        /// 
        /// </summary>
        [TestMethod()]
        [DeploymentItem("RandomForestConnection.dll")]
        public void MemorySizeTestLarge()
        {
            bool bReturnImportance = true;
            bool bReturnProximity = true;

            Console.WriteLine("large Dataset");
            int nrow = 1000;
            int ncol = 100;
            double[] temp_large = new double[nrow * ncol];
            for (int i = 0; i < nrow * ncol; i++)
                temp_large[i] = Convert.ToDouble(i);
            double[,] x_large = new MatrixVector<double>(temp_large, nrow, ncol).getMatrix();
            int[] y_large = new int[nrow];
            Random r = new Random(12345);
            for (int i = 0; i < nrow; i++)
                y_large[i] = Convert.ToInt32(Math.Floor(r.NextDouble() * 5.0));
            int ntree = 500;
            int mtry = 15;

            long time1 = GC.GetTotalMemory(true);
            RandomForest target_large = new RandomForest();
            long time2 = GC.GetTotalMemory(true);
            long size1 = GetObjectSize((object)target_large);

            Console.WriteLine(string.Format("Estimate of initial object size: {0} MB or {1} MB", (double)size1 / 1024.0 / 1024.0, (double)(time2 - time1) / 1024.0 / 1024.0));

            long time3 = GC.GetTotalMemory(true);
            target_large.Classification(x_large, y_large, ntree, mtry, bReturnImportance, bReturnProximity);
            long time4 = GC.GetTotalMemory(true);
            long size2 = GetObjectSize((object)target_large);

            Console.WriteLine(string.Format("Estimate of final object size: {0} MB or {1} MB", (double)size2 / 1024.0 / 1024.0, (double)(time4 - time3) / 1024.0 / 1024.0));

            Assert.IsTrue(true);
        }

        /// <summary>
        ///A test for GetConfustionMatrix
        ///</summary>
        [TestMethod()]
        [DeploymentItem("RandomForestConnection.dll")]
        public void GetConfusionMatrixTest()
        {
            double[] temp = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0 };
            double[,] x = new MatrixVector<double>(temp, 4, 3).getMatrix();
            int[] y = { 1, 1, 2, 2 };
            int ntree = 15;
            int mtry = 1;
            bool bReturnImportance = true;
            bool bReturnProximity = true;
            int[] outclExpected = { 1, 1, 1, 2 };

            RandomForest target = new RandomForest();
            target.Classification(x, y, ntree, mtry, bReturnImportance, bReturnProximity);
            ConfusionMatrix result = target.ConfusionMatrix;

            int[,] confusionExpected = {{2,0},{1,1}};
            IEnumerator e1 = confusionExpected.GetEnumerator();
            IEnumerator e2 = result.confusionMatrix.GetEnumerator();

            while (e1.MoveNext() && e2.MoveNext())
            {
                Assert.AreEqual(Convert.ToInt32(e1.Current), Convert.ToInt32(e2.Current));
            }
        }

        [TestMethod()]
        [DeploymentItem("RandomForestConnection.dll")]
        public void BalancedRandomForestTest()
        {
            string csvData = @"y	x1	x2	x3	x4	x5 
            A	0.917237159	2.146691864	1.789959632	TRUE	TRUE
            A	1.63369131	1.09043989	1.113032225	TRUE	FALSE
            A	2.362172543	-0.145447382	1.340008439	TRUE	TRUE
            A	-0.101674958	1.605055905	0.288163183	TRUE	FALSE
            A	0.174883417	0.40383944	0.363902165	TRUE	FALSE
            A	0.652548908	2.361084711	-0.032802169	FALSE	TRUE
            A	3.160824583	-0.132888021	1.364784063	TRUE	TRUE
            A	1.151966315	1.012262243	0.421007316	TRUE	FALSE
            A	1.405606162	3.126022489	0.389298562	TRUE	TRUE
            A	1.336533147	0.42457774	1.076929929	TRUE	FALSE
            B	2.054645795	2.484338765	4.129466579	TRUE	TRUE
            B	1.475605642	2.091799337	0.904425324	FALSE	TRUE
            B	0.070627551	3.001334187	2.37004469	TRUE	TRUE
            B	1.989493715	2.39045891	2.438459225	TRUE	TRUE
            B	1.652313377	-0.033045832	4.17295972	FALSE	FALSE
            B	-0.320422635	2.304221126	3.148064574	TRUE	TRUE
            B	3.409706581	1.982482508	1.942468178	FALSE	TRUE
            B	0.904682213	1.999048101	2.593934708	FALSE	TRUE
            B	2.149186061	2.557527385	2.95991439	TRUE	TRUE
            B	2.397628841	0.900190642	3.076183707	TRUE	TRUE
            B	3.189796432	2.482304434	1.088716631	TRUE	TRUE
            B	0.965947704	1.89352095	2.013836751	FALSE	TRUE
            B	3.19892117	1.752271613	2.264944537	TRUE	TRUE
            B	0.580687571	1.82157053	4.017504746	TRUE	TRUE
            B	0.536576492	-0.304151433	3.454106379	FALSE	TRUE
            C	1.895614516	3.674274766	2.387585051	FALSE	FALSE
            C	5.46750355	3.646879655	3.490451827	FALSE	FALSE
            C	2.829515177	1.839708019	2.993245716	FALSE	FALSE
            C	3.500804566	3.161539634	2.05351729	TRUE	FALSE
            C	2.682532811	4.012624159	2.139786366	FALSE	FALSE";

            string[] splits = {" ", "\t", "\n", "\r"};
            List<string> csvDataList = csvData.Split(splits, StringSplitOptions.RemoveEmptyEntries).ToList();

            // remove the first 6 header items
            csvDataList.RemoveRange(0, 6);
            // replace A,B,C with 1, 2, 3
            // replace TRUE with 1 and False with 0
            for (int i = 0; i < csvDataList.Count; i++)
            {
                if (csvDataList[i] == "A") csvDataList[i] = "1";
                else if (csvDataList[i] == "B") csvDataList[i] = "2";
                else if (csvDataList[i] == "C") csvDataList[i] = "3";
                else if (csvDataList[i] == "TRUE") csvDataList[i] = "1";
                else if (csvDataList[i] == "FALSE") csvDataList[i] = "0";
            }

            int finalRows = csvDataList.Count / 6;
            double[,] x = new double[finalRows, 5];
            int[] y = new int[finalRows];

            int index = 0;
            for (int i = 0; i < finalRows; i++) // Rows
            {
                y[i] = Convert.ToInt32(csvDataList[index]);
                index++;
                for (int j = 0; j < 5; j++) // Columns
                {
                    x[i, j] = Convert.ToDouble(csvDataList[index]);
                    index++;
                }
            }

            // find the smallest class size and use that as the sample size
            List<int> distincty = y.Distinct().ToList();
            List<int> sampsize = new List<int>(distincty.Count);
            for (int i = 0; i < distincty.Count; i++)
                sampsize.Add(y.Count(a => a == distincty[i]));

            RandomForest rf = new RandomForest(100, 1976);

            rf.Classification(x, y, 500, -1, true, true, null, null, null, null, y, sampsize.ToArray());
            ConfusionMatrix cm = rf.ConfusionMatrix;

            int[,] expectedElements = {{8, 2, 0}, {2, 12, 1}, {0, 1, 4}};  // R says this should be 9, 1, 0, ...  seed is the same, checking on problem = the problem is the random number generator not the same as R
            for (int i = 0; i < expectedElements.GetLength(0); i++)
            {
                for (int j = 0; j < expectedElements.GetLength(1); j++)
                {
                    Assert.AreEqual(expectedElements[i, j], cm.confusionMatrix[i, j]);
                }
            }
            double[] expectedError = { 0.2, 0.2, 0.2 };
            for (int i = 0; i < expectedError.Length; i++)
                Assert.AreEqual(expectedError[i], cm.error_rate[i], 1E-12);

            RandomForest rf2 = new RandomForest(100, 1976);

            rf2.BalancedClassification(x, y);
            ConfusionMatrix cm2 = rf2.ConfusionMatrix;

            for (int i = 0; i < expectedElements.GetLength(0); i++)
            {
                for (int j = 0; j < expectedElements.GetLength(1); j++)
                {
                    Assert.AreEqual(expectedElements[i, j], cm2.confusionMatrix[i, j]);
                }
            }
            for (int i = 0; i < expectedError.Length; i++)
                Assert.AreEqual(expectedError[i], cm2.error_rate[i], 1E-12);

            RandomForest rf3 = new RandomForest(100, 1976);

            rf3.Classification(x, y, 500, -1, true, true, null, null, null, null, y, sampsize.ToArray());
            ConfusionMatrix cm3 = rf3.ConfusionMatrix;

            for (int i = 0; i < expectedElements.GetLength(0); i++)
            {
                for (int j = 0; j < expectedElements.GetLength(1); j++)
                {
                    Assert.AreEqual(expectedElements[i, j], cm.confusionMatrix[i, j]);
                }
            }
            for (int i = 0; i < expectedError.Length; i++)
                Assert.AreEqual(expectedError[i], cm.error_rate[i], 1E-12);

            // check to see that re-ordering the columns has no effect on the predicted classes
            RandomForest rf4 = new RandomForest(100, 1976);

            double[,] z = new double[finalRows, 5];
            for (int i = 0; i < finalRows; i++)
            {
                for (int j = 0; j < 5 ; j++)
                {
                    z[i, j] = x[i, 4-j];
                }
            }

            rf4.Classification(z, y, 500, -1, true, true, null, null, null, null, y, sampsize.ToArray());
            ConfusionMatrix cm4 = rf4.ConfusionMatrix;

            for (int i = 0; i < expectedElements.GetLength(0); i++)
            {
                for (int j = 0; j < expectedElements.GetLength(1); j++)
                {
                    Assert.AreEqual(expectedElements[i, j], cm.confusionMatrix[i, j]);
                }
            }
            for (int i = 0; i < expectedError.Length; i++)
                Assert.AreEqual(expectedError[i], cm.error_rate[i], 1E-12);
        }
 
       /// <summary>
        ///A test for classRF
        ///</summary>
        [TestMethod()]
        [DeploymentItem("CsRandomForest.dll")]
        public void irisDataTest()
        {
            IrisData irisData = new IrisData();

            RandomForest rf = new RandomForest(10, 20);
          
            rf.Classification(irisData.irisVariates, irisData.irisClasses, 500, -1, true, true);
            ConfusionMatrix cm = rf.ConfusionMatrix;

            int[,] expectedElements = {{50, 0, 0}, {0, 47, 3}, {0, 3, 47}};
            for (int i = 0; i < expectedElements.GetLength(0); i++)
            {
                for (int j = 0; j < expectedElements.GetLength(1); j++)
                {
                    Assert.AreEqual(expectedElements[i, j], cm.confusionMatrix[i, j]);
                }
            }
            double[] expectedError = { 0.0, 0.06, 0.06 };
            for (int i = 0; i < expectedError.Length; i++)
                Assert.AreEqual(expectedError[i], cm.error_rate[i], 1E-12);
            
            /*
             * test the prediction for one of the rows
             */
            RandomForest rf2 = new RandomForest(10, 20);
            rf2.Classification(irisData.irisVariates, irisData.irisClasses, 500, -1, true, true, irisData.irisVariatesTest);
            ConfusionMatrix cm2 = rf2.ConfusionMatrix;

            // Because the operations are different with the test element, the random seed yields different results
            int[,] expectedElements2 = { { 50, 0, 0 }, { 0, 47, 3 }, { 0, 4, 46 } }; 
            for (int i = 0; i < expectedElements2.GetLength(0); i++)
            {
                for (int j = 0; j < expectedElements2.GetLength(1); j++)
                {
                    Assert.AreEqual(expectedElements2[i, j], cm2.confusionMatrix[i, j]);
                }
            }
            // Because the operations are different with the test element, the random seed yields different results
            double[] expectedError2 = { 0.0, 0.06, 0.08 }; 
            for (int i = 0; i < expectedError2.Length; i++)
                Assert.AreEqual(expectedError2[i], cm2.error_rate[i], 1E-12);
            double[] expectedVotes2 = { 0.0, 0.03, 0.97 };
            Assert.AreEqual(expectedVotes2.Length, rf2.VotesTest.Length);
            for (int i = 0; i < expectedVotes2.Length; i++)
                Assert.AreEqual(expectedVotes2[i], rf2.VotesTest[i]);
            int[] expectedTestClass2 = { 3 };
            Assert.AreEqual(expectedTestClass2[0], rf2.PredictedClassTest[0]);

            /*
             * test the same seed gives the same result with the test element
             */
            rf2 = new RandomForest(10, 20);
            rf2.Classification(irisData.irisVariates, irisData.irisClasses, 500, -1, true, true, irisData.irisVariatesTest);
            cm2 = rf2.ConfusionMatrix;

            for (int i = 0; i < expectedElements2.GetLength(0); i++)
            {
                for (int j = 0; j < expectedElements2.GetLength(1); j++)
                {
                    Assert.AreEqual(expectedElements2[i, j], cm2.confusionMatrix[i, j]);
                }
            }
            for (int i = 0; i < expectedError2.Length; i++)
                Assert.AreEqual(expectedError2[i], cm2.error_rate[i], 1E-12);
            Assert.AreEqual(expectedVotes2.Length, rf2.VotesTest.Length);
            for (int i = 0; i < expectedVotes2.Length; i++)
                Assert.AreEqual(expectedVotes2[i], rf2.VotesTest[i]);
            Assert.AreEqual(expectedTestClass2[0], rf2.PredictedClassTest[0]);

            /*
             * test the same seed gives the same result with the columns of the data re-ordered
             *  - it does not
             */
            RandomForest rf21 = new RandomForest(10, 20);
            rf21.Classification(irisData.irisVariatesReordered, irisData.irisClasses, 500, -1, true, true, irisData.irisVariatesTestReordered);
            ConfusionMatrix cm21 = rf21.ConfusionMatrix;

            // Because the operations are different with the test element, the random seed yields different results
            int[,] expectedElements21 = { { 50, 0, 0 }, { 0, 47, 3 }, { 0, 5, 45 } }; 
            for (int i = 0; i < expectedElements21.GetLength(0); i++)
            {
                for (int j = 0; j < expectedElements21.GetLength(1); j++)
                {
                    Assert.AreEqual(expectedElements21[i, j], cm21.confusionMatrix[i, j]);
                }
            }
            // Because the operations are different with the test element, the random seed yields different results
            double[] expectedError21 = { 0.0, 0.06, 0.10 };
            for (int i = 0; i < expectedError21.Length; i++)
                Assert.AreEqual(expectedError21[i], cm21.error_rate[i], 1E-12);
            double[] expectedVotes21 = { 0.0, 0.03, 0.97 };
            Assert.AreEqual(expectedVotes21.Length, rf21.VotesTest.Length);
            for (int i = 0; i < expectedVotes21.Length; i++)
                Assert.AreEqual(expectedVotes21[i], rf21.VotesTest[i]);
            int[] expectedTestClass21 = { 3 };
            Assert.AreEqual(expectedTestClass21[0], rf21.PredictedClassTest[0]);

            /*
             * test this is reproducible with the same seed
             */
            rf21 = new RandomForest(10, 20);
            rf21.Classification(irisData.irisVariatesReordered, irisData.irisClasses, 500, -1, true, true, irisData.irisVariatesTestReordered);
            cm21 = rf21.ConfusionMatrix;

            for (int i = 0; i < expectedElements21.GetLength(0); i++)
            {
                for (int j = 0; j < expectedElements21.GetLength(1); j++)
                {
                    Assert.AreEqual(expectedElements21[i, j], cm21.confusionMatrix[i, j]);
                }
            }
            for (int i = 0; i < expectedError21.Length; i++)
                Assert.AreEqual(expectedError21[i], cm21.error_rate[i], 1E-12);
            Assert.AreEqual(expectedVotes21.Length, rf21.VotesTest.Length);
            for (int i = 0; i < expectedVotes21.Length; i++)
                Assert.AreEqual(expectedVotes21[i], rf21.VotesTest[i]);
            Assert.AreEqual(expectedTestClass21[0], rf21.PredictedClassTest[0]);

            /*
             * stratified sampling
             */
            RandomForest rf3 = new RandomForest(100, 200);
            rf3.Classification(irisData.irisVariates, irisData.irisClasses, 500, -1, true, true, null, null, null, null, null, new int[3] { 10, 10, 50 });
            ConfusionMatrix cm3 = rf3.ConfusionMatrix;

            int[,] expectedElements3 = { { 50, 0, 0 }, { 0, 45, 5 }, { 0, 1, 49 } };
            for (int i = 0; i < expectedElements3.GetLength(0); i++)
            {
                for (int j = 0; j < expectedElements3.GetLength(1); j++)
                {
                    Assert.AreEqual(expectedElements3[i, j], cm3.confusionMatrix[i, j]);
                }
            }
            double[] expectedError3 = { 0.0, 0.1, 0.02 };
            for (int i = 0; i < expectedError3.Length; i++)
                Assert.AreEqual(expectedError3[i], cm3.error_rate[i], 1E-12);
        }

        /// <summary>
        ///A test for classRF
        ///</summary>
        [TestMethod()]
        [DeploymentItem("CsRandomForest.dll")]
        public void irisDataTestBigSample()
        {
            IrisData irisData = new IrisData();

            RandomForest rf = new RandomForest(10, 20);

            rf.Classification(irisData.irisVariates, irisData.irisClasses, 10000, -1, true, true);
            ConfusionMatrix cm = rf.ConfusionMatrix;

            int[,] expectedElements = {{50, 0, 0}, {0, 47, 3}, {0, 3, 47}};
            for (int i = 0; i < expectedElements.GetLength(0); i++)
            {
                for (int j = 0; j < expectedElements.GetLength(1); j++)
                {
                    Assert.AreEqual(expectedElements[i, j], cm.confusionMatrix[i, j]);
                }
            }
            double[] expectedError = { 0.0, 0.06, 0.06 };
            for (int i = 0; i < expectedError.Length; i++)
                Assert.AreEqual(expectedError[i], cm.error_rate[i], 1E-12);
            
            /*
             * test the prediction for one of the rows
             */
            RandomForest rf2 = new RandomForest(10, 20);
            rf2.Classification(irisData.irisVariates, irisData.irisClasses, 10000, -1, true, true, irisData.irisVariatesTest);
            ConfusionMatrix cm2 = rf2.ConfusionMatrix;

            // Because the operations are different with the test element, the random seed yields different results for a smaller sample size
            //   Here with the bigger sample size, the results are similar to the baseline case
            for (int i = 0; i < expectedElements.GetLength(0); i++)
            {
                for (int j = 0; j < expectedElements.GetLength(1); j++)
                {
                    Assert.AreEqual(expectedElements[i, j], cm2.confusionMatrix[i, j]);
                }
            }
            // Because the operations are different with the test element, the random seed yields different results
            for (int i = 0; i < expectedError.Length; i++)
                Assert.AreEqual(expectedError[i], cm2.error_rate[i], 1E-12);
            double[] expectedVotes2 = { 0.0, 0.0237, 0.9763 };
            Assert.AreEqual(expectedVotes2.Length, rf2.VotesTest.Length);
            for (int i = 0; i < expectedVotes2.Length; i++)
                Assert.AreEqual(expectedVotes2[i], rf2.VotesTest[i]);
            int[] expectedTestClass2 = { 3 };
            Assert.AreEqual(expectedTestClass2[0], rf2.PredictedClassTest[0]);

            /*
             * test the same seed gives the same result with the test element
             */
            rf2 = new RandomForest(10, 20);
            rf2.Classification(irisData.irisVariates, irisData.irisClasses, 10000, -1, true, true, irisData.irisVariatesTest);
            cm2 = rf2.ConfusionMatrix;

            for (int i = 0; i < expectedElements.GetLength(0); i++)
            {
                for (int j = 0; j < expectedElements.GetLength(1); j++)
                {
                    Assert.AreEqual(expectedElements[i, j], cm2.confusionMatrix[i, j]);
                }
            }
            for (int i = 0; i < expectedError.Length; i++)
                Assert.AreEqual(expectedError[i], cm2.error_rate[i], 1E-12);
            Assert.AreEqual(expectedVotes2.Length, rf2.VotesTest.Length);
            for (int i = 0; i < expectedVotes2.Length; i++)
                Assert.AreEqual(expectedVotes2[i], rf2.VotesTest[i]);
            Assert.AreEqual(expectedTestClass2[0], rf2.PredictedClassTest[0]);

            /*
             * test the same seed gives the same result with the columns of the data re-ordered
             *  - it does not for the probabilities, even with 10000 samples
             */
            RandomForest rf21 = new RandomForest(10, 20);
            rf21.Classification(irisData.irisVariatesReordered, irisData.irisClasses, 10000, -1, true, true, irisData.irisVariatesTestReordered);
            ConfusionMatrix cm21 = rf21.ConfusionMatrix;

            // Because the operations are different with the test element, the random seed yields different results
            for (int i = 0; i < expectedElements.GetLength(0); i++)
            {
                for (int j = 0; j < expectedElements.GetLength(1); j++)
                {
                    Assert.AreEqual(expectedElements[i, j], cm21.confusionMatrix[i, j]);
                }
            }
            // Because the operations are different with the test element, the random seed yields different results
            for (int i = 0; i < expectedError.Length; i++)
                Assert.AreEqual(expectedError[i], cm21.error_rate[i], 1E-12);
            double[] expectedVotes21 = { 0.0, 0.0256, 0.9744 };
            Assert.AreEqual(expectedVotes21.Length, rf21.VotesTest.Length);
            for (int i = 0; i < expectedVotes21.Length; i++)
                Assert.AreEqual(expectedVotes21[i], rf21.VotesTest[i]);
            Assert.AreEqual(expectedTestClass2[0], rf21.PredictedClassTest[0]);

            /*
             * test this is reproducible with the same seed
             */
            rf21 = new RandomForest(10, 20);
            rf21.Classification(irisData.irisVariatesReordered, irisData.irisClasses, 10000, -1, true, true, irisData.irisVariatesTestReordered);
            cm21 = rf21.ConfusionMatrix;

            for (int i = 0; i < expectedElements.GetLength(0); i++)
            {
                for (int j = 0; j < expectedElements.GetLength(1); j++)
                {
                    Assert.AreEqual(expectedElements[i, j], cm21.confusionMatrix[i, j]);
                }
            }
            for (int i = 0; i < expectedError.Length; i++)
                Assert.AreEqual(expectedError[i], cm21.error_rate[i], 1E-12);
            Assert.AreEqual(expectedVotes21.Length, rf21.VotesTest.Length);
            for (int i = 0; i < expectedVotes21.Length; i++)
                Assert.AreEqual(expectedVotes21[i], rf21.VotesTest[i]);
            Assert.AreEqual(expectedTestClass2[0], rf21.PredictedClassTest[0]);
        }
    }

    class IrisData
    {
        public int[] irisClasses;
        public double[,] irisVariates;
        public double[,] irisVariatesReordered;
        public double[,] irisVariatesTest;
        public double[,] irisVariatesTestReordered;

        public IrisData()
        {
            string csvData = @"Sepal.Length	Sepal.Width	Petal.Length	Petal.Width
                    1	5.1	3.5	1.4	0.2
                    2	4.9	3.0	1.4	0.2
                    3	4.7	3.2	1.3	0.2
                    4	4.6	3.1	1.5	0.2
                    5	5.0	3.6	1.4	0.2
                    6	5.4	3.9	1.7	0.4
                    7	4.6	3.4	1.4	0.3
                    8	5.0	3.4	1.5	0.2
                    9	4.4	2.9	1.4	0.2
                    10	4.9	3.1	1.5	0.1
                    11	5.4	3.7	1.5	0.2
                    12	4.8	3.4	1.6	0.2
                    13	4.8	3.0	1.4	0.1
                    14	4.3	3.0	1.1	0.1
                    15	5.8	4.0	1.2	0.2
                    16	5.7	4.4	1.5	0.4
                    17	5.4	3.9	1.3	0.4
                    18	5.1	3.5	1.4	0.3
                    19	5.7	3.8	1.7	0.3
                    20	5.1	3.8	1.5	0.3
                    21	5.4	3.4	1.7	0.2
                    22	5.1	3.7	1.5	0.4
                    23	4.6	3.6	1.0	0.2
                    24	5.1	3.3	1.7	0.5
                    25	4.8	3.4	1.9	0.2
                    26	5.0	3.0	1.6	0.2
                    27	5.0	3.4	1.6	0.4
                    28	5.2	3.5	1.5	0.2
                    29	5.2	3.4	1.4	0.2
                    30	4.7	3.2	1.6	0.2
                    31	4.8	3.1	1.6	0.2
                    32	5.4	3.4	1.5	0.4
                    33	5.2	4.1	1.5	0.1
                    34	5.5	4.2	1.4	0.2
                    35	4.9	3.1	1.5	0.2
                    36	5.0	3.2	1.2	0.2
                    37	5.5	3.5	1.3	0.2
                    38	4.9	3.6	1.4	0.1
                    39	4.4	3.0	1.3	0.2
                    40	5.1	3.4	1.5	0.2
                    41	5.0	3.5	1.3	0.3
                    42	4.5	2.3	1.3	0.3
                    43	4.4	3.2	1.3	0.2
                    44	5.0	3.5	1.6	0.6
                    45	5.1	3.8	1.9	0.4
                    46	4.8	3.0	1.4	0.3
                    47	5.1	3.8	1.6	0.2
                    48	4.6	3.2	1.4	0.2
                    49	5.3	3.7	1.5	0.2
                    50	5.0	3.3	1.4	0.2
                    51	7.0	3.2	4.7	1.4
                    52	6.4	3.2	4.5	1.5
                    53	6.9	3.1	4.9	1.5
                    54	5.5	2.3	4.0	1.3
                    55	6.5	2.8	4.6	1.5
                    56	5.7	2.8	4.5	1.3
                    57	6.3	3.3	4.7	1.6
                    58	4.9	2.4	3.3	1.0
                    59	6.6	2.9	4.6	1.3
                    60	5.2	2.7	3.9	1.4
                    61	5.0	2.0	3.5	1.0
                    62	5.9	3.0	4.2	1.5
                    63	6.0	2.2	4.0	1.0
                    64	6.1	2.9	4.7	1.4
                    65	5.6	2.9	3.6	1.3
                    66	6.7	3.1	4.4	1.4
                    67	5.6	3.0	4.5	1.5
                    68	5.8	2.7	4.1	1.0
                    69	6.2	2.2	4.5	1.5
                    70	5.6	2.5	3.9	1.1
                    71	5.9	3.2	4.8	1.8
                    72	6.1	2.8	4.0	1.3
                    73	6.3	2.5	4.9	1.5
                    74	6.1	2.8	4.7	1.2
                    75	6.4	2.9	4.3	1.3
                    76	6.6	3.0	4.4	1.4
                    77	6.8	2.8	4.8	1.4
                    78	6.7	3.0	5.0	1.7
                    79	6.0	2.9	4.5	1.5
                    80	5.7	2.6	3.5	1.0
                    81	5.5	2.4	3.8	1.1
                    82	5.5	2.4	3.7	1.0
                    83	5.8	2.7	3.9	1.2
                    84	6.0	2.7	5.1	1.6
                    85	5.4	3.0	4.5	1.5
                    86	6.0	3.4	4.5	1.6
                    87	6.7	3.1	4.7	1.5
                    88	6.3	2.3	4.4	1.3
                    89	5.6	3.0	4.1	1.3
                    90	5.5	2.5	4.0	1.3
                    91	5.5	2.6	4.4	1.2
                    92	6.1	3.0	4.6	1.4
                    93	5.8	2.6	4.0	1.2
                    94	5.0	2.3	3.3	1.0
                    95	5.6	2.7	4.2	1.3
                    96	5.7	3.0	4.2	1.2
                    97	5.7	2.9	4.2	1.3
                    98	6.2	2.9	4.3	1.3
                    99	5.1	2.5	3.0	1.1
                    100	5.7	2.8	4.1	1.3
                    101	6.3	3.3	6.0	2.5
                    102	5.8	2.7	5.1	1.9
                    103	7.1	3.0	5.9	2.1
                    104	6.3	2.9	5.6	1.8
                    105	6.5	3.0	5.8	2.2
                    106	7.6	3.0	6.6	2.1
                    107	4.9	2.5	4.5	1.7
                    108	7.3	2.9	6.3	1.8
                    109	6.7	2.5	5.8	1.8
                    110	7.2	3.6	6.1	2.5
                    111	6.5	3.2	5.1	2.0
                    112	6.4	2.7	5.3	1.9
                    113	6.8	3.0	5.5	2.1
                    114	5.7	2.5	5.0	2.0
                    115	5.8	2.8	5.1	2.4
                    116	6.4	3.2	5.3	2.3
                    117	6.5	3.0	5.5	1.8
                    118	7.7	3.8	6.7	2.2
                    119	7.7	2.6	6.9	2.3
                    120	6.0	2.2	5.0	1.5
                    121	6.9	3.2	5.7	2.3
                    122	5.6	2.8	4.9	2.0
                    123	7.7	2.8	6.7	2.0
                    124	6.3	2.7	4.9	1.8
                    125	6.7	3.3	5.7	2.1
                    126	7.2	3.2	6.0	1.8
                    127	6.2	2.8	4.8	1.8
                    128	6.1	3.0	4.9	1.8
                    129	6.4	2.8	5.6	2.1
                    130	7.2	3.0	5.8	1.6
                    131	7.4	2.8	6.1	1.9
                    132	7.9	3.8	6.4	2.0
                    133	6.4	2.8	5.6	2.2
                    134	6.3	2.8	5.1	1.5
                    135	6.1	2.6	5.6	1.4
                    136	7.7	3.0	6.1	2.3
                    137	6.3	3.4	5.6	2.4
                    138	6.4	3.1	5.5	1.8
                    139	6.0	3.0	4.8	1.8
                    140	6.9	3.1	5.4	2.1
                    141	6.7	3.1	5.6	2.4
                    142	6.9	3.1	5.1	2.3
                    143	5.8	2.7	5.1	1.9
                    144	6.8	3.2	5.9	2.3
                    145	6.7	3.3	5.7	2.5
                    146	6.7	3.0	5.2	2.3
                    147	6.3	2.5	5.0	1.9
                    148	6.5	3.0	5.2	2.0
                    149	6.2	3.4	5.4	2.3
                    150	5.9	3.0	5.1	1.8";

            string classes = @"  setosa     setosa     setosa     setosa     setosa     setosa     setosa     setosa    
                         setosa     setosa     setosa     setosa     setosa     setosa     setosa     setosa    
                         setosa     setosa     setosa     setosa     setosa     setosa     setosa     setosa    
                         setosa     setosa     setosa     setosa     setosa     setosa     setosa     setosa    
                         setosa     setosa     setosa     setosa     setosa     setosa     setosa     setosa    
                         setosa     setosa     setosa     setosa     setosa     setosa     setosa     setosa    
                         setosa     setosa     versicolor versicolor versicolor versicolor versicolor versicolor
                         versicolor versicolor versicolor versicolor versicolor versicolor versicolor versicolor
                         versicolor versicolor versicolor versicolor versicolor versicolor versicolor versicolor
                         versicolor versicolor versicolor versicolor versicolor versicolor versicolor versicolor
                         versicolor versicolor versicolor versicolor versicolor versicolor versicolor versicolor
                         versicolor versicolor versicolor versicolor versicolor versicolor versicolor versicolor
                         versicolor versicolor versicolor versicolor virginica  virginica  virginica  virginica 
                        virginica  virginica  virginica  virginica  virginica  virginica  virginica  virginica 
                        virginica  virginica  virginica  virginica  virginica  virginica  virginica  virginica 
                        virginica  virginica  virginica  virginica  virginica  virginica  virginica  virginica 
                        virginica  virginica  virginica  virginica  virginica  virginica  virginica  virginica 
                        virginica  virginica  virginica  virginica  virginica  virginica  virginica  virginica 
                        virginica  virginica  virginica  virginica  virginica  virginica";

            string[] splits = { " ", "\t", "\n", "\r" };
            List<string> csvDataList = csvData.Split(splits, StringSplitOptions.RemoveEmptyEntries).ToList();
            List<string> classesList = classes.Split(splits, StringSplitOptions.RemoveEmptyEntries).ToList();

            // remove the first 4 header items
            csvDataList.RemoveRange(0, 4);

            int finalRows = csvDataList.Count / 5;
            irisVariates = new double[finalRows, 4];
            irisClasses = new int[finalRows];

            int index = 0;
            for (int i = 0; i < finalRows; i++) // Rows
            {
                index++; // skip past the row number
                for (int j = 0; j < 4; j++) // Columns
                {
                    irisVariates[i, j] = Convert.ToDouble(csvDataList[index]);
                    index++;
                }
            }

            for (int i = 0; i < finalRows; i++)
            {
                if (classesList[i] == "setosa") irisClasses[i] = 1;
                else if (classesList[i] == "versicolor") irisClasses[i] = 2;
                else if (classesList[i] == "virginica") irisClasses[i] = 3;
            }

            irisVariatesTest = new double[1, 4];
            irisVariatesTest[0, 0] = 6.3;
            irisVariatesTest[0, 1] = 2.5;
            irisVariatesTest[0, 2] = 5.0;
            irisVariatesTest[0, 3] = 1.9;

            irisVariatesReordered = new double[finalRows, 4];
            for (int i = 0; i < finalRows; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    irisVariatesReordered[i, j] = irisVariates[i, 4 - 1 - j];
                }
            }
            irisVariatesTestReordered = new double[1, 4];
            for (int j = 0; j < 4; j++)
            {
                irisVariatesTestReordered[0, j] = irisVariatesTest[0, 4 - 1 - j];
            }

        }
    }
}
