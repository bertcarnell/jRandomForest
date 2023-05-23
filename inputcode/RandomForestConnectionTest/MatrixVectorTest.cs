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
    ///This is a test class for MatrixVectorTest and is intended
    ///to contain all MatrixVectorTest Unit Tests
    ///</summary>
    [TestClass()]
    public class MatrixVectorTest
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
        ///A test for MatrixVector`1 Constructor
        ///</summary>
        public void MatrixVectorConstructorTestHelper<T>()
        {
            int nrows = 4; 
            int ncols = 3; 
            MatrixVector<T> target = new MatrixVector<T>(nrows, ncols);
            Assert.AreEqual(nrows, target.nrow);
            Assert.AreEqual(ncols, target.ncol);
            Assert.AreEqual(12, target.values.Length);

            ExceptionAssert.Throws<ArgumentException>(() => new MatrixVector<T>(3, 0));
            ExceptionAssert.Throws<ArgumentException>(() => new MatrixVector<T>(0, 4));
            ExceptionAssert.Throws<ArgumentException>(() => new MatrixVector<T>(0, 0));
        }

        [TestMethod()]
        public void MatrixVectorConstructorTest()
        {
            MatrixVectorConstructorTestHelper<GenericParameterHelper>();
            MatrixVectorConstructorTestHelper<double>();
            MatrixVectorConstructorTestHelper<int>();
        }

        /// <summary>
        ///A test for MatrixVector`1 Constructor
        ///</summary>
        [TestMethod()]
        public void MatrixVectorConstructorTest1()
        {
            int[,] expected = { { 1, 2, 3 }, { 4, 5, 6 } };
            MatrixVector<int> target = new MatrixVector<int>(expected);
            Assert.AreEqual(expected.GetLength(0), target.nrow);
            Assert.AreEqual(expected.GetLength(1), target.ncol);
            Assert.AreEqual(expected.GetLength(0) * expected.GetLength(1), target.values.Length);
            for (int i = 0; i < expected.GetLength(0); i++)
            {
                for (int j = 0; j < expected.GetLength(1); j++)
                {
                    Assert.AreEqual(expected[i, j], target[i, j]);
                }
            }
        }

        /// <summary>
        ///A test for fill
        ///</summary>
        public void fillTestHelper<T>()
        {
            MatrixVector<T> target = new MatrixVector<T>(3, 4);
            T x = default(T);
            target.fill(x);
            for (int i = 0; i < target.values.Length; i++)
                Assert.AreEqual(default(T), target.values[i]);
        }

        [TestMethod()]
        public void fillTest()
        {
            fillTestHelper<GenericParameterHelper>();
            fillTestHelper<double>();
            fillTestHelper<int>();

            MatrixVector<double> target = new MatrixVector<double>(2, 2);
            double expected = 5.0;
            target.fill(expected);
            foreach (double d in target.values)
                Assert.AreEqual(expected, d);
        }

        /// <summary>
        ///A test for getMatrix
        ///</summary>
        public void getMatrixTestHelper<T>()
        {
            int nrow = 3;
            int ncol = 4;
            MatrixVector<T> target = new MatrixVector<T>(nrow, ncol);
            target.fill(default(T));
            Assert.AreEqual(nrow, target.getMatrix().GetLength(0));
            Assert.AreEqual(ncol, target.getMatrix().GetLength(1));
            for (int i = 0; i < nrow; i++)
            {
                for (int j = 0; j < ncol; j++)
                {
                    Assert.AreEqual(default(T), target.getMatrix()[i, j]);
                }
            }
            IEnumerator ie = target.getMatrix().GetEnumerator();
            while (ie.MoveNext())
                Assert.AreEqual(default(T), (T) ie.Current);
        }

        [TestMethod()]
        public void getMatrixTest()
        {
            getMatrixTestHelper<GenericParameterHelper>();
            getMatrixTestHelper<double>();
            getMatrixTestHelper<int>();

            int nrow = 10;
            int ncol = 2;
            double expected = 1.2;
            MatrixVector<double> target = new MatrixVector<double>(nrow, ncol);
            target.fill(expected);
            double[,] actual = target.getMatrix();
            Assert.AreEqual(nrow, actual.GetLength(0));
            Assert.AreEqual(ncol, actual.GetLength(1));
            IEnumerator i = actual.GetEnumerator();
            while (i.MoveNext())
            {
                Assert.AreEqual(expected, Convert.ToDouble(i.Current));
            }

            MatrixVector<int> test = new MatrixVector<int>(2, 3, 4);
            ExceptionAssert.Throws<InvalidCastException>(() => test.getMatrix());
        }

        /// <summary>
        ///A test for getVector
        ///</summary>
        public void getVectorTestHelper<T>()
        {
            int nrow = 2;
            int ncol = 3;
            MatrixVector<T> target = new MatrixVector<T>(nrow, ncol);
            T[] actual = target.getVector();
            Assert.AreEqual(nrow*ncol, target.getVector().Length);
        }

        [TestMethod()]
        public void getVectorTest()
        {
            getVectorTestHelper<GenericParameterHelper>();
            getVectorTestHelper<double>();
            getVectorTestHelper<int>();
        }

        /// <summary>
        ///A test for transpose
        ///</summary>
        [TestMethod()]
        public void transposeTest()
        {
            // 1 2 3 4 5 6 7 8 9 10 11 12
            ///////////
            // 1 5 9
            // 2 6 10
            // 3 7 11
            // 4 8 12
            ///////////
            // 1 2 3 4
            // 5 6 7 8
            // 9 10 11 12
            /////////////
            // 1 5 9 2 6 10 3 7 11 4 8 12
            int[,] orig = { { 1, 5, 9 }, { 2, 6, 10 }, { 3, 7, 11 }, { 4, 8, 12 } };
            int[] expected1 = {1,2,3,4,5,6,7,8,9,10,11,12};
            MatrixVector<int> test = new MatrixVector<int>(orig);
            for (int i = 0; i < test.values.Length; i++)
                Assert.AreEqual(expected1[i], test[i]);
            test.transpose();
            int[] expected2 = { 1, 5, 9, 2, 6, 10, 3, 7, 11, 4, 8, 12 };
            for (int i = 0; i < test.values.Length; i++)
                Assert.AreEqual(expected2[i], test[i]);
            for (int i = 0; i < orig.GetLength(0); i++)
            {
                for (int j = 0; j < orig.GetLength(1); j++)
                {
                    Assert.AreEqual(orig[i, j], test[j, i]);
                }
            }
            Assert.AreEqual(orig.GetLength(0), test.ncol);
            Assert.AreEqual(orig.GetLength(1), test.nrow);

            MatrixVector<int> test2 = new MatrixVector<int>(1, 2, 3);
            test2.fill(0);
            ExceptionAssert.Throws<InvalidCastException>(() => test2.transpose());
        }

        /// <summary>
        ///A test for []
        ///</summary>
        [TestMethod()]
        public void ItemTest()
        {
            int[,] orig = { { 1, 5, 9 }, { 2, 6, 10 }, { 3, 7, 11 }, { 4, 8, 12 } };
            MatrixVector<int> test = new MatrixVector<int>(orig);
            Assert.AreEqual(orig[0, 0], test[0]);
            Assert.AreEqual(orig[2, 2], test[10]);

            int[] expected = new MatrixVector<int>(orig).values;
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], test[i]);
        }

        /// <summary>
        ///A test for Item
        ///</summary>
        [TestMethod()]
        public void ItemTest1()
        {
            int[,] orig = { { 1, 5, 9 }, { 2, 6, 10 }, { 3, 7, 11 }, { 4, 8, 12 } };
            MatrixVector<int> test = new MatrixVector<int>(orig);

            for (int i = 0; i < orig.GetLength(0); i++)
            {
                for (int j = 0; j < orig.GetLength(1); j++)
                {
                    Assert.AreEqual(orig[i,j], test[i,j]);
                }
            }
        }

        /// <summary>
        ///A test for rowSum
        ///</summary>
        [TestMethod()]
        public void rowSumTest()
        {
            int nrow = 2;
            int ncol = 3;
            double value = 3.3;
            MatrixVector<double> test = new MatrixVector<double>(nrow, ncol);
            test.fill(value);
            double[] actual = test.rowSum();

            foreach (double d in actual)
                Assert.AreEqual(ncol * value, d);
            Assert.AreEqual(nrow, actual.Length);

            int[,] orig = { { 1, 5, 9 }, { 2, 6, 10 }, { 3, 7, 11 }, { 4, 8, 12 } };
            double[] expected = { 15.0, 18.0, 21.0, 24.0 };
            MatrixVector<int> test2 = new MatrixVector<int>(orig);
            double[] actual2 = test2.rowSum();
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], actual2[i]);
        }

        /// <summary>
        ///A test for normalizeByRow
        ///</summary>
        [TestMethod()]
        public void normalizeByRowTest()
        {
            double[,] orig = { { 1.0, 5.0, 9.0 }, { 2.0, 6.0, 10.0 }, { 3.0, 7.0, 11.0 }, { 4.0, 8.0, 12.0 } };
            double[,] expected = { { 1.0 / 15.0, 5.0 / 15.0, 9.0 / 15.0 }, { 2.0 / 18.0, 6.0 / 18.0, 10.0 / 18.0 }, { 3.0 / 21.0, 7.0 / 21.0, 11.0 / 21.0 }, { 4.0 / 24.0, 8.0 / 24.0, 12.0 / 24.0 } };
            MatrixVector<double> test = new MatrixVector<double>(orig);
            test.normalizeByRow();
            double[,] actual = test.getMatrix();
            IEnumerator itest = actual.GetEnumerator();
            IEnumerator iexpected = expected.GetEnumerator();
            while (itest.MoveNext() && iexpected.MoveNext())
                Assert.AreEqual(Convert.ToDouble(iexpected.Current), Convert.ToDouble(itest.Current), 1E-12);

            // test throws for non-double or float
            MatrixVector<int> test2 = new MatrixVector<int>(2, 3);
            test2.fill(7);
            ExceptionAssert.Throws<Exception>(() => test2.normalizeByRow());

            // exercise where rowsum is zero
            MatrixVector<double> test3 = new MatrixVector<double>(2, 3);
            test3.fill(0.0);
            test3.normalizeByRow();
            foreach (double a in test3.values)
                Assert.AreEqual(0.0, a);
        }

        /// <summary>
        ///A test for MatrixVector`1 Constructor
        ///</summary>
        public void MatrixVectorConstructorTest2Helper<T>()
        {
            int nrows = 4;
            int ncols = 5;
            T[] array = new T[4*5];
            MatrixVector<T> target = new MatrixVector<T>(array, nrows, ncols);
            Assert.AreEqual(nrows, target.nrow);
            Assert.AreEqual(ncols, target.ncol);

            ExceptionAssert.Throws<Exception>(() => new MatrixVector<T>(array, 6, 7));
        }

        [TestMethod()]
        public void MatrixVectorConstructorTest2()
        {
            MatrixVectorConstructorTest2Helper<GenericParameterHelper>();
            MatrixVectorConstructorTest2Helper<int>();
            MatrixVectorConstructorTest2Helper<double>();
            int[] array1 = { 1, 2, 3, 4, 5, 6 };
            int nrow = 2;
            int ncol = 3;
            MatrixVector<int> test = new MatrixVector<int>(array1, nrow, ncol);
            for (int i = 0; i<array1.Length; i++)
                Assert.AreEqual(array1[i], test[i]);
        }

        /// <summary>
        ///A test for MatrixVector`1 Constructor
        ///</summary>
        public void MatrixVectorConstructorTest3Helper<T>()
        {
            int nrows = 3;
            int ncols = 5;
            T fillValue = default(T);
            MatrixVector<T> target = new MatrixVector<T>(nrows, ncols, fillValue);
            foreach (T i in target.getVector())
                Assert.AreEqual(fillValue, i);
        }

        [TestMethod()]
        public void MatrixVectorConstructorTest3()
        {
            MatrixVectorConstructorTest3Helper<GenericParameterHelper>();
            MatrixVectorConstructorTest3Helper<int>();
            MatrixVectorConstructorTest3Helper<double>();
        }

        /// <summary>
        ///A test for permute
        ///</summary>
        [TestMethod()]
        public void permuteTest()
        {            
            /*123456789 10 11 12 13 14 15 ... 24

            1 5 9     13 17 21
            2 6 10    14 18 22
            3 7 11    15 19 23
            4 8 12    16 20 24

            3 1 2 =>

            1  2  3  4      5  6  7  8    9  10 11 12
            13 14 15 16     17 18 19 20   21 22 23 24

            1 13 2 14 3 15 4 16 5 17 6 18 7 19 8 20 9 21 10 22 11 23 12 24*/

            int[] temp = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24 };
            MatrixVector<int> target = new MatrixVector<int>(temp, 4, 3, 2);
            int[] dimensionOrder = {0, 1, 2};
            target.permute(dimensionOrder);
            Assert.AreEqual(4, target.nrow);
            Assert.AreEqual(3, target.ncol);
            Assert.AreEqual(2, target.dim3);
            for (int i = 0; i < temp.Length; i++)
                Assert.AreEqual(temp[i], target.values[i]);

            int[] dimensionOrder2 = {2, 0, 1};
            int[] expected = {1, 13, 2, 14, 3, 15, 4, 16, 5, 17, 6, 18, 7, 19, 8, 20, 9, 21, 10, 22, 11, 23, 12, 24};
            target.permute(dimensionOrder2);
            Assert.AreEqual(2, target.nrow);
            Assert.AreEqual(4, target.ncol);
            Assert.AreEqual(3, target.dim3);
            for (int i = 0; i < temp.Length; i++)
                Assert.AreEqual(expected[i], target.values[i]);

            MatrixVector<int> test2 = new MatrixVector<int>(2, 3);
            test2.fill(1);
            ExceptionAssert.Throws<NotImplementedException>(() => test2.permute(new int[3]));

            MatrixVector<int> test3 = new MatrixVector<int>(2, 3, 4);
            test3.fill(1);
            int[] dimensionOrderFail = { 4, 5, 6 };
            ExceptionAssert.Throws<ArgumentException>(() => test3.permute(dimensionOrderFail));
            dimensionOrderFail = new int[]{ 0, 5, 6 };
            ExceptionAssert.Throws<ArgumentException>(() => test3.permute(dimensionOrderFail));
            dimensionOrderFail = new int[] { 0, 0, 2 };
            ExceptionAssert.Throws<ArgumentException>(() => test3.permute(dimensionOrderFail));
        }

        /// <summary>
        /// Test for the GetColumn method
        /// </summary>
        [TestMethod()]
        public void GetColumnTest()
        {
            int[,] matrix = {{1,2,3},{4,5,6}};
            MatrixVector<int> target = new MatrixVector<int>(matrix);
            int[] actual;
            for (int ci = 0; ci < matrix.GetLength(1); ci++)
            {
                actual = target.GetColumn(ci);
                for (int r = 0; r < matrix.GetLength(0); r++)
                {
                    Assert.AreEqual(matrix[r, ci], actual[r]);
                }
            }
            // Action is a lambda expression with no arguments "()" that uses target.GetColumn(4)
            ExceptionAssert.Throws<ArgumentException>(() => target.GetColumn(4));
        }

        [TestMethod()]
        public void get3dMatrixTest()
        {
            int[] dims = {2, 3, 4};
            MatrixVector<int> test = new MatrixVector<int>(dims[0], dims[1], dims[2]);
            int expected = 2;
            test.fill(expected);
            int[,,] actual = test.get3dMatrix();
            Assert.AreEqual(3, actual.Rank);
            for (int i = 0; i < 3; i++)
            {
                Assert.AreEqual(dims[i], actual.GetLength(i));
            }

            MatrixVector<double> test2 = new MatrixVector<double>(2, 3);
            ExceptionAssert.Throws<InvalidCastException>(() => test2.get3dMatrix());            
        }

        /// <summary>
        ///A test for MatrixVector`1 Constructor
        ///</summary>
        public void MatrixVectorConstructorTest4Helper<T>()
        {
            int x = 4;
            int y = 5;
            int z = 6;
            T[] array = new T[x * y * z];
            MatrixVector<T> target = new MatrixVector<T>(array, x, y, z);
            Assert.AreEqual(x, target.nrow);
            Assert.AreEqual(y, target.ncol);
            Assert.AreEqual(z, target.dim3);

            ExceptionAssert.Throws<Exception>(() => new MatrixVector<T>(array, 2, 2, 2));
        }

        [TestMethod()]
        public void MatrixVectorConstructorTest4()
        {
            MatrixVectorConstructorTest4Helper<GenericParameterHelper>();
            MatrixVectorConstructorTest4Helper<int>();
            MatrixVectorConstructorTest4Helper<double>();
        }

        /// <summary>
        ///A test for MatrixVector`1 Constructor
        ///</summary>
        public void MatrixVectorConstructorTest5Helper<T>()
        {
            int dim1 = 1;
            int dim2 = 2;
            int dim3 = 3;
            MatrixVector<T> target = new MatrixVector<T>(dim1, dim2, dim3);

            Assert.AreEqual(dim1, target.get3dMatrix().GetLength(0));
            Assert.AreEqual(dim2, target.get3dMatrix().GetLength(1));
            Assert.AreEqual(dim3, target.get3dMatrix().GetLength(2));
            Assert.AreEqual(dim1 * dim2 * dim3, target.values.Length);

            ExceptionAssert.Throws<ArgumentException>(() => new MatrixVector<T>(-1, 2, 3));
        }

        [TestMethod()]
        public void MatrixVectorConstructorTest5()
        {
            MatrixVectorConstructorTest5Helper<GenericParameterHelper>();
            MatrixVectorConstructorTest5Helper<int>();
            MatrixVectorConstructorTest5Helper<double>();
        }

        /// <summary>
        ///A test for MatrixVector`1 Constructor
        ///</summary>
        public void MatrixVectorConstructorTest6Helper<T>()
        {
            T[, ,] matrix = new T[2,3,4];
            MatrixVector<T> target = new MatrixVector<T>(matrix);
            Assert.AreEqual(2, target.nrow);
            Assert.AreEqual(3, target.ncol);
            Assert.AreEqual(4, target.dim3);
        }

        [TestMethod()]
        public void MatrixVectorConstructorTest6()
        {
            MatrixVectorConstructorTest6Helper<GenericParameterHelper>();
            MatrixVectorConstructorTest6Helper<int>();
            MatrixVectorConstructorTest6Helper<double>();

            int[, ,] matrix = { { { 1 }, { 2 } }, { { 3 }, { 4 } }, { { 5 }, { 6 } } };
            MatrixVector<int> target = new MatrixVector<int>(matrix);
            Assert.AreEqual(3, target.nrow);
            Assert.AreEqual(2, target.ncol);
            Assert.AreEqual(1, target.dim3);

            Assert.AreEqual(matrix[0, 0, 0], target[0, 0, 0]);
            Assert.AreEqual(matrix[0, 1, 0], target[0, 1, 0]);
            Assert.AreEqual(matrix[1, 0, 0], target[1, 0, 0]);
            Assert.AreEqual(matrix[1, 1, 0], target[1, 1, 0]);
            Assert.AreEqual(matrix[2, 0, 0], target[2, 0, 0]);
            Assert.AreEqual(matrix[2, 1, 0], target[2, 1, 0]);
        }

        [TestMethod()]
        public void GetAndSetTest()
        {
            MatrixVector<int> a = new MatrixVector<int>(2,3);
            a.fill(5);
            // Get
            Assert.AreEqual(5, a[0, 0]);
            // set
            a[1, 1] = 7;
            Assert.AreEqual(7, a[1, 1]);

            MatrixVector<int> b = new MatrixVector<int>(2, 3, 4);
            b.fill(5);
            // Get
            Assert.AreEqual(5, b[0, 0, 0]);
            // set
            b[1, 1, 1] = 7;
            Assert.AreEqual(7, b[1, 1, 1]);

            // as array..
            // Get
            Assert.AreEqual(5, a[0]);
            // set
            a[5] = 7;
            Assert.AreEqual(7, a[5]);
        }
    }
}
