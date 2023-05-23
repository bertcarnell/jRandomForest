using CsRandomForest;
using CsUtility;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace CsRandomForestTest
{
    /// <summary>
    ///This is a test class for CArrayTest and is intended
    ///to contain all CArrayTest Unit Tests
    ///</summary>
    [TestClass()]
    public class CArrayTest
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
        ///A test for CArray Constructor
        ///</summary>
        public void CArrayConstructorTestHelper<T>()
        {
            int length = 5; // TODO: Initialize to an appropriate value
            T[] basea = new T[length];
            CArray<T> target = new CArray<T>(basea);
            Assert.AreEqual(length, target.Length);
        }

        [TestMethod()]
        [DeploymentItem("CsRandomForest.dll")]
        public void CArrayConstructorTest()
        {
            CArrayConstructorTestHelper<GenericParameterHelper>();
            CArrayConstructorTestHelper<int>();
            CArrayConstructorTestHelper<double>();
        }

        /// <summary>
        ///A test for CArray Copy Constructor
        ///</summary>
        [TestMethod()]
        [DeploymentItem("CsRandomForest.dll")]
        public void CArrayConstructorTest1()
        {
            int length = 5;
            int[] basea = new int[length];
            CArray<int> oldCArray = new CArray<int>(basea);
            for (int i = 0; i < length; i++)
                oldCArray[i] = i;
            oldCArray++;
            CArray<int> target = new CArray<int>(oldCArray);

            Assert.AreEqual(1, target.Offset);
            Assert.AreEqual(1, oldCArray.Offset);
            Assert.AreEqual(length, target.Length);
            for (int i = 0; i < length - 1; i++)
            {
                Assert.AreEqual(oldCArray[i], target[i]);
                Assert.AreEqual(oldCArray[i], basea[i+1]);
            }
        }

        /// <summary>
        ///A test for op_Addition
        ///</summary>
        [TestMethod()]
        [DeploymentItem("CsRandomForest.dll")]
        public void op_AdditionTest()
        {
            int length = 4;
            int[] basea = new int[length];
            CArray<int> x = new CArray<int>(basea);
            for (int i = 0; i < length; i++) 
                x[i] = i;
            x += 2; // x = x + 2
            Assert.AreEqual(2, x[0]);
            Assert.AreEqual(3, x[1]);

            int temp;
            ExceptionAssert.Throws<IndexOutOfRangeException>(() => temp = x[2]);

            // test that a new CArray is created in addition and that it is independent of the offset of the original
            CArray<int> y = new CArray<int>(basea);
            CArray<int> z = y + 1;
            Assert.AreEqual(0, basea[0]);
            Assert.AreEqual(1, basea[1]);
            Assert.AreEqual(2, basea[2]);
            Assert.AreEqual(3, basea[3]);
            Assert.AreEqual(0, z[-1]);
            Assert.AreEqual(1, z[0]);
            Assert.AreEqual(2, z[1]);
            Assert.AreEqual(3, z[2]);
            Assert.AreEqual(1, z.Offset);
            Assert.AreEqual(0, y.Offset);

            ExceptionAssert.Throws<IndexOutOfRangeException>(() => temp = z[3]);

            z++;
            Assert.AreEqual(2, z.Offset);
            Assert.AreEqual(0, y.Offset);
        }

        /// <summary>
        ///A test for op_Decrement
        ///</summary>
        [TestMethod()]
        [DeploymentItem("CsRandomForest.dll")]
        public void op_DecrementTest()
        {
            int length = 4;
            int[] basea = new int[length];
            int temp;
            CArray<int> x = new CArray<int>(basea);
            for (int i = 0; i < length; i++)
                x[i] = i;
            x--;
            Assert.AreEqual(0, x[1]);
            Assert.AreEqual(1, x[2]);
            Assert.AreEqual(2, x[3]);
            ExceptionAssert.NoThrow(() => temp = x[4]);
            Assert.AreEqual(3, x[4]);

            ExceptionAssert.Throws<IndexOutOfRangeException>(() => temp = x[5]);
            ExceptionAssert.Throws<IndexOutOfRangeException>(() => temp = x[0]);

            --x;
            Assert.AreEqual(0, x[2]);
            Assert.AreEqual(1, x[3]);

            Assert.AreEqual(0, basea[0]);
            Assert.AreEqual(1, basea[1]);
            Assert.AreEqual(2, basea[2]);
            Assert.AreEqual(3, basea[3]);
        }

        /// <summary>
        ///A test for op_Increment
        ///</summary>
        [TestMethod()]
        [DeploymentItem("CsRandomForest.dll")]
        public void op_IncrementTest()
        {
            int length = 4;
            int[] basea = new int[length];
            int temp;
            CArray<int> x = new CArray<int>(basea);
            for (int i = 0; i < length; i++)
                x[i] = i;
            x++;
            Assert.AreEqual(1, x[0]);
            Assert.AreEqual(2, x[1]);
            Assert.AreEqual(3, x[2]);
            ExceptionAssert.NoThrow(() => temp = x[-1]);
            Assert.AreEqual(0, x[-1]);

            ExceptionAssert.Throws<IndexOutOfRangeException>(() => temp = x[-2]);
            ExceptionAssert.Throws<IndexOutOfRangeException>(() => temp = x[3]);
        }

        /// <summary>
        ///A test for op_Subtraction
        ///</summary>
        [TestMethod()]
        [DeploymentItem("CsRandomForest.dll")]
        public void op_SubtractionTest()
        {
            int length = 4;
            int[] basea = new int[length];
            int temp;
            CArray<int> x = new CArray<int>(basea);
            for (int i = 0; i < length; i++)
                x[i] = i;
            x -= 2;
            Assert.AreEqual(0, x[2]);
            Assert.AreEqual(1, x[3]);
            Assert.AreEqual(2, x[4]);
            ExceptionAssert.NoThrow(() => temp = x[5]);
            Assert.AreEqual(3, x[5]);

            ExceptionAssert.Throws<IndexOutOfRangeException>(() => temp = x[0]);
            ExceptionAssert.Throws<IndexOutOfRangeException>(() => temp = x[6]);
        }

        /// <summary>
        ///A test for Item
        ///</summary>
        [TestMethod()]
        [DeploymentItem("CsRandomForest.dll")]
        public void ItemTest()
        {
            double[] basea = new double[4];
            CArray<double> x = new CArray<double>(basea);
            x[0] = 6.0;
            double test = x[0];
            Assert.AreEqual(6.0, test);
            Assert.AreEqual(6.0, x[0]);
        }

        /// <summary>
        ///A test for Length
        ///</summary>
        public void LengthTestHelper<T>()
        {
            T[] basea = new T[10];
            CArray<T> target = new CArray<T>(basea);
            int actual = target.Length;
            Assert.AreEqual(10, actual);
        }

        [TestMethod()]
        [DeploymentItem("CsRandomForest.dll")]
        public void LengthTest()
        {
            LengthTestHelper<GenericParameterHelper>();
            LengthTestHelper<double>();
            LengthTestHelper<int>();
        }

        /// <summary>
        ///A test for Offset
        ///</summary>
        public void OffsetTestHelper<T>()
        {
            T[] basea = new T[5];
            CArray<T> x = new CArray<T>(basea);
            x--;
            x += 2;
            x++;
            x -= 5;
            Assert.AreEqual(0 - 1 + 2 + 1 - 5, x.Offset);
        }

        [TestMethod()]
        [DeploymentItem("CsRandomForest.dll")]
        public void OffsetTest()
        {
            OffsetTestHelper<GenericParameterHelper>();
            OffsetTestHelper<double>();
            OffsetTestHelper<int>();
        }

        /// <summary>
        ///A test for reset
        ///</summary>
        public void resetTestHelper<T>()
        {
            T[] basea = new T[5];
            CArray<T> x = new CArray<T>(basea);
            x++;
            Assert.AreEqual(1, x.Offset);
            x.reset();
            Assert.AreEqual(0, x.Offset);
        }

        [TestMethod()]
        public void resetTest()
        {
            resetTestHelper<GenericParameterHelper>();
            resetTestHelper<int>();
            resetTestHelper<double>();
        }

        [TestMethod()]
        public void functionTest()
        {
            int[] basea = new int[10];
            CArray<int> y = new CArray<int>(basea);
            for (int i = 0; i < 10; i++)
                y[i] = i * 3;
            simulateFunction(y);
            Assert.AreEqual(-1, y.Offset);
            Assert.AreEqual(10, y[1]);
            Assert.AreEqual(12, y[3]);
            y.reset();
            Assert.AreEqual(10, y[0]);
            Assert.AreEqual(3, y[1]);
            Assert.AreEqual(12, y[2]);
            Assert.AreEqual(9, y[3]);
        }

        public void simulateFunction(CArray<int> x)
        {
            x--;
            x[1] = 10;
            x[3] = 12;
        }

        [TestMethod()]
        public void innerOuterFunctionTest()
        {
            int[] mem = new int[10];
            for (int i = 0; i < 10; i++)
                mem[i] = i;
            CArray<int> x = new CArray<int>(mem);
            for (int i = 0; i < 10; i++)
                Assert.AreEqual(i, x[i]);

            outerFunction(x + 1);
            Assert.AreEqual(0, x.Offset);
            Assert.AreEqual(999, mem[2]);
            Assert.AreEqual(999, x[2]);
            Assert.AreEqual(88, mem[4]);
            Assert.AreEqual(88, x[4]);
            Assert.AreEqual(7, mem[5]);
            Assert.AreEqual(7, x[5]);

            for (int i = 0; i < 10; i++)
                mem[i] = i;
            outerFunction(++x);
            Assert.AreEqual(2, x.Offset); // x got incremented here and passed in to get incremented again
            Assert.AreEqual(999, mem[2]);
            Assert.AreEqual(999, x[0]);
            Assert.AreEqual(88, mem[4]);
            Assert.AreEqual(88, x[2]);
            Assert.AreEqual(7, mem[5]);
            Assert.AreEqual(7, x[3]);

            for (int i = 0; i < 10; i++)
                mem[i] = i;
            outerFunction(x.copy()); // x does not get passed in as a reference here
            Assert.AreEqual(2, x.Offset);
            Assert.AreEqual(999, mem[3]);
            Assert.AreEqual(999, x[1]);
            Assert.AreEqual(88, mem[5]);
            Assert.AreEqual(88, x[3]);
            Assert.AreEqual(7, mem[6]);
            Assert.AreEqual(7, x[4]);
        }

        public void outerFunction(CArray<int> a)
        {
            a++;
            a[0] = 999;
            innerFunction(a + 2);
        }

        public void innerFunction(CArray<int> b)
        {
            b[0] = 88;
            b[1] = 7;
        }
    }
}
