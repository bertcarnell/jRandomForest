using RandomForestConnection;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace RandomForestConnectionTest
{
    
    
    /// <summary>
    ///This is a test class for ArrayPointerTest and is intended
    ///to contain all ArrayPointerTest Unit Tests
    ///</summary>
    [TestClass()]
    public class ArrayPointerTest
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
        ///A test for ArrayPointer`1 Constructor
        ///</summary>
        public void ArrayPointerConstructorTestHelper<T>()
        {
            T expected = default(T);
            ArrayPointer<T> target = expected;
            Assert.AreEqual(expected, target.Value);
            Assert.AreEqual(expected, target.m_ptr[0]);

            ArrayPointer<T> target2 = new ArrayPointer<T>(expected);
            Assert.AreEqual(expected, target2.Value);
            Assert.AreEqual(expected, target2.m_ptr[0]);
        }

        [TestMethod()]
        public void ArrayPointerConstructorTest()
        {
            ArrayPointerConstructorTestHelper<GenericParameterHelper>();
            ArrayPointerConstructorTestHelper<double>();
            ArrayPointerConstructorTestHelper<int>();

            int x = 5;
            ArrayPointer<int> target = x;
            Assert.AreEqual(x, target.Value);
            Assert.AreEqual(x, target.m_ptr[0]);
        }

        /// <summary>
        ///A test for ArrayPointer`1 Constructor
        ///</summary>
        public void ArrayPointerConstructorTest1Helper<T>()
        {
            ArrayPointer<T> target = new ArrayPointer<T>();
            Assert.AreEqual(1, target.m_ptr.Length);
        }

        [TestMethod()]
        public void ArrayPointerConstructorTest1()
        {
            ArrayPointerConstructorTest1Helper<GenericParameterHelper>();
            ArrayPointerConstructorTest1Helper<int>();
            ArrayPointerConstructorTest1Helper<double>();
        }

        /// <summary>
        ///A test for op_Explicit
        ///</summary>
        [TestMethod()]
        public void op_ExplicitTest()
        {
            double expected = 8.0;
            ArrayPointer<int> value = Convert.ToInt32(expected);
            Assert.AreEqual(expected, (double)value);
            Assert.AreEqual(expected, Convert.ToDouble(value));
        }

        /// <summary>
        ///A test for op_Implicit
        ///</summary>
        [TestMethod()]
        public void op_ImplicitTest()
        {
            int expected = 9;
            ArrayPointer<int> value = expected; 
            Assert.AreEqual(expected, value.Value);

            ArrayPointer<int> value2 = 10;  
            Assert.AreEqual(10, value2.Value);
            value2 = 11;
            Assert.AreEqual(11, value2.Value);
            value2 = expected;
            Assert.AreEqual(expected, value2.Value);
        }

        /// <summary>
        ///A test for op_Implicit
        ///</summary>
        [TestMethod()]
        public void op_ImplicitTest1()
        {
            ArrayPointer<int> test = 2;
            // ArrayPointer test is implicitly converted to an int for these tests
            Assert.IsTrue(test == 2);
            Assert.IsTrue(test < 5);
            Assert.IsTrue(test > -1);
            Assert.IsTrue(test <= 5);
            Assert.IsTrue(test >= -1);
            Assert.AreEqual(4, test * 2);
            Assert.AreEqual(2, test.Value); // Note, implicit conversion is not done in Assert.AreEqual
            Assert.AreEqual(1, test - 1);
            Assert.AreEqual(3, test + 1);
            Assert.AreEqual(1, test / 2);
            int x = test;
            Assert.AreEqual(2, x);

        }

        /// <summary>
        ///A test for Value
        ///</summary>
        public void ValueTestHelper<T>()
        {
            ArrayPointer<T> target = new ArrayPointer<T>();
            T expected = default(T); 
            T actual;
            target.Value = expected;
            actual = target.Value;
            Assert.AreEqual(expected, actual);
        }

        [TestMethod()]
        public void ValueTest()
        {
            ValueTestHelper<GenericParameterHelper>();
            ValueTestHelper<int>();
            ValueTestHelper<double>();
        }
    }
}
