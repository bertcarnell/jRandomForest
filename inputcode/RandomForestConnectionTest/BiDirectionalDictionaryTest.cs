/*
 * Copyright 2013 Rob Carnell
 * 
 */

using RandomForestConnection;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Collections;

namespace RandomForestConnectionTest
{
    
    
    /// <summary>
    ///This is a test class for BiDirectionalDictionaryTest and is intended
    ///to contain all BiDirectionalDictionaryTest Unit Tests
    ///</summary>
    [TestClass()]
    public class BiDirectionalDictionaryTest
    {


        private TestContext testContextInstance;
        private BiDirectionalDictionary<string, int> bdd_empty, bdd_four;
        private KeyValuePair<string, int> kvpe, kvpf;

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
        [TestInitialize()]
        public void MyTestInitialize()
        {
            bdd_empty = new BiDirectionalDictionary<string, int>();
            bdd_four = new BiDirectionalDictionary<string,int>(4);
            bdd_four.Add("a", 1);
            bdd_four.Add("b", 2);
            bdd_four.Add("c", 3);
            bdd_four.Add("d", 4);
            kvpe = new KeyValuePair<string, int>("e", 5);
            kvpf = new KeyValuePair<string, int>("f", 6);
        }
        //
        //Use TestCleanup to run code after each test has run
        //[TestCleanup()]
        //public void MyTestCleanup()
        //{
        //}
        //
        #endregion


        /// <summary>
        ///A test for BiDirectionalDictionary`2 Constructor
        ///</summary>
        public void BiDirectionalDictionaryConstructorTestHelper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>(8);
            Assert.AreEqual(0, target.Count);
            Assert.IsNotNull(target);
        }

        [TestMethod()]
        public void BiDirectionalDictionaryConstructorTest()
        {
            BiDirectionalDictionaryConstructorTestHelper<GenericParameterHelper, GenericParameterHelper>();
        }

        /// <summary>
        ///A test for BiDirectionalDictionary`2 Constructor
        ///</summary>
        public void BiDirectionalDictionaryConstructorTest1Helper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>();
            Assert.AreEqual(0, target.Count);
            Assert.IsNotNull(target);
        }

        [TestMethod()]
        public void BiDirectionalDictionaryConstructorTest1()
        {
            BiDirectionalDictionaryConstructorTest1Helper<GenericParameterHelper, GenericParameterHelper>();
        }

        /// <summary>
        ///A test for Add
        ///</summary>
        public void AddTestHelper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>();
            TKey key = default(TKey);
            TValue val = default(TValue);
            KeyValuePair<TKey, TValue> kvp = new KeyValuePair<TKey, TValue>(key, val);
            target.Add(kvp);
            Assert.AreEqual(1, target.Count);
        }

        [TestMethod()]
        public void AddTest()
        {
            AddTestHelper<int, int>();
            AddTestHelper<int, double>();

            bdd_empty.Add(kvpe);
            Assert.AreEqual(1, bdd_empty.Count);
            Assert.AreEqual(kvpe.Value, bdd_empty[kvpe.Key]);
            Assert.AreEqual(kvpe.Key, bdd_empty.getKeyWithValue(kvpe.Value));
            bdd_four.Add(kvpf);
            Assert.AreEqual(5, bdd_four.Count);
        }

        /// <summary>
        ///A test for Add
        ///</summary>
        public void AddTest1Helper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>();
            TKey key = default(TKey);
            TValue val = default(TValue);
            target.Add(key, val);
            Assert.AreEqual(1, target.Count);
            Assert.AreEqual(val, target[key]);
        }

        [TestMethod()]
        public void AddTest1()
        {
            AddTest1Helper<int, int>();
            AddTest1Helper<int, double>();

            bdd_empty.Add("h", 11);
            Assert.AreEqual(11, bdd_empty["h"]);

            // test that you are prevented from add another pair with the same key
            ExceptionAssert.Throws<ArgumentException>(() => bdd_empty.Add("h", 10));

            // test that you are prevented from adding nother pair with the same value
            ExceptionAssert.Throws<ArgumentException>(() => bdd_empty.Add("h", 11));

        }

        /// <summary>
        ///A test for BackwardContainsKey
        ///</summary>
        public void BackwardContainsKeyTestHelper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>(1);
            target.Add(default(TKey), default(TValue));
            Assert.IsTrue(target.BackwardContainsKey(default(TValue)));
        }

        [TestMethod()]
        public void BackwardContainsKeyTest()
        {
            BackwardContainsKeyTestHelper<int, int>();
            BackwardContainsKeyTestHelper<int, double>();

            Assert.IsFalse(bdd_empty.BackwardContainsKey(1));
            Assert.IsTrue(bdd_four.BackwardContainsKey(3));
            Assert.IsFalse(bdd_four.BackwardContainsKey(101));
        }

        /// <summary>
        ///A test for BackwardTryGetValue
        ///</summary>
        public void BackwardTryGetValueTestHelper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>();
        }

        [TestMethod()]
        public void BackwardTryGetValueTest()
        {
            BackwardTryGetValueTestHelper<GenericParameterHelper, GenericParameterHelper>();

            string temp;
            Assert.IsTrue(bdd_four.BackwardTryGetValue(1, out temp));
            Assert.AreEqual("a", temp);
            Assert.IsFalse(bdd_four.BackwardTryGetValue(1001, out temp));
        }

        /// <summary>
        ///A test for BackwardsCopyTo
        ///</summary>
        public void BackwardsCopyToTestHelper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>(1);
            target.Add(default(TKey), default(TValue));
            KeyValuePair<TValue, TKey>[] kvp_array = new KeyValuePair<TValue, TKey>[1];
            int arrayIndex = 0;
            target.BackwardsCopyTo(kvp_array, arrayIndex);
            Assert.AreEqual(default(TValue), kvp_array[0].Key);
            Assert.AreEqual(default(TKey), kvp_array[0].Value);
        }

        [TestMethod()]
        public void BackwardsCopyToTest()
        {
            BackwardsCopyToTestHelper<int, int>();
            BackwardsCopyToTestHelper<int, Double>();
            BackwardsCopyToTestHelper<double, int>();

            KeyValuePair<int, string>[] kvp_array = new KeyValuePair<int, string>[6];
            bdd_four.BackwardsCopyTo(kvp_array, 2);
            Assert.AreEqual(bdd_four["a"], kvp_array[2].Key);

            // fail when you copy too far
            ExceptionAssert.Throws<IndexOutOfRangeException>(() => bdd_four.BackwardsCopyTo(kvp_array, 5));
        }

        /// <summary>
        ///A test for Clear
        ///</summary>
        public void ClearTestHelper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>(10);
            target.Add(default(TKey), default(TValue));
            Assert.AreEqual(1, target.Count);
            target.Clear();
            Assert.AreEqual(0, target.Count);
        }

        [TestMethod()]
        public void ClearTest()
        {
            ClearTestHelper<int, double>();

            Assert.AreEqual(4, bdd_four.Count);
            bdd_four.Clear();
            Assert.AreEqual(0, bdd_four.Count);
        }

        /// <summary>
        ///A test for Contains
        ///</summary>
        public void ContainsTestHelper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>(1);
            KeyValuePair<TKey, TValue> kvp = new KeyValuePair<TKey, TValue>(default(TKey), default(TValue));
            target.Add(kvp);
            Assert.IsTrue(target.Contains(kvp));
        }

        [TestMethod()]
        public void ContainsTest()
        {
            ContainsTestHelper<int, int>();
            ContainsTestHelper<int, double>();

            Assert.IsTrue(bdd_four.Contains(new KeyValuePair<string, int>("a", 1)));
            Assert.IsFalse(bdd_empty.Contains(kvpe));
        }

        /// <summary>
        ///A test for ContainsKey
        ///</summary>
        public void ContainsKeyTestHelper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>(1);
            target.Add(default(TKey), default(TValue));
            Assert.IsTrue(target.ContainsKey(default(TKey)));
        }

        [TestMethod()]
        public void ContainsKeyTest()
        {
            ContainsKeyTestHelper<int, int>();
            ContainsKeyTestHelper<double, int>();

            Assert.IsFalse(bdd_empty.ContainsKey("a"));
            Assert.IsTrue(bdd_four.ContainsKey("a"));
            Assert.IsFalse(bdd_four.ContainsKey("z"));
        }

        /// <summary>
        ///A test for CopyTo
        ///</summary>
        public void CopyToTestHelper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>(1);
            target.Add(default(TKey), default(TValue));
            KeyValuePair<TKey, TValue>[] kvp_array = new KeyValuePair<TKey,TValue>[1];
            int arrayIndex = 0;
            target.CopyTo(kvp_array, arrayIndex);
            Assert.AreEqual(default(TKey), kvp_array[0].Key);
            Assert.AreEqual(default(TValue), kvp_array[0].Value);
        }

        [TestMethod()]
        public void CopyToTest()
        {
            CopyToTestHelper<int, int>();
            CopyToTestHelper<Double, int>();

            KeyValuePair<string, int>[] kvp_array = new KeyValuePair<string, int>[6];
            bdd_four.CopyTo(kvp_array, 2);
            Assert.AreEqual(bdd_four["a"], kvp_array[2].Value);

            // fail when you copy too far
            ExceptionAssert.Throws<IndexOutOfRangeException>(() => bdd_four.CopyTo(kvp_array, 5));
        }

        /// <summary>
        ///A test for GetBackwardEnumerator
        ///</summary>
        public void GetBackwardEnumeratorTestHelper<TKey, TValue>()
        {
            /*BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>();
            TKey key = default(TKey);
            TValue val = default(TValue);
            target.Add(key, val);
            Dictionary<TValue, TKey>.Enumerator actual;
            actual = target.GetBackwardEnumerator();
            Assert.IsNotNull(actual);
            Assert.IsTrue(actual.MoveNext());
            Assert.AreEqual(key, actual.Current.Value);
            Assert.AreEqual(val, actual.Current.Key);*/
        }

        [TestMethod()]
        public void GetBackwardEnumeratorTest()
        {
            GetBackwardEnumeratorTestHelper<GenericParameterHelper, GenericParameterHelper>();

            Assert.IsNotNull(bdd_four.GetBackwardEnumerator());
        }

        /// <summary>
        ///A test for GetEnumerator
        ///</summary>
        public void GetEnumeratorTestHelper<TKey, TValue>()
        {
            /*BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>();
            TKey key = default(TKey);
            TValue val = default(TValue);
            target.Add(key, val);
            Dictionary<TKey, TValue>.Enumerator actual;
            actual = target.GetEnumerator();
            Assert.IsNotNull(actual);
            Assert.IsTrue(actual.MoveNext());
            Assert.AreEqual(val, actual.Current.Value);
            Assert.AreEqual(key, actual.Current.Key);*/
        }

        [TestMethod()]
        public void GetEnumeratorTest()
        {
            GetEnumeratorTestHelper<GenericParameterHelper, GenericParameterHelper>();

            Assert.IsNotNull(bdd_four.GetEnumerator());
            IEnumerator ie = bdd_four.GetEnumerator();
            int count = 0;
            while (ie.MoveNext())
                count++;
            Assert.AreEqual(4, count);
            ie = bdd_four.GetEnumerator();
            ie.MoveNext();
            kvpe = (KeyValuePair<string, int>) ie.Current;
            Assert.AreEqual(1, kvpe.Value);
            ie.MoveNext();
            kvpe = (KeyValuePair<string, int>)ie.Current;
            Assert.AreEqual(2, kvpe.Value);
        }

        /// <summary>
        ///A test for Remove
        ///</summary>
        public void RemoveTestHelper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>();
        }

        [TestMethod()]
        public void RemoveTest()
        {
            RemoveTestHelper<GenericParameterHelper, GenericParameterHelper>();

            bdd_four.Add(kvpe);
            bdd_four.Add(kvpf);
            Assert.AreEqual(6, bdd_four.Count);
            Assert.IsTrue(bdd_four.Contains(kvpe));
            bdd_four.Remove(kvpe);
            Assert.AreEqual(5, bdd_four.Count);
            Assert.IsFalse(bdd_four.Contains(kvpe));
        }

        /// <summary>
        ///A test for Remove
        ///</summary>
        public void RemoveTest1Helper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>();
        }

        [TestMethod()]
        public void RemoveTest1()
        {
            RemoveTest1Helper<GenericParameterHelper, GenericParameterHelper>();

            bdd_four.Add(kvpe);
            bdd_four.Add(kvpf);
            Assert.AreEqual(6, bdd_four.Count);
            Assert.AreEqual(kvpe.Value, bdd_four[kvpe.Key]);
            bdd_four.Remove(kvpe.Key);
            Assert.AreEqual(5, bdd_four.Count);
            Assert.IsFalse(bdd_four.ContainsKey(kvpe.Key));
        }

        /// <summary>
        ///A test for System.Collections.Generic.IEnumerable<System.Collections.Generic.KeyValuePair<TKey,TValue>>.GetEnumerator
        ///</summary>
        public void GetEnumeratorTest1Helper<TKey, TValue>()
        {
            /*IEnumerable<KeyValuePair<TKey, TValue>> target = new BiDirectionalDictionary<TKey, TValue>(); // TODO: Initialize to an appropriate value
            IEnumerator<KeyValuePair<KeyValuePair<TKey, TValue>, TValue>> expected = null; // TODO: Initialize to an appropriate value
            IEnumerator<KeyValuePair<TKey, TValue>> actual;
            actual = target.GetEnumerator();
            Assert.AreEqual(expected, actual);*/
        }

        [TestMethod()]
        [DeploymentItem("RandomForestConnection.dll")]
        public void GetEnumeratorTest1()
        {
            GetEnumeratorTest1Helper<GenericParameterHelper, GenericParameterHelper>();

            Assert.IsNotNull(bdd_four.GetEnumerator());
        }

        /// <summary>
        ///A test for System.Collections.IEnumerable.GetEnumerator
        ///</summary>
        public void GetEnumeratorTest2Helper<TKey, TValue>()
        {
            /*IEnumerable target = new BiDirectionalDictionary<TKey, TValue>(); // TODO: Initialize to an appropriate value
            IEnumerator expected = null; // TODO: Initialize to an appropriate value
            IEnumerator actual;
            actual = target.GetEnumerator();
            Assert.AreEqual(expected, actual);*/
        }

        [TestMethod()]
        [DeploymentItem("RandomForestConnection.dll")]
        public void GetEnumeratorTest2()
        {
            GetEnumeratorTest2Helper<GenericParameterHelper, GenericParameterHelper>();

            Assert.IsNotNull(bdd_four.GetEnumerator());
        }

        /// <summary>
        ///A test for TryGetValue
        ///</summary>
        public void TryGetValueTestHelper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>();
            TKey key = default(TKey);
            TValue value = default(TValue);
            Assert.IsFalse(target.TryGetValue(key, out value));
        }

        [TestMethod()]
        public void TryGetValueTest()
        {
            TryGetValueTestHelper<int, int>();
            TryGetValueTestHelper<int, Double>();
            TryGetValueTestHelper<double, int>();

            int outValue;

            Assert.IsTrue(bdd_four.TryGetValue("a", out outValue));
            Assert.AreEqual(1, outValue);
            Assert.IsFalse(bdd_four.TryGetValue("z", out outValue));
        }

        /// <summary>
        ///A test for getKeyWithValue
        ///</summary>
        public void getKeyWithValueTestHelper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>();
        }

        [TestMethod()]
        public void getKeyWithValueTest()
        {
            getKeyWithValueTestHelper<GenericParameterHelper, GenericParameterHelper>();

            Assert.AreEqual("a", bdd_four.getKeyWithValue(1));
            Assert.AreEqual("b", bdd_four.getKeyWithValue(2));
            Assert.AreEqual("c", bdd_four.getKeyWithValue(3));
            Assert.AreEqual("d", bdd_four.getKeyWithValue(4));
        }

        /// <summary>
        ///A test for getKeyWithValue
        ///</summary>
        public void getKeyWithValueTest1Helper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>();
        }

        [TestMethod()]
        public void getKeyWithValueTest1()
        {
            getKeyWithValueTest1Helper<GenericParameterHelper, GenericParameterHelper>();

            int[] expectedVals = { 1, 2, 3, 4 };
            string[] expectedKeys = { "a", "b", "c", "d" };

            Assert.AreEqual(expectedKeys.Length, bdd_four.getKeyWithValue(expectedVals).Length);
            for (int i = 0; i < expectedKeys.Length; i++)
                Assert.AreEqual(expectedKeys[i], bdd_four.getKeyWithValue(expectedVals)[i]);

            expectedVals = null;
            Assert.IsNull(bdd_four.getKeyWithValue(expectedVals));
        }

        /// <summary>
        ///A test for getValueWithKey
        ///</summary>
        public void getValueWithKeyTestHelper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>();
        }

        [TestMethod()]
        public void getValueWithKeyTest()
        {
            getValueWithKeyTestHelper<GenericParameterHelper, GenericParameterHelper>();

            int[] expectedVals = { 1, 2, 3, 4 };
            string[] expectedKeys = { "a", "b", "c", "d" };

            for (int i = 0; i < expectedVals.Length; i++)
                Assert.AreEqual(expectedVals[i], bdd_four.getValueWithKey(expectedKeys)[i]);

            expectedKeys = null;
            Assert.IsNull(bdd_four.getValueWithKey(expectedKeys));
        }

        /// <summary>
        ///A test for BackwardKeys
        ///</summary>
        public void BackwardKeysTestHelper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>();
            target.Add(default(TKey), default(TValue));
            Assert.AreEqual(1, target.BackwardKeys.Count);
            Assert.IsTrue(target.BackwardKeys.Contains(default(TValue)));
        }

        [TestMethod()]
        public void BackwardKeysTest()
        {
            BackwardKeysTestHelper<int, int>();
            BackwardKeysTestHelper<int, double>();

            Assert.AreEqual(4, bdd_four.BackwardKeys.Count);
            Assert.IsTrue(bdd_four.BackwardKeys.Contains(1));
            Assert.IsTrue(bdd_four.BackwardKeys.Contains(2));
            Assert.IsTrue(bdd_four.BackwardKeys.Contains(3));
            Assert.IsTrue(bdd_four.BackwardKeys.Contains(4));
            Assert.IsFalse(bdd_four.BackwardKeys.Contains(5));
        }

        /// <summary>
        ///A test for BackwardValues
        ///</summary>
        public void BackwardValuesTestHelper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>();
        }

        [TestMethod()]
        public void BackwardValuesTest()
        {
            BackwardValuesTestHelper<GenericParameterHelper, GenericParameterHelper>();

            string[] expectedKeys = { "a", "b", "c", "d" };
            for (int i = 0; i < expectedKeys.Length; i++)
                Assert.IsTrue(bdd_four.BackwardValues.Contains(expectedKeys[i]));
        }

        /// <summary>
        ///A test for Count
        ///</summary>
        public void CountTestHelper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>(10);
            Assert.AreEqual(0, target.Count);
        }

        [TestMethod()]
        public void CountTest()
        {
            CountTestHelper<GenericParameterHelper, GenericParameterHelper>();

            Assert.AreEqual(4, bdd_four.Count);
            Assert.AreEqual(0, bdd_empty.Count);
        }

        /// <summary>
        ///A test for IsReadOnly
        ///</summary>
        public void IsReadOnlyTestHelper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>();
            Assert.IsFalse(target.IsReadOnly);
        }

        [TestMethod()]
        public void IsReadOnlyTest()
        {
            IsReadOnlyTestHelper<GenericParameterHelper, GenericParameterHelper>();
        }

        /// <summary>
        ///A test for Item
        ///</summary>
        public void ItemTestHelper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>();
            target.Add(default(TKey), default(TValue));
            TKey key = default(TKey);
            TValue expected = default(TValue);
            TValue actual;
            target[key] = expected;
            actual = target[key];
            Assert.AreEqual(expected, actual);
        }

        [TestMethod()]
        public void ItemTest()
        {
            ItemTestHelper<int, int>();
            ItemTestHelper<int, Double>();
            ItemTestHelper<double, int>();

            int test = bdd_four["a"];
            Assert.AreEqual(1, test);
            bdd_four["a"] = 7;
            test = bdd_four["a"];
            Assert.AreEqual(7, test);
        }

        /// <summary>
        ///A test for Keys
        ///</summary>
        public void KeysTestHelper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>();
        }

        [TestMethod()]
        public void KeysTest()
        {
            KeysTestHelper<GenericParameterHelper, GenericParameterHelper>();

            string[] expectedKeys = { "a", "b", "c", "d" };

            Assert.AreEqual(4, bdd_four.Keys.Count);
            for (int i = 0; i < expectedKeys.Length; i++)
                Assert.IsTrue(bdd_four.Keys.Contains(expectedKeys[i]));
        }

        /// <summary>
        ///A test for Values
        ///</summary>
        public void ValuesTestHelper<TKey, TValue>()
        {
            BiDirectionalDictionary<TKey, TValue> target = new BiDirectionalDictionary<TKey, TValue>();
        }

        [TestMethod()]
        public void ValuesTest()
        {
            ValuesTestHelper<GenericParameterHelper, GenericParameterHelper>();

            int[] expectedVals = { 1, 2, 3, 4 };

            for (int i = 0; i < expectedVals.Length; i++)
                Assert.IsTrue(bdd_four.Values.Contains(expectedVals[i]));
        }
    }
}
