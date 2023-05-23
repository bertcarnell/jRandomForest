/*
 * Copyright 2013 Rob Carnell
 * 
 */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace CsUtility
{
    /// <summary>
    /// Contains assertion types that are not provided with the standard MSTest assertions.   
    /// </summary>   
    public static class ExceptionAssert
    {   
        /// <summary>   
        /// Checks to make sure that the input delegate throws a exception of type TException.   
        /// </summary>   
        /// <typeparam name="TException">The type of exception expected.</typeparam>  
        /// <param name="blockToExecute">The block of code to execute to generate the exception.</param>  
        public static void Throws<TException>(Action blockToExecute) where TException : System.Exception  
        {  
            try  
            {
                blockToExecute();
            }
            catch (Exception ex)
            {
                Assert.IsTrue(ex.GetType() == typeof(TException), "Expected exception of type " + typeof(TException) + " but type of " + ex.GetType() + " was thrown instead.");
                return;
            }
            Assert.Fail("Expected exception of type " + typeof(TException) + " but no exception was thrown.");
        }
        
        /// <summary>
        /// Checks to make sure that the input delegate throws a exception of type TException.
        /// </summary>
        /// <typeparam name="TException">The type of exception expected.</typeparam>
        /// <param name="blockToExecute">The block of code to execute to generate the exception.</param>
        public static void Throws<TException>(string expectedMessage, Action blockToExecute) where TException : System.Exception
        {
            try
            {
                blockToExecute();
            }
            catch (Exception ex)
            {
                Assert.IsTrue(ex.GetType() == typeof(TException), "Expected exception of type " + typeof(TException) + " but type of " + ex.GetType() + " was thrown instead.");
                Assert.AreEqual(expectedMessage, ex.Message, "Expected exception with a message of '" + expectedMessage + "' but exception with message of '" + ex.Message + "' was thrown instead.");
                return;
            }
            Assert.Fail("Expected exception of type " + typeof(TException) + " but no exception was thrown.");
        }

        public static void NoThrow(Action blockToExecute)
        {
            try
            {
                blockToExecute();
            }
            catch (Exception ex)
            {
                Assert.Fail("Expected no exception, but exception was thrown: " + ex.Message);
            }
            Assert.IsTrue(true);
        }
    }
}
