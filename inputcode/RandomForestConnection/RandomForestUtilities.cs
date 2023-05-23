/*
 * Copyright 2013 Rob Carnell
 * 
 */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace RandomForestConnection
{
    /// <summary>
    /// Utilities for the Random Forest class
    /// </summary>
    public static class RandomForestUtilities
    {
        /// <summary>
        /// Create an array and fill it with a value
        /// </summary>
        /// <typeparam name="T">Type</typeparam>
        /// <param name="length">Length of the array</param>
        /// <param name="value">value to fill the array with</param>
        /// <returns>An array of Type T</returns>
        public static T[] CreateAndFill<T>(int length, T value)
        {
            T[] result = new T[length];
            for (int i = 0; i < length; i++)
            {
                result[i] = value;
            }
            return result;
        }
    }
}
