/*
 *  Copyright (C) 2012 Rob Carnell
 *  
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *  
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA,
 * or visit http://www.gnu.org/licenses/gpl-2.0.html * 
 */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace CsRandomForest
{
    /// <summary>
    /// A class that operates like a C-style array pointer
    /// </summary>
    /// <typeparam name="T">Type of the Array</typeparam>
    public class CArray<T>
    {
        private T[] m_array;
        private int m_offset;

        /// <summary>
        /// Full length of underlying Array
        /// </summary>
        public int Length { get { return m_array.Length; } }
        /// <summary>
        /// The current offset of the Array pointer
        /// </summary>
        public int Offset { get { return m_offset; } }

        /// <summary>
        /// Copy Constructor
        /// </summary>
        /// <param name="oldCArray">Array to be copied</param>
        public CArray(CArray<T> oldCArray)
        {
            m_array = oldCArray.m_array;
            m_offset = oldCArray.m_offset;
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="x">the base array that will be pointed to</param>
        public CArray(T[] x)
        {
            m_array = x;
            m_offset = 0;
        }

        /// <summary>
        /// Reset the pointer to zero
        /// </summary>
        public void reset()
        {
            m_offset = 0;
        }

        /// <summary>
        /// Get and set array elements using the internal pointer offset
        /// </summary>
        /// <param name="loc">the array element (using the internal offset)</param>
        /// <returns>an element</returns>
        public T this[int loc]
        {
            get 
            {
                return m_array[loc + m_offset]; 
            }
            set 
            {
                m_array[loc + m_offset] = value; 
            }
        }

        /// <summary>
        /// Increment the pointer and return a new pointer
        /// </summary>
        /// <param name="x">Existing pointer class</param>
        /// <param name="i">Amount to increment the pointer by</param>
        /// <returns>a new pointer</returns>
        public static CArray<T> operator +(CArray<T> x, int i)
        {
            // y points to the same base array as x, with a zero offset
            CArray<T> y = new CArray<T>(x.m_array);
            // y has an offset equal to where x was plus i
            y.m_offset = x.m_offset + i;
            return y;
        }

        /// <summary>
        /// Decrement the pointer and return a new pointer
        /// </summary>
        /// <param name="x">Existing pointer class</param>
        /// <param name="i">Amount to increment the pointer by</param>
        /// <returns>a new pointer</returns>
        public static CArray<T> operator -(CArray<T> x, int i)
        {
            CArray<T> y = new CArray<T>(x.m_array);
            y.m_offset = x.m_offset - i;
            return y;
        }

        /// <summary>
        /// Increment the existing pointer by one.
        /// </summary>
        /// x++ is similar to x = x + 1.
        /// x++ and ++x are both covered.
        /// <param name="x">Existing pointer class</param>
        /// <returns>The incremented pointer</returns>
        public static CArray<T> operator ++(CArray<T> x)
        {
            CArray<T> y = x;
            y.m_offset += 1;
            return y;
        }

        /// <summary>
        /// Decrement the existing pointer by one.
        /// </summary>
        /// x-- is similar to x = x - 1.
        /// x-- and --x are both covered.
        /// <param name="x">Existing pointer class</param>
        /// <returns>The decremented pointer</returns>
        public static CArray<T> operator --(CArray<T> x)
        {
            CArray<T> y = x;
            y.m_offset -= 1;
            return y;
        }

        /// <summary>
        /// Create a new copy that points to the same memory, but with an independent offset
        /// </summary>
        /// Used when another CArray is needed to the same memory
        /// <returns>A separate copy</returns>
        public CArray<T> copy()
        {
            // point to the same memory location
            CArray<T> y = new CArray<T>(this.m_array);
            // copy the current offset value (copies the int)
            y.m_offset = this.m_offset;
            return y;
        }
    }
}
