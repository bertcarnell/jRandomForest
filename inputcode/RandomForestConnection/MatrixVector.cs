/*
 * Copyright 2013 Rob Carnell
 * 
 */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Diagnostics;

namespace RandomForestConnection
{
    /// <summary>
    /// Creates a flexible matrix which is stored as vector for passing to R
    /// </summary>
    /// <typeparam name="T"></typeparam>
    [Serializable]
    public class MatrixVector<T>
    {
        // values are organized by column in the matrix
        // 1 2 3 4 5 6 7 8 9 10 11 12
        // 1 5 9
        // 2 6 10
        // 3 7 11
        // 4 8 12

        /// <summary>
        /// Single dimensional array storage
        /// </summary>
        public T[] values;
        /// <summary>
        /// row and column counts
        /// </summary>
        public int nrow, ncol, dim3;
        /// <summary>
        /// Column Names
        /// </summary>
        public string[] columnNames;
        private int m_nDims;

        private int accessElement(int rowIndex, int colIndex)
        {
            return colIndex * nrow + rowIndex;
        }

        private int accessElement3d(int rowIndex, int colIndex, int dim3Index)
        {
            return dim3Index * nrow * ncol + colIndex * nrow + rowIndex;
        }

        /// <summary>
        /// Allocates a MatrixVector
        /// </summary>
        /// <param name="nrows">Number of Rows</param>
        /// <param name="ncols">Number of Columns</param>
        public MatrixVector(int nrows, int ncols)
        {
            if (nrows < 1 || ncols < 1)
                throw new ArgumentException(string.Format("Array dimension less than one:  dim1={0}, dim2={1}", nrows, ncols)); 
            values = new T[nrows * ncols];
            nrow = nrows;
            ncol = ncols;
            m_nDims = 2;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="dim1"></param>
        /// <param name="dim2"></param>
        /// <param name="dim3"></param>
        public MatrixVector(int dim1, int dim2, int dim3)
        {
            if (dim1 < 1 || dim2 < 1 || dim3 < 1)
                throw new ArgumentException(string.Format("Array dimension less than one:  dim1={0}, dim2={1}, dim3={2}", dim1, dim2, dim3)); 
            m_nDims = 3;
            nrow = dim1;
            ncol = dim2;
            this.dim3 = dim3;
            int product = dim1 * dim2 * dim3;
            values = new T[product];
        }
        /// <summary>
        /// Creates a MatrixVector from an existing two dimensional array
        /// </summary>
        /// <param name="matrix">A two dimensional array</param>
        public MatrixVector(T[,] matrix)
        {
            nrow = matrix.GetLength(0);
            ncol = matrix.GetLength(1);
            values = new T[nrow * ncol];
            for (int j = 0; j < ncol; j++)
            {
                for (int i = 0; i < nrow; i++)
                {
                    values[accessElement(i,j)] = matrix[i, j];
                }
            }
            m_nDims = 2;
        }
        /// <summary>
        /// Creates a MatrixVector from an existing three dimensional array
        /// </summary>
        /// <param name="matrix">A two dimensional array</param>
        public MatrixVector(T[,,] matrix)
        {
            nrow = matrix.GetLength(0);
            ncol = matrix.GetLength(1);
            dim3 = matrix.GetLength(2);
            int product = nrow * ncol * dim3;
            values = new T[product];
            for (int k = 0; k < dim3; k++)
            {
                for (int j = 0; j < ncol; j++)
                {
                    for (int i = 0; i < nrow; i++)
                    {
                        values[accessElement3d(i, j, k)] = matrix[i, j, k];
                    }
                }
            }
            m_nDims = 3;
        }
        /// <summary>
        /// Construct a matrix vector from an array
        /// </summary>
        /// <param name="array">array of type</param>
        /// <param name="nrows">rows for the matrix vector</param>
        /// <param name="ncols">columns for the matrix vector</param>
        public MatrixVector(T[] array, int nrows, int ncols)
        {
            if (array.Length != nrows * ncols)
                throw new Exception(string.Format("Can not create MatrixVector with incompatible array (length = {0}) and rows/columns ({1} and {2}", array.Length, nrows, ncols));
            nrow = nrows;
            ncol = ncols;
            values = array;
            m_nDims = 2;
        }
        /// <summary>
        /// Construct a matrix vector from an array
        /// </summary>
        /// <param name="array">array of type</param>
        /// <param name="nrows">rows for the matrix vector</param>
        /// <param name="ncols">columns for the matrix vector</param>
        /// <param name="dim3">3rd dimension of the matrix</param>
        public MatrixVector(T[] array, int nrows, int ncols, int dim3)
        {
            if (array.Length != nrows * ncols * dim3)
                throw new Exception(string.Format("Can not create MatrixVector with incompatible array (length = {0}), rows={1}, columbs={2}, and dim3={3}", array.Length, nrows, ncols, dim3));
            nrow = nrows;
            ncol = ncols;
            this.dim3 = dim3;
            m_nDims = 3;
            values = array;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="nrows"></param>
        /// <param name="ncols"></param>
        /// <param name="fillValue"></param>
        public MatrixVector(int nrows, int ncols, T fillValue)
        {
            values = new T[nrows * ncols];
            nrow = nrows;
            ncol = ncols;
            for (int i = 0; i < values.Length; i++)
                values[i] = fillValue;
            m_nDims = 2;
        }
        /// <summary>
        /// Returns a value at a matrix location
        /// </summary>
        /// <param name="rowLoc">matrix row for the value</param>
        /// <param name="colLoc">matrix column for the value</param>
        /// <returns>value</returns>
        public T this[int rowLoc, int colLoc]
        {
            get { return (values[accessElement(rowLoc, colLoc)]); }
            set { values[accessElement(rowLoc, colLoc)] = value; }
        }
        /// <summary>
        /// Returns a value at a matrix location
        /// </summary>
        /// <param name="rowLoc">matrix row for the value</param>
        /// <param name="colLoc">matrix column for the value</param>
        /// <param name="dim3Loc">matrix 3rd dimension for the value</param>
        /// <returns>value</returns>
        public T this[int rowLoc, int colLoc, int dim3Loc]
        {
            get { return (values[accessElement3d(rowLoc, colLoc, dim3Loc)]); }
            set { values[accessElement3d(rowLoc, colLoc, dim3Loc)] = value; }
        }
        /// <summary>
        /// Returns a value at an array location
        /// </summary>
        /// <param name="loc">array index</param>
        /// <returns>value</returns>
        public T this[int loc]
        {
            get { return (values[loc]); }
            set { values[loc] = value; }
        }
        /// <summary>
        /// Gets the matrix representation of the MatrixVector
        /// </summary>
        /// <returns>Two dimensional array</returns>
        public T[,] getMatrix()
        {
            if (m_nDims != 2)
                throw new InvalidCastException(String.Format("A two dimensional matrix cannot be returned from a MatrixVector of dimension {0}", m_nDims));
            T[,] result = new T[nrow, ncol];
            for (int j = 0; j < ncol; j++)
            {
                for (int i = 0; i < nrow; i++)
                {
                    result[i, j] = values[accessElement(i,j)];
                }
            }
            return result;
        }
        /// <summary>
        /// Gets the matrix representation of the MatrixVector
        /// </summary>
        /// <returns>Three dimensional array</returns>
        public T[,,] get3dMatrix()
        {
            if (m_nDims != 3)
                throw new InvalidCastException(String.Format("A three dimensional matrix cannot be returned from a MatrixVector of dimension {0}", m_nDims));
            T[,,] result = new T[nrow, ncol, dim3];
            for (int k = 0; k < dim3; k++)
            {
                for (int j = 0; j < ncol; j++)
                {
                    for (int i = 0; i < nrow; i++)
                    {
                        result[i, j, k] = values[accessElement3d(i, j, k)];
                    }
                }
            }
            return result;
        }
        /// <summary>
        /// Gets the vector implementation of the MatrixVector
        /// </summary>
        /// <returns>One dimensional array</returns>
        public T[] getVector()
        {
            return values;
        }
        /// <summary>
        /// Fills the matrix vector with a value to initialize it
        /// </summary>
        /// <param name="x">fill parameter</param>
        public void fill(T x)
        {
            for (int i = 0; i < values.Length; i++)
                values[i] = x;
        }
        /// <summary>
        /// Transpose the MatrixVector
        /// </summary>
        public void transpose()
        {
            if (m_nDims != 2)
                throw new InvalidCastException("Only a two dimensional matrix can be transposed");
            T[] temp = new T[values.Length];
            int newRow = ncol;
            int newCol = nrow;
            for (int j = 0; j < ncol; j++)
            {
                for (int i = 0; i < nrow; i++)
                {
                    temp[i * ncol + j] = values[accessElement(i, j)];
                }
            }
            values = temp;
            nrow = newRow;
            ncol = newCol;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="dimensionOrder"></param>
        public void permute(int[] dimensionOrder)
        {
            if (m_nDims != 3)
                throw new NotImplementedException("permute is only implemented for 3 dimensional arrays");
            if (dimensionOrder.Min() != 0 || dimensionOrder.Max() != 2 || !dimensionOrder.Contains(1))
                throw new ArgumentException("dimensionOrder must be 3 zero-based indexes", "dimensionOrder");
            T[] temp = new T[values.Length];
            int[] tempLengths = new int[3];
            tempLengths[0] = nrow;
            tempLengths[1] = ncol;
            tempLengths[2] = dim3;
            int newRow = tempLengths[dimensionOrder[0]];
            int newCol = tempLengths[dimensionOrder[1]];
            int newDim3 = tempLengths[dimensionOrder[2]];
            int[] tempIndexes = new int[3];

            for (int k = 0; k < dim3; k++)
            {
                for (int j = 0; j < ncol; j++)
                {
                    for (int i = 0; i < nrow; i++)
                    {
                        tempIndexes[0] = i;
                        tempIndexes[1] = j;
                        tempIndexes[2] = k;
                        temp[tempIndexes[dimensionOrder[2]] * newRow * newCol + 
                             tempIndexes[dimensionOrder[1]] * newRow + 
                             tempIndexes[dimensionOrder[0]]] = values[accessElement3d(i, j, k)];
                    }
                }
            }
            values = temp;
            nrow = newRow;
            ncol = newCol;
            dim3 = newDim3;
        }

        /// <summary>
        /// Calculate the row sums of the Matrix vector
        /// </summary>
        /// <returns>An array of row sums</returns>
        public double[] rowSum()
        {
            // need to cast everything to a type that can do the summing...  T+T will not work in .NET3.5
            double[] temp = Array.ConvertAll<T, double>(values, new Converter<T,double>(a => Convert.ToDouble(a))).ToArray();
            double[] result = new double[nrow];
            for (int i = 0; i < nrow; i++)
            {
                for (int j = 0; j < ncol; j++)
                {
                    result[i] += temp[accessElement(i,j)];
                }
            }
            return result;
        }

        /// <summary>
        /// normalize the maxtric vector across the rows such that the row sums are all equal to 1
        /// </summary>
        public void normalizeByRow()
        {
            // can only be performed on MatrixVectors of type double and float
            if (typeof(T) != typeof(double) && typeof(T) != typeof(float))
                throw new Exception("Normalize By Row cannot be performed on a MatrixVector that is not of type double or float");
            double[] rowSum = this.rowSum().Cast<double>().ToArray();
            double[] temp = values.Cast<double>().ToArray();
            Debug.Assert(rowSum.Length == nrow);

            for (int i = 0; i < nrow; i++) // rows of out.votes
            {
                for (int j = 0; j < ncol; j++)
                {
                    if (rowSum[i] != 0.0)
                        temp[accessElement(i, j)] /= rowSum[i];
                    else
                        temp[accessElement(i, j)] = 0.0;
                }
            }
            values = temp.Cast<T>().ToArray();
        }

        /// <summary>
        /// Returns a column of a MatrixVector
        /// </summary>
        /// <param name="colIndex">Column Index</param>
        /// <returns>An array for the column</returns>
        public T[] GetColumn(int colIndex)
        {
            if (colIndex >= ncol)
                throw new ArgumentException("Column index is outside the Matrix bounds", "colIndex");
            T[] result = new T[nrow];
            for (int i = 0; i < nrow; i++)
                result[i] = this[i, colIndex];
            return result;
        }
    }
}
