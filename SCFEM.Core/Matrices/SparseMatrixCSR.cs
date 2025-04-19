using System;
using System.Collections.Generic;
using System.Linq;

namespace SCFEM.Core.Matrices
{
    public class SparseMatrixCSR
    {
        private readonly double[] _values;        // Non-zero values
        private readonly int[] _columnIndices;    // Column indices of non-zero values
        private readonly int[] _rowPointers;      // Starting index of each row in values array
        private readonly int _rows;
        private readonly int _cols;

        public SparseMatrixCSR(double[] values, int[] columnIndices, int[] rowPointers, int rows, int cols)
        {
            _values = values;
            _columnIndices = columnIndices;
            _rowPointers = rowPointers;
            _rows = rows;
            _cols = cols;
        }

        public int RowCount => _rows;
        public int ColumnCount => _cols;

        public double[] Multiply(double[] vector)
        {
            if (vector.Length != _cols)
                throw new ArgumentException("Vector length must match matrix columns");

            var result = new double[_rows];
            for (int i = 0; i < _rows; i++)
            {
                var rowStart = _rowPointers[i];
                var rowEnd = _rowPointers[i + 1];
                for (int j = rowStart; j < rowEnd; j++)
                {
                    result[i] += _values[j] * vector[_columnIndices[j]];
                }
            }
            return result;
        }

        public void Solve(double[] rhs, double[] solution)
        {
            if (rhs.Length != _rows || solution.Length != _rows)
                throw new ArgumentException("Vector lengths must match matrix dimensions");

            // Simple Gauss-Seidel iteration
            var maxIterations = 1000;
            var tolerance = 1e-10;
            Array.Copy(rhs, solution, rhs.Length);

            for (int iter = 0; iter < maxIterations; iter++)
            {
                var maxDiff = 0.0;
                for (int i = 0; i < _rows; i++)
                {
                    var sum = 0.0;
                    var diag = 0.0;
                    var rowStart = _rowPointers[i];
                    var rowEnd = _rowPointers[i + 1];

                    for (int j = rowStart; j < rowEnd; j++)
                    {
                        var col = _columnIndices[j];
                        if (col == i)
                        {
                            diag = _values[j];
                        }
                        else
                        {
                            sum += _values[j] * solution[col];
                        }
                    }

                    if (System.Math.Abs(diag) < 1e-10)
                        throw new InvalidOperationException("Zero diagonal element found");

                    var newValue = (rhs[i] - sum) / diag;
                    maxDiff = System.Math.Max(maxDiff, System.Math.Abs(newValue - solution[i]));
                    solution[i] = newValue;
                }

                if (maxDiff < tolerance)
                    break;
            }
        }

        public static SparseMatrixCSR FromSparseMatrix(SparseMatrix matrix)
        {
            var values = new List<double>();
            var columnIndices = new List<int>();
            var rowPointers = new List<int> { 0 };

            for (int i = 0; i < matrix.RowCount; i++)
            {
                var rowStart = values.Count;
                for (int j = 0; j < matrix.ColumnCount; j++)
                {
                    var value = matrix[i, j];
                    if (System.Math.Abs(value) > 1e-10)
                    {
                        values.Add(value);
                        columnIndices.Add(j);
                    }
                }
                rowPointers.Add(values.Count);
            }

            return new SparseMatrixCSR(
                values.ToArray(),
                columnIndices.ToArray(),
                rowPointers.ToArray(),
                matrix.RowCount,
                matrix.ColumnCount
            );
        }
    }
} 