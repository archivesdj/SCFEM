using System;
using System.Collections.Generic;
using System.Linq;

namespace SCFEM.Core.Matrices
{
    public class SparseMatrix
    {
        private readonly Dictionary<(int, int), double> _values;
        private readonly int _rows;
        private readonly int _columns;

        public int Rows => _rows;
        public int Columns => _columns;

        public SparseMatrix(int rows, int columns)
        {
            _rows = rows;
            _columns = columns;
            _values = new Dictionary<(int, int), double>();
        }

        public double this[int row, int col]
        {
            get => _values.TryGetValue((row, col), out double value) ? value : 0.0;
            set
            {
                if (value != 0.0)
                {
                    _values[(row, col)] = value;
                }
                else
                {
                    _values.Remove((row, col));
                }
            }
        }

        public int RowCount => _rows;
        public int ColumnCount => _columns;

        public void Clear()
        {
            _values.Clear();
        }

        public SparseMatrixCSR ToCSR()
        {
            return SparseMatrixCSR.FromSparseMatrix(this);
        }
    }
} 