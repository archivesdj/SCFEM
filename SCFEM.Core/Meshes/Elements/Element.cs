using System;
using System.Collections.Generic;
using SCFEM.Core.Matrices;
using SCFEM.Core.Meshes;

namespace SCFEM.Core.Meshes.Elements
{
    public abstract class Element
    {
        public ElementType Type { get; }
        public Node[] Nodes { get; }
        public string PhysicalGroup { get; }

        protected Element(ElementType type, Node[] nodes, string physicalGroup)
        {
            Type = type;
            Nodes = nodes;
            PhysicalGroup = physicalGroup;
        }

        public abstract double[,] CalculateStiffnessMatrix(MaterialProperties materialProperties);
        public abstract double[] CalculateShapeFunctions(double xi, double eta, double zeta);
        public abstract double[] CalculateShapeFunctionGradients(double xi, double eta, double zeta);
        public abstract double CalculateJacobian(double xi, double eta, double zeta);
        public abstract double[] CalculateNaturalCoordinates(double[] point);

        protected static double[,] Multiply(double[,] a, double[,] b)
        {
            var rows = a.GetLength(0);
            var cols = b.GetLength(1);
            var result = new double[rows, cols];

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    for (int k = 0; k < a.GetLength(1); k++)
                    {
                        result[i, j] += a[i, k] * b[k, j];
                    }
                }
            }

            return result;
        }

        protected static double[,] Multiply(double[,] a, double b)
        {
            var rows = a.GetLength(0);
            var cols = a.GetLength(1);
            var result = new double[rows, cols];

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    for (int k = 0; k < a.GetLength(1); k++)
                    {
                        result[i, j] += a[i, k] * b;
                    }
                }
            }

            return result;
        }

        protected static double[,] Transpose(double[,] matrix)
        {
            var rows = matrix.GetLength(0);
            var cols = matrix.GetLength(1);
            var result = new double[cols, rows];

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    result[j, i] = matrix[i, j];
                }
            }

            return result;
        }

        protected static double[,] Add(double[,] a, double[,] b)
        {
            var rows = a.GetLength(0);
            var cols = a.GetLength(1);
            var result = new double[rows, cols];

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    result[i, j] = a[i, j] + b[i, j];
                }
            }

            return result;
        }

        protected static double Determinant(double[,] matrix)
        {
            if (matrix.GetLength(0) == 1)
            {
                return matrix[0, 0];
            }
            else if (matrix.GetLength(0) == 2)
            {
                return matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0];
            }
            else if (matrix.GetLength(0) == 3)
            {
                return matrix[0, 0] * (matrix[1, 1] * matrix[2, 2] - matrix[1, 2] * matrix[2, 1]) -
                       matrix[0, 1] * (matrix[1, 0] * matrix[2, 2] - matrix[1, 2] * matrix[2, 0]) +
                       matrix[0, 2] * (matrix[1, 0] * matrix[2, 1] - matrix[1, 1] * matrix[2, 0]);
            }
            else
            {
                throw new NotImplementedException("Determinant calculation for matrices larger than 3x3 is not implemented");
            }
        }

        protected static double[,] Inverse(double[,] matrix)
        {
            var det = Determinant(matrix);
            if (Math.Abs(det) < 1e-10)
            {
                throw new Exception("Matrix is singular");
            }

            if (matrix.GetLength(0) == 1)
            {
                return new double[,] { { 1.0 / matrix[0, 0] } };
            }
            else if (matrix.GetLength(0) == 2)
            {
                var invDet = 1.0 / det;
                return new double[,]
                {
                    { matrix[1, 1] * invDet, -matrix[0, 1] * invDet },
                    { -matrix[1, 0] * invDet, matrix[0, 0] * invDet }
                };
            }
            else if (matrix.GetLength(0) == 3)
            {
                var invDet = 1.0 / det;
                return new double[,]
                {
                    {
                        (matrix[1, 1] * matrix[2, 2] - matrix[1, 2] * matrix[2, 1]) * invDet,
                        (matrix[0, 2] * matrix[2, 1] - matrix[0, 1] * matrix[2, 2]) * invDet,
                        (matrix[0, 1] * matrix[1, 2] - matrix[0, 2] * matrix[1, 1]) * invDet
                    },
                    {
                        (matrix[1, 2] * matrix[2, 0] - matrix[1, 0] * matrix[2, 2]) * invDet,
                        (matrix[0, 0] * matrix[2, 2] - matrix[0, 2] * matrix[2, 0]) * invDet,
                        (matrix[0, 2] * matrix[1, 0] - matrix[0, 0] * matrix[1, 2]) * invDet
                    },
                    {
                        (matrix[1, 0] * matrix[2, 1] - matrix[1, 1] * matrix[2, 0]) * invDet,
                        (matrix[0, 1] * matrix[2, 0] - matrix[0, 0] * matrix[2, 1]) * invDet,
                        (matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0]) * invDet
                    }
                };
            }
            else
            {
                throw new NotImplementedException("Inverse calculation for matrices larger than 3x3 is not implemented");
            }
        }
    }
} 