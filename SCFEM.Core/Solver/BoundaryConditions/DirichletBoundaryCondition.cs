using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace SCFEM.Core.Solver.BoundaryConditions
{
    public class DirichletBoundaryCondition : BoundaryCondition
    {
        public double Value { get; set; }

        public DirichletBoundaryCondition(string physicalGroup, double value) : base(physicalGroup)
        {
            Value = value;
        }

        public override void Apply(MathNet.Numerics.LinearAlgebra.Double.SparseMatrix stiffnessMatrix, MathNet.Numerics.LinearAlgebra.Vector<double> rightHandSide)
        {
            foreach (var nodeIndex in NodeIndices)
            {
                // Set diagonal element to 1 and other elements in the row to 0
                for (int j = 0; j < stiffnessMatrix.ColumnCount; j++)
                {
                    if (j == nodeIndex)
                    {
                        stiffnessMatrix[nodeIndex, j] = 1.0;
                    }
                    else
                    {
                        stiffnessMatrix[nodeIndex, j] = 0.0;
                    }
                }

                // Set right-hand side to the prescribed value
                rightHandSide[nodeIndex] = Value;
            }
        }
    }
} 