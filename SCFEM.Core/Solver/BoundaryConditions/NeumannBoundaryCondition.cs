using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace SCFEM.Core.Solver.BoundaryConditions
{
    public class NeumannBoundaryCondition : BoundaryCondition
    {
        public double Value { get; set; }

        public NeumannBoundaryCondition(string physicalGroup, double value) : base(physicalGroup)
        {
            Value = value;
        }

        public override void Apply(MathNet.Numerics.LinearAlgebra.Double.SparseMatrix stiffnessMatrix, MathNet.Numerics.LinearAlgebra.Vector<double> rightHandSide)
        {
            foreach (var nodeIndex in NodeIndices)
            {
                // Add the Neumann boundary condition value to the right-hand side
                rightHandSide[nodeIndex] += Value;
            }
        }
    }
} 