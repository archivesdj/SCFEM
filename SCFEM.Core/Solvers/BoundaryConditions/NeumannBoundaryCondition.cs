using System;
using SCFEM.Core.Matrices;

namespace SCFEM.Core.Solver.BoundaryConditions
{
    public class NeumannBoundaryCondition : BoundaryCondition
    {
        public NeumannBoundaryCondition(string physicalGroup, double value) 
            : base(physicalGroup, value, BoundaryConditionType.Neumann)
        {
        }

        public override void Apply(SparseMatrix globalStiffnessMatrix, double[] globalForceVector)
        {
            // For Neumann boundary conditions, we only modify the force vector
            foreach (var nodeIndex in NodeIndices)
            {
                globalForceVector[nodeIndex] += Value;
            }
        }
    }
} 