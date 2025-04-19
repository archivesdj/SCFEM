using System;
using SCFEM.Core.Matrices;

namespace SCFEM.Core.Solver.BoundaryConditions
{
    public class DirichletBoundaryCondition : BoundaryCondition
    {
        public DirichletBoundaryCondition(string physicalGroup, double value) 
            : base(physicalGroup, value, BoundaryConditionType.Dirichlet)
        {
        }

        public override void Apply(SparseMatrix globalStiffnessMatrix, double[] globalForceVector)
        {
            // For Dirichlet boundary conditions, we modify both the stiffness matrix and force vector
            foreach (var nodeIndex in NodeIndices)
            {
                // Set the row to zero except for the diagonal
                for (int j = 0; j < globalStiffnessMatrix.Columns; j++)
                {
                    if (j != nodeIndex)
                    {
                        globalStiffnessMatrix[nodeIndex, j] = 0.0;
                    }
                }

                // Set the diagonal to 1
                globalStiffnessMatrix[nodeIndex, nodeIndex] = 1.0;

                // Set the force vector entry to the prescribed value
                globalForceVector[nodeIndex] = Value;
            }
        }
    }
} 