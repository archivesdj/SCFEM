using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace SCFEM.Core.Solver.BoundaryConditions
{
    public abstract class BoundaryCondition
    {
        public string PhysicalGroup { get; }
        public List<int> NodeIndices { get; set; }

        protected BoundaryCondition(string physicalGroup)
        {
            PhysicalGroup = physicalGroup;
            NodeIndices = new List<int>();
        }

        public abstract void Apply(MathNet.Numerics.LinearAlgebra.Double.SparseMatrix stiffnessMatrix, MathNet.Numerics.LinearAlgebra.Vector<double> rightHandSide);
    }
} 