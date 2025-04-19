using System;
using System.Collections.Generic;
using SCFEM.Core.Matrices;

namespace SCFEM.Core.Solver.BoundaryConditions
{
    public abstract class BoundaryCondition : IBoundaryCondition
    {
        public string PhysicalGroup { get; }
        public double Value { get; set; }
        public BoundaryConditionType Type { get; }
        public List<int> NodeIndices { get; set; }

        protected BoundaryCondition(string physicalGroup, double value, BoundaryConditionType type)
        {
            PhysicalGroup = physicalGroup;
            Value = value;
            Type = type;
            NodeIndices = new List<int>();
        }

        public abstract void Apply(SparseMatrix globalStiffnessMatrix, double[] globalForceVector);
    }

    public enum BoundaryConditionType
    {
        Dirichlet,
        Neumann
    }
} 