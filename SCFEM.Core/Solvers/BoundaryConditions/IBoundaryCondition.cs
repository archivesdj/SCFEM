using SCFEM.Core.Matrices;

namespace SCFEM.Core.Solver.BoundaryConditions
{
    public interface IBoundaryCondition
    {
        string PhysicalGroup { get; }
        void Apply(SparseMatrix globalStiffnessMatrix, double[] globalForceVector);
    }
} 