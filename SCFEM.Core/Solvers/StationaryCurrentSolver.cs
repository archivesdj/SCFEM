using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using SCFEM.Core.Meshes;
using SCFEM.Core.Solver.BoundaryConditions;
using SCFEM.Core.Matrices;
using SCFEM.Core.Meshes.Elements;

namespace SCFEM.Core.Solvers
{
    public class StationaryCurrentSolver
    {
        private readonly Mesh _mesh;
        private SparseMatrix _globalStiffnessMatrix;
        private readonly double[] _forceVector;
        private readonly double[] _solution;
        private readonly Dictionary<string, MaterialProperties> _materialPropertiesByGroup;

        public StationaryCurrentSolver(Meshes.Mesh mesh)
        {
            _mesh = mesh;
            _materialPropertiesByGroup = new Dictionary<string, MaterialProperties>();

            var numNodes = _mesh.Nodes.Count;
            _globalStiffnessMatrix = new SparseMatrix(numNodes, numNodes);
             _forceVector = new double[numNodes];
            _solution = new double[numNodes];
        }

        public void SetMaterialProperties(string physicalGroup, MaterialProperties properties)
        {
            _materialPropertiesByGroup[physicalGroup] = properties;
        }

        public void Assemble()
        {
            foreach (var element in _mesh.Elements)
            {
                var materialProperties = _materialPropertiesByGroup[element.PhysicalGroup];
                var elementStiffness = element.CalculateStiffnessMatrix(materialProperties);

                for (int i = 0; i < element.Nodes.Length; i++)
                {
                    var nodeI = element.Nodes[i];
                    for (int j = 0; j < element.Nodes.Length; j++)
                    {
                        var nodeJ = element.Nodes[j];
                        _globalStiffnessMatrix[nodeI.Id, nodeJ.Id] += elementStiffness[i, j];
                    }
                }
            }

            foreach (var bc in _mesh.BoundaryConditions)
            {
                bc.Apply(_globalStiffnessMatrix, _forceVector);
            }
        }

        public void Solve()
        {
            // TODO: Implement solver
        }

        public double[] GetSolution()
        {
            return _solution;
        }

        public double[] GetCurrentDensity(Element element, MaterialProperties materialProperties)
        {
            // TODO: Implement current density calculation
            return new double[3];
        }

        public void SaveSolution(string filePath)
        {
            // TODO: Implement solution saving
        }
    }
} 