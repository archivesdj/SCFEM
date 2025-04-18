using System;
using System.Collections.Generic;
using System.IO;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using SCFEM.Core.Mesh.Elements;
using SCFEM.Core.Solver.BoundaryConditions;

namespace SCFEM.Core.Solver
{
    public class StationaryCurrentSolver
    {
        private readonly SCFEM.Core.Mesh.Mesh _mesh;
        private SparseMatrix _stiffnessMatrix;
        private MathNet.Numerics.LinearAlgebra.Vector<double> _rightHandSide;
        private MathNet.Numerics.LinearAlgebra.Vector<double> _solution;
        private List<BoundaryCondition> _boundaryConditions;

        public StationaryCurrentSolver(SCFEM.Core.Mesh.Mesh mesh)
        {
            _mesh = mesh;
            _boundaryConditions = new List<BoundaryCondition>();
        }

        public void AddBoundaryCondition(BoundaryCondition boundaryCondition)
        {
            // Find nodes belonging to the physical group
            foreach (var element in _mesh.Elements)
            {
                if (element.PhysicalGroup == boundaryCondition.PhysicalGroup)
                {
                    foreach (var nodeIndex in element.NodeIndices)
                    {
                        if (!boundaryCondition.NodeIndices.Contains(nodeIndex))
                        {
                            boundaryCondition.NodeIndices.Add(nodeIndex);
                        }
                    }
                }
            }
            _boundaryConditions.Add(boundaryCondition);
        }

        public void AssembleSystem()
        {
            int numNodes = _mesh.Nodes.Count;
            _stiffnessMatrix = SparseMatrix.Create(numNodes, numNodes, 0.0);
            _rightHandSide = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(numNodes, 0.0);

            // Assemble element matrices
            foreach (var element in _mesh.Elements)
            {
                // Get element nodes
                var elementNodes = new Vector3[element.NodeIndices.Length];
                for (int i = 0; i < element.NodeIndices.Length; i++)
                {
                    elementNodes[i] = _mesh.Nodes[element.NodeIndices[i]];
                }

                // Get material property (conductivity)
                double conductivity = 1.0; // Default value
                if (!string.IsNullOrEmpty(element.PhysicalGroup) && 
                    _mesh.MaterialProperties.ContainsKey(element.PhysicalGroup))
                {
                    conductivity = _mesh.MaterialProperties[element.PhysicalGroup];
                }

                // Calculate element stiffness matrix
                var elementStiffness = element.CalculateStiffnessMatrix(elementNodes, conductivity);

                // Assemble into global matrix
                for (int i = 0; i < element.NodeIndices.Length; i++)
                {
                    int globalI = element.NodeIndices[i];
                    for (int j = 0; j < element.NodeIndices.Length; j++)
                    {
                        int globalJ = element.NodeIndices[j];
                        _stiffnessMatrix[globalI, globalJ] += elementStiffness[i, j];
                    }
                }
            }
        }

        public void ApplyBoundaryConditions()
        {
            foreach (var boundaryCondition in _boundaryConditions)
            {
                boundaryCondition.Apply(_stiffnessMatrix, _rightHandSide);
            }
        }

        public void Solve()
        {
            // Create a solver using LU decomposition
            var solver = _stiffnessMatrix.LU();

            // Solve the system
            _solution = solver.Solve(_rightHandSide);
        }

        public void ExportToVTK(string filePath)
        {
            using (var writer = new StreamWriter(filePath))
            {
                // Write VTK header
                writer.WriteLine("# vtk DataFile Version 2.0");
                writer.WriteLine("SCFEM Solution");
                writer.WriteLine("ASCII");
                writer.WriteLine("DATASET UNSTRUCTURED_GRID");

                // Write points
                writer.WriteLine($"POINTS {_mesh.Nodes.Count} double");
                foreach (var node in _mesh.Nodes)
                {
                    writer.WriteLine($"{node.X} {node.Y} {node.Z}");
                }

                // Write cells
                int totalCellDataSize = 0;
                foreach (var element in _mesh.Elements)
                {
                    totalCellDataSize += 1 + element.NodeIndices.Length; // 1 for cell type + number of nodes
                }

                writer.WriteLine($"CELLS {_mesh.Elements.Count} {totalCellDataSize}");
                foreach (var element in _mesh.Elements)
                {
                    writer.Write($"{element.NodeIndices.Length}");
                    foreach (var nodeIndex in element.NodeIndices)
                    {
                        writer.Write($" {nodeIndex}");
                    }
                    writer.WriteLine();
                }

                // Write cell types
                writer.WriteLine($"CELL_TYPES {_mesh.Elements.Count}");
                foreach (var element in _mesh.Elements)
                {
                    int vtkCellType = element.Type switch
                    {
                        ElementType.Tetrahedron => 10, // VTK_TETRA
                        ElementType.Hexahedron => 12,  // VTK_HEXAHEDRON
                        ElementType.Prism => 13,       // VTK_WEDGE
                        _ => throw new NotSupportedException($"Unsupported element type: {element.Type}")
                    };
                    writer.WriteLine(vtkCellType);
                }

                // Write point data (solution)
                writer.WriteLine($"POINT_DATA {_mesh.Nodes.Count}");
                writer.WriteLine("SCALARS potential double 1");
                writer.WriteLine("LOOKUP_TABLE default");
                foreach (var value in _solution)
                {
                    writer.WriteLine(value);
                }

                // Write cell data (physical groups)
                writer.WriteLine($"CELL_DATA {_mesh.Elements.Count}");
                writer.WriteLine("SCALARS physical_group int 1");
                writer.WriteLine("LOOKUP_TABLE default");
                foreach (var element in _mesh.Elements)
                {
                    writer.WriteLine(element.PhysicalGroup ?? "0");
                }
            }
        }
    }
} 