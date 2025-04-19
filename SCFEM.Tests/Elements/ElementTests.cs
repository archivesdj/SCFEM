using System;
using System.Numerics;
using System.Collections.Generic;
using SCFEM.Core.Meshes;
using SCFEM.Core.Matrices;
using SCFEM.Core.Meshes.Elements;
using Xunit;

namespace SCFEM.Tests.Elements
{
    public class ElementTests
    {
        private readonly Mesh _mesh;
        private readonly MaterialProperties _materialProperties;

        public ElementTests()
        {
            _mesh = new Mesh();
            _materialProperties = new MaterialProperties { Conductivity = 1.0 };
        }

        [Fact]
        public void TetrahedronElement_StiffnessMatrix_IsSymmetric()
        {
            var nodes = new Node[]
            {
                new Node(0, new[] { 0.0, 0.0, 0.0 }),
                new Node(1, new[] { 1.0, 0.0, 0.0 }),
                new Node(2, new[] { 0.0, 1.0, 0.0 }),
                new Node(3, new[] { 0.0, 0.0, 1.0 })
            };

            var element = new TetrahedronElement(nodes, "test");
            var K = element.CalculateStiffnessMatrix(_materialProperties);
            Assert.True(IsSymmetric(K));
        }

        [Fact]
        public void TriangleElement_StiffnessMatrix_IsSymmetric()
        {
            var nodes = new Node[]
            {
                new Node(0, new[] { 0.0, 0.0, 0.0 }),
                new Node(1, new[] { 1.0, 0.0, 0.0 }),
                new Node(2, new[] { 0.0, 1.0, 0.0 })
            };

            var element = new TriangleElement(nodes, "test");
            var K = element.CalculateStiffnessMatrix(_materialProperties);
            Assert.True(IsSymmetric(K));
        }

        [Fact]
        public void LineElement_StiffnessMatrix_IsCorrect()
        {
            // Arrange
            var nodes = new Node[]
            {
                new Node(1, new[] { 0.0, 0.0, 0.0 }),
                new Node(2, new[] { 1.0, 0.0, 0.0 })
            };
            var element = new LineElement(nodes, "line");

            // Act
            var stiffnessMatrix = element.CalculateStiffnessMatrix(_materialProperties);

            // Assert
            Assert.Equal(2, stiffnessMatrix.GetLength(0));
            Assert.Equal(2, stiffnessMatrix.GetLength(1));
            Assert.Equal(1.0, stiffnessMatrix[0, 0], 6);
            Assert.Equal(-1.0, stiffnessMatrix[0, 1], 6);
            Assert.Equal(-1.0, stiffnessMatrix[1, 0], 6);
            Assert.Equal(1.0, stiffnessMatrix[1, 1], 6);
        }

        [Fact]
        public void TriangleElement_StiffnessMatrix_IsCorrect()
        {
            // Arrange
            var nodes = new Node[]
            {
                new Node(1, new[] { 0.0, 0.0, 0.0 }),
                new Node(2, new[] { 1.0, 0.0, 0.0 }),
                new Node(3, new[] { 0.0, 1.0, 0.0 })
            };
            var element = new TriangleElement(nodes, "triangle");

            // Act
            var stiffnessMatrix = element.CalculateStiffnessMatrix(_materialProperties);

            // Assert
            Assert.Equal(3, stiffnessMatrix.GetLength(0));
            Assert.Equal(3, stiffnessMatrix.GetLength(1));
            Assert.Equal(0.5, stiffnessMatrix[0, 0], 6);
            Assert.Equal(-0.25, stiffnessMatrix[0, 1], 6);
            Assert.Equal(-0.25, stiffnessMatrix[0, 2], 6);
            Assert.Equal(-0.25, stiffnessMatrix[1, 0], 6);
            Assert.Equal(0.5, stiffnessMatrix[1, 1], 6);
            Assert.Equal(-0.25, stiffnessMatrix[1, 2], 6);
            Assert.Equal(-0.25, stiffnessMatrix[2, 0], 6);
            Assert.Equal(-0.25, stiffnessMatrix[2, 1], 6);
            Assert.Equal(0.5, stiffnessMatrix[2, 2], 6);
        }

        [Fact]
        public void TetrahedronElement_StiffnessMatrix_IsCorrect()
        {
            // Arrange
            var nodes = new Node[]
            {
                new Node(1, new[] { 0.0, 0.0, 0.0 }),
                new Node(2, new[] { 1.0, 0.0, 0.0 }),
                new Node(3, new[] { 0.0, 1.0, 0.0 }),
                new Node(4, new[] { 0.0, 0.0, 1.0 })
            };
            var element = new TetrahedronElement(nodes, "tetrahedron");

            // Act
            var stiffnessMatrix = element.CalculateStiffnessMatrix(_materialProperties);

            // Assert
            Assert.Equal(4, stiffnessMatrix.GetLength(0));
            Assert.Equal(4, stiffnessMatrix.GetLength(1));
            Assert.Equal(1.0/6.0, stiffnessMatrix[0, 0], 6);
            Assert.Equal(-1.0/12.0, stiffnessMatrix[0, 1], 6);
            Assert.Equal(-1.0/12.0, stiffnessMatrix[0, 2], 6);
            Assert.Equal(-1.0/12.0, stiffnessMatrix[0, 3], 6);
            Assert.Equal(-1.0/12.0, stiffnessMatrix[1, 0], 6);
            Assert.Equal(1.0/6.0, stiffnessMatrix[1, 1], 6);
            Assert.Equal(-1.0/12.0, stiffnessMatrix[1, 2], 6);
            Assert.Equal(-1.0/12.0, stiffnessMatrix[1, 3], 6);
            Assert.Equal(-1.0/12.0, stiffnessMatrix[2, 0], 6);
            Assert.Equal(-1.0/12.0, stiffnessMatrix[2, 1], 6);
            Assert.Equal(1.0/6.0, stiffnessMatrix[2, 2], 6);
            Assert.Equal(-1.0/12.0, stiffnessMatrix[2, 3], 6);
            Assert.Equal(-1.0/12.0, stiffnessMatrix[3, 0], 6);
            Assert.Equal(-1.0/12.0, stiffnessMatrix[3, 1], 6);
            Assert.Equal(-1.0/12.0, stiffnessMatrix[3, 2], 6);
            Assert.Equal(1.0/6.0, stiffnessMatrix[3, 3], 6);
        }

        [Fact]
        public void PrismElement_StiffnessMatrix_IsCorrect()
        {
            // Arrange
            var nodes = new Node[]
            {
                new Node(1, new[] { 0.0, 0.0, 0.0 }),
                new Node(2, new[] { 1.0, 0.0, 0.0 }),
                new Node(3, new[] { 0.0, 1.0, 0.0 }),
                new Node(4, new[] { 0.0, 0.0, 1.0 }),
                new Node(5, new[] { 1.0, 0.0, 1.0 }),
                new Node(6, new[] { 0.0, 1.0, 1.0 })
            };
            var element = new PrismElement(nodes, "prism");

            // Act
            var stiffnessMatrix = element.CalculateStiffnessMatrix(_materialProperties);

            // Assert
            Assert.Equal(6, stiffnessMatrix.GetLength(0));
            Assert.Equal(6, stiffnessMatrix.GetLength(1));
            // Add specific assertions for prism element stiffness matrix values
        }

        [Fact]
        public void HexahedronElement_StiffnessMatrix_IsCorrect()
        {
            // Arrange
            var nodes = new Node[]
            {
                new Node(1, new[] { 0.0, 0.0, 0.0 }),
                new Node(2, new[] { 1.0, 0.0, 0.0 }),
                new Node(3, new[] { 1.0, 1.0, 0.0 }),
                new Node(4, new[] { 0.0, 1.0, 0.0 }),
                new Node(5, new[] { 0.0, 0.0, 1.0 }),
                new Node(6, new[] { 1.0, 0.0, 1.0 }),
                new Node(7, new[] { 1.0, 1.0, 1.0 }),
                new Node(8, new[] { 0.0, 1.0, 1.0 })
            };
            var element = new HexahedronElement(nodes, "hexahedron");

            // Act
            var stiffnessMatrix = element.CalculateStiffnessMatrix(_materialProperties);

            // Assert
            Assert.Equal(8, stiffnessMatrix.GetLength(0));
            Assert.Equal(8, stiffnessMatrix.GetLength(1));
            // Add specific assertions for hexahedron element stiffness matrix values
        }

        private bool IsSymmetric(double[,] matrix)
        {
            var n = matrix.GetLength(0);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (Math.Abs(matrix[i, j] - matrix[j, i]) > 1e-10)
                    {
                        return false;
                    }
                }
            }
            return true;
        }
    }
} 