using System;
using System.Linq;
using Xunit;
using SCFEM.Core.Meshes;
using SCFEM.Core.Meshes.Elements;
using SCFEM.Core.Meshes.Nodes;
using SCFEM.Core.Materials;
using SCFEM.Core.Matrix;

namespace SCFEM.Core.Tests.Mesh.Elements
{
    public class ElementTests
    {
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
            var materialProperties = new MaterialProperties { Conductivity = 1.0 };

            // Act
            var stiffnessMatrix = element.CalculateStiffnessMatrix(materialProperties);

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
            var materialProperties = new MaterialProperties { Conductivity = 1.0 };

            // Act
            var stiffnessMatrix = element.CalculateStiffnessMatrix(materialProperties);

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
            var materialProperties = new MaterialProperties { Conductivity = 1.0 };

            // Act
            var stiffnessMatrix = element.CalculateStiffnessMatrix(materialProperties);

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
            var materialProperties = new MaterialProperties { Conductivity = 1.0 };

            // Act
            var stiffnessMatrix = element.CalculateStiffnessMatrix(materialProperties);

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
            var materialProperties = new MaterialProperties { Conductivity = 1.0 };

            // Act
            var stiffnessMatrix = element.CalculateStiffnessMatrix(materialProperties);

            // Assert
            Assert.Equal(8, stiffnessMatrix.GetLength(0));
            Assert.Equal(8, stiffnessMatrix.GetLength(1));
            // Add specific assertions for hexahedron element stiffness matrix values
        }

        [Fact]
        public void LineElement_ShapeFunctions_AreCorrect()
        {
            // Arrange
            var nodes = new Node[]
            {
                new Node(1, new[] { 0.0, 0.0, 0.0 }),
                new Node(2, new[] { 1.0, 0.0, 0.0 })
            };
            var element = new LineElement(nodes, "line");

            // Act
            var shapeFunctions = element.CalculateShapeFunctions(0.5, 0, 0);

            // Assert
            Assert.Equal(2, shapeFunctions.Length);
            Assert.Equal(0.5, shapeFunctions[0], 6);
            Assert.Equal(0.5, shapeFunctions[1], 6);
        }

        [Fact]
        public void TriangleElement_ShapeFunctions_AreCorrect()
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
            var shapeFunctions = element.CalculateShapeFunctions(1.0/3.0, 1.0/3.0, 0);

            // Assert
            Assert.Equal(3, shapeFunctions.Length);
            Assert.Equal(1.0/3.0, shapeFunctions[0], 6);
            Assert.Equal(1.0/3.0, shapeFunctions[1], 6);
            Assert.Equal(1.0/3.0, shapeFunctions[2], 6);
        }

        [Fact]
        public void TetrahedronElement_ShapeFunctions_AreCorrect()
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
            var shapeFunctions = element.CalculateShapeFunctions(0.25, 0.25, 0.25);

            // Assert
            Assert.Equal(4, shapeFunctions.Length);
            Assert.Equal(0.25, shapeFunctions[0], 6);
            Assert.Equal(0.25, shapeFunctions[1], 6);
            Assert.Equal(0.25, shapeFunctions[2], 6);
            Assert.Equal(0.25, shapeFunctions[3], 6);
        }

        [Fact]
        public void PrismElement_ShapeFunctions_AreCorrect()
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
            var shapeFunctions = element.CalculateShapeFunctions(1.0/3.0, 1.0/3.0, 0.5);

            // Assert
            Assert.Equal(6, shapeFunctions.Length);
            Assert.Equal(1.0/6.0, shapeFunctions[0], 6);
            Assert.Equal(1.0/6.0, shapeFunctions[1], 6);
            Assert.Equal(1.0/6.0, shapeFunctions[2], 6);
            Assert.Equal(1.0/6.0, shapeFunctions[3], 6);
            Assert.Equal(1.0/6.0, shapeFunctions[4], 6);
            Assert.Equal(1.0/6.0, shapeFunctions[5], 6);
        }

        [Fact]
        public void HexahedronElement_ShapeFunctions_AreCorrect()
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
            var shapeFunctions = element.CalculateShapeFunctions(0, 0, 0);

            // Assert
            Assert.Equal(8, shapeFunctions.Length);
            Assert.Equal(1.0, shapeFunctions[0], 6);
            Assert.Equal(0.0, shapeFunctions[1], 6);
            Assert.Equal(0.0, shapeFunctions[2], 6);
            Assert.Equal(0.0, shapeFunctions[3], 6);
            Assert.Equal(0.0, shapeFunctions[4], 6);
            Assert.Equal(0.0, shapeFunctions[5], 6);
            Assert.Equal(0.0, shapeFunctions[6], 6);
            Assert.Equal(0.0, shapeFunctions[7], 6);
        }

        [Fact]
        public void Element_CalculateShapeFunctions_ReturnsCorrectValues()
        {
            // Arrange
            var nodes = new[]
            {
                new Node(0, new[] { 0.0, 0.0, 0.0 }),
                new Node(1, new[] { 1.0, 0.0, 0.0 }),
                new Node(2, new[] { 0.0, 1.0, 0.0 })
            };
            var element = new TriangleElement(nodes, "test");

            // Act
            var shapeFunctions = element.CalculateShapeFunctions(0.5, 0.5, 0.0);

            // Assert
            Assert.Equal(3, shapeFunctions.Length);
            Assert.Equal(0.0, shapeFunctions[0], 6);
            Assert.Equal(0.5, shapeFunctions[1], 6);
            Assert.Equal(0.5, shapeFunctions[2], 6);
        }

        [Fact]
        public void Element_CalculateShapeFunctionGradients_ReturnsCorrectValues()
        {
            // Arrange
            var nodes = new[]
            {
                new Node(0, new[] { 0.0, 0.0, 0.0 }),
                new Node(1, new[] { 1.0, 0.0, 0.0 }),
                new Node(2, new[] { 0.0, 1.0, 0.0 })
            };
            var element = new TriangleElement(nodes, "test");

            // Act
            var gradients = element.CalculateShapeFunctionGradients(0.5, 0.5, 0.0);

            // Assert
            Assert.Equal(6, gradients.Length);
            Assert.Equal(-1.0, gradients[0], 6);
            Assert.Equal(-1.0, gradients[1], 6);
            Assert.Equal(1.0, gradients[2], 6);
            Assert.Equal(0.0, gradients[3], 6);
            Assert.Equal(0.0, gradients[4], 6);
            Assert.Equal(1.0, gradients[5], 6);
        }

        [Fact]
        public void Element_CalculateJacobian_ReturnsCorrectValue()
        {
            // Arrange
            var nodes = new[]
            {
                new Node(0, new[] { 0.0, 0.0, 0.0 }),
                new Node(1, new[] { 1.0, 0.0, 0.0 }),
                new Node(2, new[] { 0.0, 1.0, 0.0 })
            };
            var element = new TriangleElement(nodes, "test");

            // Act
            var jacobian = element.CalculateJacobian(0.5, 0.5, 0.0);

            // Assert
            Assert.Equal(0.5, jacobian, 6);
        }

        [Fact]
        public void Element_CalculateNaturalCoordinates_ReturnsCorrectValues()
        {
            // Arrange
            var nodes = new[]
            {
                new Node(0, new[] { 0.0, 0.0, 0.0 }),
                new Node(1, new[] { 1.0, 0.0, 0.0 }),
                new Node(2, new[] { 0.0, 1.0, 0.0 })
            };
            var element = new TriangleElement(nodes, "test");
            var point = new[] { 0.5, 0.5, 0.0 };

            // Act
            var naturalCoordinates = element.CalculateNaturalCoordinates(point);

            // Assert
            Assert.Equal(3, naturalCoordinates.Length);
            Assert.Equal(0.0, naturalCoordinates[0], 6);
            Assert.Equal(0.5, naturalCoordinates[1], 6);
            Assert.Equal(0.5, naturalCoordinates[2], 6);
        }
    }
} 