using System;
using System.Linq;
using SCFEM.Core.Matrices;
using SCFEM.Core.Meshes;

namespace SCFEM.Core.Meshes.Elements
{
    public class LineElement : Element
    {
        private static readonly double[] GaussPoints = { -0.577350269189626, 0.577350269189626 };
        private static readonly double[] GaussWeights = { 1.0, 1.0 };

        public LineElement(Node[] nodes, string physicalGroup) : base(ElementType.Line, nodes, physicalGroup)
        {
        }

        public override double[,] CalculateStiffnessMatrix(MaterialProperties materialProperties)
        {
            var stiffnessMatrix = new double[2, 2]; // 2 nodes * 1 DOF
            var length = CalculateLength();
            var conductivity = materialProperties.Conductivity;

            // For 1D elements, the stiffness matrix is:
            // k = (A/L) * [ 1 -1 ]
            //             [ -1 1 ]
            stiffnessMatrix[0, 0] = conductivity / length;
            stiffnessMatrix[0, 1] = -conductivity / length;
            stiffnessMatrix[1, 0] = -conductivity / length;
            stiffnessMatrix[1, 1] = conductivity / length;

            return stiffnessMatrix;
        }

        public override double[] CalculateShapeFunctions(double xi, double eta, double zeta)
        {
            return new[]
            {
                (1 - xi) / 2,
                (1 + xi) / 2
            };
        }

        public override double[] CalculateShapeFunctionGradients(double xi, double eta, double zeta)
        {
            var length = CalculateLength();
            return new[]
            {
                -1.0 / length,  // dN1/dx
                1.0 / length    // dN2/dx
            };
        }

        public override double CalculateJacobian(double xi, double eta, double zeta)
        {
            return CalculateLength() / 2;
        }

        public override double[] CalculateNaturalCoordinates(double[] point)
        {
            var x = point[0];
            var x1 = Nodes[0].Coordinates[0];
            var x2 = Nodes[1].Coordinates[0];

            var length = CalculateLength();
            var xi = 2 * (x - x1) / length - 1;

            return new[] { xi };
        }

        private double CalculateLength()
        {
            var x1 = Nodes[0].Coordinates[0];
            var x2 = Nodes[1].Coordinates[0];

            return Math.Abs(x2 - x1);
        }

        private double[,] CalculateJacobianMatrix(double xi)
        {
            var J = new double[1, 1];
            var dN_dxi = new[]
            {
                -0.5,
                0.5
            };

            for (int i = 0; i < 2; i++)
            {
                J[0, 0] += dN_dxi[i] * Nodes[i].Coordinates[0];
            }

            return J;
        }

        private double[,] CalculateBMatrix(double xi)
        {
            var length = CalculateLength();
            var B = new double[1, 2];
            
            // dN/dx = (dN/dξ) * (dξ/dx)
            // dξ/dx = 2/L
            B[0, 0] = -1.0 / length;  // dN1/dx
            B[0, 1] = 1.0 / length;   // dN2/dx
            
            return B;
        }

        private double[,] CalculateConductivityMatrix(MaterialProperties materialProperties)
        {
            var conductivity = materialProperties.Conductivity;
            return new double[,]
            {
                { conductivity }
            };
        }
    }
} 