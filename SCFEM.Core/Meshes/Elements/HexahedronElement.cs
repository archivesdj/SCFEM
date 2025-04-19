using System;
using SCFEM.Core.Matrices;

namespace SCFEM.Core.Meshes.Elements
{
    public class HexahedronElement : Element
    {
        private static readonly double[] GaussPoints = { -1.0 / Math.Sqrt(3), 1.0 / Math.Sqrt(3) };
        private static readonly double[] GaussWeights = { 1.0, 1.0 };

        public HexahedronElement(Node[] nodes, string physicalGroup) : base(ElementType.Hexahedron, nodes, physicalGroup)
        {
        }

        public override double[,] CalculateStiffnessMatrix(MaterialProperties materialProperties)
        {
            var stiffnessMatrix = new double[8, 8]; // 8 nodes * 1 DOF

            for (int i = 0; i < GaussPoints.Length; i++)
            {
                var xi = GaussPoints[i];
                var eta = GaussPoints[i];
                var zeta = GaussPoints[i];
                var weight = GaussWeights[i];

                var B = CalculateBMatrix(xi, eta, zeta);
                var D = CalculateConductivityMatrix(materialProperties);
                var J = CalculateJacobian(xi, eta, zeta);

                var ke = Multiply(Transpose(B), Multiply(D, B));
                ke = Multiply(ke, J * weight);
                stiffnessMatrix = Add(stiffnessMatrix, ke);
            }

            return stiffnessMatrix;
        }

        public override double[] CalculateShapeFunctions(double xi, double eta, double zeta)
        {
            return new[]
            {
                (1 - xi) * (1 - eta) * (1 - zeta) / 8,
                (1 + xi) * (1 - eta) * (1 - zeta) / 8,
                (1 + xi) * (1 + eta) * (1 - zeta) / 8,
                (1 - xi) * (1 + eta) * (1 - zeta) / 8,
                (1 - xi) * (1 - eta) * (1 + zeta) / 8,
                (1 + xi) * (1 - eta) * (1 + zeta) / 8,
                (1 + xi) * (1 + eta) * (1 + zeta) / 8,
                (1 - xi) * (1 + eta) * (1 + zeta) / 8
            };
        }

        public override double[] CalculateShapeFunctionGradients(double xi, double eta, double zeta)
        {
            var J = CalculateJacobianMatrix(xi, eta, zeta);
            var JInv = Inverse(J);

            var dN_dxi = new[]
            {
                -(1 - eta) * (1 - zeta) / 8,
                (1 - eta) * (1 - zeta) / 8,
                (1 + eta) * (1 - zeta) / 8,
                -(1 + eta) * (1 - zeta) / 8,
                -(1 - eta) * (1 + zeta) / 8,
                (1 - eta) * (1 + zeta) / 8,
                (1 + eta) * (1 + zeta) / 8,
                -(1 + eta) * (1 + zeta) / 8
            };

            var dN_deta = new[]
            {
                -(1 - xi) * (1 - zeta) / 8,
                -(1 + xi) * (1 - zeta) / 8,
                (1 + xi) * (1 - zeta) / 8,
                (1 - xi) * (1 - zeta) / 8,
                -(1 - xi) * (1 + zeta) / 8,
                -(1 + xi) * (1 + zeta) / 8,
                (1 + xi) * (1 + zeta) / 8,
                (1 - xi) * (1 + zeta) / 8
            };

            var dN_dzeta = new[]
            {
                -(1 - xi) * (1 - eta) / 8,
                -(1 + xi) * (1 - eta) / 8,
                -(1 + xi) * (1 + eta) / 8,
                -(1 - xi) * (1 + eta) / 8,
                (1 - xi) * (1 - eta) / 8,
                (1 + xi) * (1 - eta) / 8,
                (1 + xi) * (1 + eta) / 8,
                (1 - xi) * (1 + eta) / 8
            };

            var gradients = new double[24];
            for (int i = 0; i < 8; i++)
            {
                gradients[i * 3] = JInv[0, 0] * dN_dxi[i] + JInv[0, 1] * dN_deta[i] + JInv[0, 2] * dN_dzeta[i];
                gradients[i * 3 + 1] = JInv[1, 0] * dN_dxi[i] + JInv[1, 1] * dN_deta[i] + JInv[1, 2] * dN_dzeta[i];
                gradients[i * 3 + 2] = JInv[2, 0] * dN_dxi[i] + JInv[2, 1] * dN_deta[i] + JInv[2, 2] * dN_dzeta[i];
            }

            return gradients;
        }

        public override double CalculateJacobian(double xi, double eta, double zeta)
        {
            var J = CalculateJacobianMatrix(xi, eta, zeta);
            return Determinant(J);
        }

        public override double[] CalculateNaturalCoordinates(double[] point)
        {
            var x = point[0];
            var y = point[1];
            var z = point[2];
            var x1 = Nodes[0].Coordinates[0];
            var y1 = Nodes[0].Coordinates[1];
            var z1 = Nodes[0].Coordinates[2];
            var x2 = Nodes[1].Coordinates[0];
            var y2 = Nodes[1].Coordinates[1];
            var z2 = Nodes[1].Coordinates[2];
            var x3 = Nodes[2].Coordinates[0];
            var y3 = Nodes[2].Coordinates[1];
            var z3 = Nodes[2].Coordinates[2];
            var x4 = Nodes[3].Coordinates[0];
            var y4 = Nodes[3].Coordinates[1];
            var z4 = Nodes[3].Coordinates[2];
            var x5 = Nodes[4].Coordinates[0];
            var y5 = Nodes[4].Coordinates[1];
            var z5 = Nodes[4].Coordinates[2];
            var x6 = Nodes[5].Coordinates[0];
            var y6 = Nodes[5].Coordinates[1];
            var z6 = Nodes[5].Coordinates[2];
            var x7 = Nodes[6].Coordinates[0];
            var y7 = Nodes[6].Coordinates[1];
            var z7 = Nodes[6].Coordinates[2];
            var x8 = Nodes[7].Coordinates[0];
            var y8 = Nodes[7].Coordinates[1];
            var z8 = Nodes[7].Coordinates[2];

            var xi = 2 * (x - x1) / (x2 - x1) - 1;
            var eta = 2 * (y - y1) / (y4 - y1) - 1;
            var zeta = 2 * (z - z1) / (z5 - z1) - 1;

            return new[]
            {
                (1 - xi) * (1 - eta) * (1 - zeta) / 8,
                (1 + xi) * (1 - eta) * (1 - zeta) / 8,
                (1 + xi) * (1 + eta) * (1 - zeta) / 8,
                (1 - xi) * (1 + eta) * (1 - zeta) / 8,
                (1 - xi) * (1 - eta) * (1 + zeta) / 8,
                (1 + xi) * (1 - eta) * (1 + zeta) / 8,
                (1 + xi) * (1 + eta) * (1 + zeta) / 8,
                (1 - xi) * (1 + eta) * (1 + zeta) / 8
            };
        }

        private double[,] CalculateJacobianMatrix(double xi, double eta, double zeta)
        {
            var J = new double[3, 3];
            var dN_dxi = new[]
            {
                -(1 - eta) * (1 - zeta) / 8,
                (1 - eta) * (1 - zeta) / 8,
                (1 + eta) * (1 - zeta) / 8,
                -(1 + eta) * (1 - zeta) / 8,
                -(1 - eta) * (1 + zeta) / 8,
                (1 - eta) * (1 + zeta) / 8,
                (1 + eta) * (1 + zeta) / 8,
                -(1 + eta) * (1 + zeta) / 8
            };

            var dN_deta = new[]
            {
                -(1 - xi) * (1 - zeta) / 8,
                -(1 + xi) * (1 - zeta) / 8,
                (1 + xi) * (1 - zeta) / 8,
                (1 - xi) * (1 - zeta) / 8,
                -(1 - xi) * (1 + zeta) / 8,
                -(1 + xi) * (1 + zeta) / 8,
                (1 + xi) * (1 + zeta) / 8,
                (1 - xi) * (1 + zeta) / 8
            };

            var dN_dzeta = new[]
            {
                -(1 - xi) * (1 - eta) / 8,
                -(1 + xi) * (1 - eta) / 8,
                -(1 + xi) * (1 + eta) / 8,
                -(1 - xi) * (1 + eta) / 8,
                (1 - xi) * (1 - eta) / 8,
                (1 + xi) * (1 - eta) / 8,
                (1 + xi) * (1 + eta) / 8,
                (1 - xi) * (1 + eta) / 8
            };

            for (int i = 0; i < 8; i++)
            {
                J[0, 0] += dN_dxi[i] * Nodes[i].Coordinates[0];
                J[0, 1] += dN_deta[i] * Nodes[i].Coordinates[0];
                J[0, 2] += dN_dzeta[i] * Nodes[i].Coordinates[0];
                J[1, 0] += dN_dxi[i] * Nodes[i].Coordinates[1];
                J[1, 1] += dN_deta[i] * Nodes[i].Coordinates[1];
                J[1, 2] += dN_dzeta[i] * Nodes[i].Coordinates[1];
                J[2, 0] += dN_dxi[i] * Nodes[i].Coordinates[2];
                J[2, 1] += dN_deta[i] * Nodes[i].Coordinates[2];
                J[2, 2] += dN_dzeta[i] * Nodes[i].Coordinates[2];
            }

            return J;
        }

        private double[,] CalculateBMatrix(double xi, double eta, double zeta)
        {
            var gradients = CalculateShapeFunctionGradients(xi, eta, zeta);
            var B = new double[3, 8];

            for (int i = 0; i < 8; i++)
            {
                B[0, i] = gradients[i * 3];
                B[1, i] = gradients[i * 3 + 1];
                B[2, i] = gradients[i * 3 + 2];
            }

            return B;
        }

        private double[,] CalculateConductivityMatrix(MaterialProperties materialProperties)
        {
            var conductivity = materialProperties.Conductivity;
            return new double[,]
            {
                { conductivity, 0, 0 },
                { 0, conductivity, 0 },
                { 0, 0, conductivity }
            };
        }
    }
} 