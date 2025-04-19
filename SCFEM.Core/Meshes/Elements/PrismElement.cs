using System;
using SCFEM.Core.Matrices;

namespace SCFEM.Core.Meshes.Elements
{
    public class PrismElement : Element
    {
        private static readonly double[] GaussPoints = { -1.0 / Math.Sqrt(3), 1.0 / Math.Sqrt(3) };
        private static readonly double[] GaussWeights = { 1.0, 1.0 };

        public PrismElement(Node[] nodes, string physicalGroup) : base(ElementType.Prism, nodes, physicalGroup)
        {
        }

        public override double[,] CalculateStiffnessMatrix(MaterialProperties materialProperties)
        {
            var stiffnessMatrix = new double[6, 6]; // 6 nodes * 1 DOF

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
                (1 - xi - eta) * (1 - zeta) / 2,
                xi * (1 - zeta) / 2,
                eta * (1 - zeta) / 2,
                (1 - xi - eta) * (1 + zeta) / 2,
                xi * (1 + zeta) / 2,
                eta * (1 + zeta) / 2
            };
        }

        public override double[] CalculateShapeFunctionGradients(double xi, double eta, double zeta)
        {
            var J = CalculateJacobianMatrix(xi, eta, zeta);
            var JInv = Inverse(J);

            var dN_dxi = new[]
            {
                -(1 - zeta) / 2,
                (1 - zeta) / 2,
                0,
                -(1 + zeta) / 2,
                (1 + zeta) / 2,
                0
            };

            var dN_deta = new[]
            {
                -(1 - zeta) / 2,
                0,
                (1 - zeta) / 2,
                -(1 + zeta) / 2,
                0,
                (1 + zeta) / 2
            };

            var dN_dzeta = new[]
            {
                -(1 - xi - eta) / 2,
                -xi / 2,
                -eta / 2,
                (1 - xi - eta) / 2,
                xi / 2,
                eta / 2
            };

            var gradients = new double[18];
            for (int i = 0; i < 6; i++)
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

            var area = CalculateArea();
            var height = CalculateHeight();
            var area1 = CalculateArea(x, y, x2, y2, x3, y3);
            var area2 = CalculateArea(x1, y1, x, y, x3, y3);
            var area3 = CalculateArea(x1, y1, x2, y2, x, y);
            var zeta = 2 * (z - z1) / height - 1;

            return new[]
            {
                area1 / area * (1 - zeta) / 2,
                area2 / area * (1 - zeta) / 2,
                area3 / area * (1 - zeta) / 2,
                area1 / area * (1 + zeta) / 2,
                area2 / area * (1 + zeta) / 2,
                area3 / area * (1 + zeta) / 2
            };
        }

        private double CalculateArea()
        {
            var x1 = Nodes[0].Coordinates[0];
            var y1 = Nodes[0].Coordinates[1];
            var x2 = Nodes[1].Coordinates[0];
            var y2 = Nodes[1].Coordinates[1];
            var x3 = Nodes[2].Coordinates[0];
            var y3 = Nodes[2].Coordinates[1];

            return CalculateArea(x1, y1, x2, y2, x3, y3);
        }

        private double CalculateArea(double x1, double y1, double x2, double y2, double x3, double y3)
        {
            return Math.Abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) / 2;
        }

        private double CalculateHeight()
        {
            var z1 = Nodes[0].Coordinates[2];
            var z4 = Nodes[3].Coordinates[2];
            return Math.Abs(z4 - z1);
        }

        private double[,] CalculateJacobianMatrix(double xi, double eta, double zeta)
        {
            var J = new double[3, 3];
            var dN_dxi = new[]
            {
                -(1 - zeta) / 2,
                (1 - zeta) / 2,
                0,
                -(1 + zeta) / 2,
                (1 + zeta) / 2,
                0
            };

            var dN_deta = new[]
            {
                -(1 - zeta) / 2,
                0,
                (1 - zeta) / 2,
                -(1 + zeta) / 2,
                0,
                (1 + zeta) / 2
            };

            var dN_dzeta = new[]
            {
                -(1 - xi - eta) / 2,
                -xi / 2,
                -eta / 2,
                (1 - xi - eta) / 2,
                xi / 2,
                eta / 2
            };

            for (int i = 0; i < 6; i++)
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
            var B = new double[3, 6];

            for (int i = 0; i < 6; i++)
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