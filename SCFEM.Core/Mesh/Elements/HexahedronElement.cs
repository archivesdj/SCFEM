using System;
using System.Numerics;

namespace SCFEM.Core.Mesh.Elements
{
    public class HexahedronElement : Element
    {
        private static readonly double[] GaussPoints = { -1.0 / Math.Sqrt(3), 1.0 / Math.Sqrt(3) };
        private static readonly double[] GaussWeights = { 1.0, 1.0 };

        public HexahedronElement(int[] nodeIndices) : base(nodeIndices, ElementType.Hexahedron)
        {
            if (nodeIndices.Length != 8)
                throw new ArgumentException("Hexahedron element must have exactly 8 nodes");
        }

        public override double[,] CalculateStiffnessMatrix(Vector3[] nodes, double conductivity)
        {
            var stiffness = new double[8, 8];
            double totalVolume = 0.0;

            // 2x2x2 Gauss quadrature
            foreach (var xi in GaussPoints)
            {
                foreach (var eta in GaussPoints)
                {
                    foreach (var zeta in GaussPoints)
                    {
                        var point = new Vector3((float)xi, (float)eta, (float)zeta);
                        var gradients = CalculateShapeFunctionGradients(point, nodes);
                        var jacobian = CalculateJacobian(nodes);
                        var weight = GaussWeights[0] * GaussWeights[0] * GaussWeights[0] * Math.Abs(jacobian);

                        for (int i = 0; i < 8; i++)
                        {
                            for (int j = 0; j < 8; j++)
                            {
                                stiffness[i, j] += conductivity * weight * 
                                    Vector3.Dot(gradients[i], gradients[j]);
                            }
                        }
                        totalVolume += weight;
                    }
                }
            }

            return stiffness;
        }

        public override double[] CalculateShapeFunctions(Vector3 point, Vector3[] nodes)
        {
            double xi = point.X;
            double eta = point.Y;
            double zeta = point.Z;

            var N = new double[8];
            N[0] = (1 - xi) * (1 - eta) * (1 - zeta) / 8.0;
            N[1] = (1 + xi) * (1 - eta) * (1 - zeta) / 8.0;
            N[2] = (1 + xi) * (1 + eta) * (1 - zeta) / 8.0;
            N[3] = (1 - xi) * (1 + eta) * (1 - zeta) / 8.0;
            N[4] = (1 - xi) * (1 - eta) * (1 + zeta) / 8.0;
            N[5] = (1 + xi) * (1 - eta) * (1 + zeta) / 8.0;
            N[6] = (1 + xi) * (1 + eta) * (1 + zeta) / 8.0;
            N[7] = (1 - xi) * (1 + eta) * (1 + zeta) / 8.0;

            return N;
        }

        public override Vector3[] CalculateShapeFunctionGradients(Vector3 point, Vector3[] nodes)
        {
            double xi = point.X;
            double eta = point.Y;
            double zeta = point.Z;

            // Natural coordinate derivatives
            var dN_dxi = new double[8];
            var dN_deta = new double[8];
            var dN_dzeta = new double[8];

            dN_dxi[0] = -(1 - eta) * (1 - zeta) / 8.0;
            dN_dxi[1] = (1 - eta) * (1 - zeta) / 8.0;
            dN_dxi[2] = (1 + eta) * (1 - zeta) / 8.0;
            dN_dxi[3] = -(1 + eta) * (1 - zeta) / 8.0;
            dN_dxi[4] = -(1 - eta) * (1 + zeta) / 8.0;
            dN_dxi[5] = (1 - eta) * (1 + zeta) / 8.0;
            dN_dxi[6] = (1 + eta) * (1 + zeta) / 8.0;
            dN_dxi[7] = -(1 + eta) * (1 + zeta) / 8.0;

            dN_deta[0] = -(1 - xi) * (1 - zeta) / 8.0;
            dN_deta[1] = -(1 + xi) * (1 - zeta) / 8.0;
            dN_deta[2] = (1 + xi) * (1 - zeta) / 8.0;
            dN_deta[3] = (1 - xi) * (1 - zeta) / 8.0;
            dN_deta[4] = -(1 - xi) * (1 + zeta) / 8.0;
            dN_deta[5] = -(1 + xi) * (1 + zeta) / 8.0;
            dN_deta[6] = (1 + xi) * (1 + zeta) / 8.0;
            dN_deta[7] = (1 - xi) * (1 + zeta) / 8.0;

            dN_dzeta[0] = -(1 - xi) * (1 - eta) / 8.0;
            dN_dzeta[1] = -(1 + xi) * (1 - eta) / 8.0;
            dN_dzeta[2] = -(1 + xi) * (1 + eta) / 8.0;
            dN_dzeta[3] = -(1 - xi) * (1 + eta) / 8.0;
            dN_dzeta[4] = (1 - xi) * (1 - eta) / 8.0;
            dN_dzeta[5] = (1 + xi) * (1 - eta) / 8.0;
            dN_dzeta[6] = (1 + xi) * (1 + eta) / 8.0;
            dN_dzeta[7] = (1 - xi) * (1 + eta) / 8.0;

            // Calculate Jacobian matrix
            var J = new double[3, 3];
            for (int i = 0; i < 8; i++)
            {
                J[0, 0] += dN_dxi[i] * nodes[i].X;
                J[0, 1] += dN_dxi[i] * nodes[i].Y;
                J[0, 2] += dN_dxi[i] * nodes[i].Z;
                J[1, 0] += dN_deta[i] * nodes[i].X;
                J[1, 1] += dN_deta[i] * nodes[i].Y;
                J[1, 2] += dN_deta[i] * nodes[i].Z;
                J[2, 0] += dN_dzeta[i] * nodes[i].X;
                J[2, 1] += dN_dzeta[i] * nodes[i].Y;
                J[2, 2] += dN_dzeta[i] * nodes[i].Z;
            }

            // Calculate inverse of Jacobian
            var detJ = J[0, 0] * (J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1]) -
                      J[0, 1] * (J[1, 0] * J[2, 2] - J[1, 2] * J[2, 0]) +
                      J[0, 2] * (J[1, 0] * J[2, 1] - J[1, 1] * J[2, 0]);

            var invJ = new double[3, 3];
            invJ[0, 0] = (J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1]) / detJ;
            invJ[0, 1] = (J[0, 2] * J[2, 1] - J[0, 1] * J[2, 2]) / detJ;
            invJ[0, 2] = (J[0, 1] * J[1, 2] - J[0, 2] * J[1, 1]) / detJ;
            invJ[1, 0] = (J[1, 2] * J[2, 0] - J[1, 0] * J[2, 2]) / detJ;
            invJ[1, 1] = (J[0, 0] * J[2, 2] - J[0, 2] * J[2, 0]) / detJ;
            invJ[1, 2] = (J[0, 2] * J[1, 0] - J[0, 0] * J[1, 2]) / detJ;
            invJ[2, 0] = (J[1, 0] * J[2, 1] - J[1, 1] * J[2, 0]) / detJ;
            invJ[2, 1] = (J[0, 1] * J[2, 0] - J[0, 0] * J[2, 1]) / detJ;
            invJ[2, 2] = (J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]) / detJ;

            // Transform to global coordinates
            var gradients = new Vector3[8];
            for (int i = 0; i < 8; i++)
            {
                double dx = dN_dxi[i] * invJ[0, 0] + dN_deta[i] * invJ[1, 0] + dN_dzeta[i] * invJ[2, 0];
                double dy = dN_dxi[i] * invJ[0, 1] + dN_deta[i] * invJ[1, 1] + dN_dzeta[i] * invJ[2, 1];
                double dz = dN_dxi[i] * invJ[0, 2] + dN_deta[i] * invJ[1, 2] + dN_dzeta[i] * invJ[2, 2];
                gradients[i] = new Vector3((float)dx, (float)dy, (float)dz);
            }

            return gradients;
        }

        public override double CalculateJacobian(Vector3[] nodes)
        {
            var point = new Vector3(0, 0, 0); // Default point for Jacobian calculation
            var gradients = CalculateShapeFunctionGradients(point, nodes);
            var J = new double[3, 3];
            for (int i = 0; i < 8; i++)
            {
                J[0, 0] += gradients[i].X * nodes[i].X;
                J[0, 1] += gradients[i].X * nodes[i].Y;
                J[0, 2] += gradients[i].X * nodes[i].Z;
                J[1, 0] += gradients[i].Y * nodes[i].X;
                J[1, 1] += gradients[i].Y * nodes[i].Y;
                J[1, 2] += gradients[i].Y * nodes[i].Z;
                J[2, 0] += gradients[i].Z * nodes[i].X;
                J[2, 1] += gradients[i].Z * nodes[i].Y;
                J[2, 2] += gradients[i].Z * nodes[i].Z;
            }
            return J[0, 0] * (J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1]) -
                   J[0, 1] * (J[1, 0] * J[2, 2] - J[1, 2] * J[2, 0]) +
                   J[0, 2] * (J[1, 0] * J[2, 1] - J[1, 1] * J[2, 0]);
        }
    }
} 