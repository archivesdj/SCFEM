using System;
using System.Linq;
using SCFEM.Core.Matrices;

namespace SCFEM.Core.Meshes.Elements
{
    public class TetrahedronElement : Element
    {
        // Gauss points for tetrahedral elements
        private static readonly double[][] GaussPoints = new double[][]
        {
            new double[] { 0.585410196624969, 0.138196601125011, 0.138196601125011 },
            new double[] { 0.138196601125011, 0.585410196624969, 0.138196601125011 },
            new double[] { 0.138196601125011, 0.138196601125011, 0.585410196624969 },
            new double[] { 0.138196601125011, 0.138196601125011, 0.138196601125011 }
        };
        private static readonly double[] GaussWeights = { 1.0/24.0, 1.0/24.0, 1.0/24.0, 1.0/24.0 };

        public TetrahedronElement(Node[] nodes, string physicalGroup) : base(ElementType.Tetrahedron, nodes, physicalGroup)
        {
        }

        public override double[,] CalculateStiffnessMatrix(MaterialProperties materialProperties)
        {
            var stiffnessMatrix = new double[4, 4]; // 4 nodes * 1 DOF
            var conductivity = materialProperties.Conductivity;

            // For 3D elements, the stiffness matrix is:
            // k = ∫∫∫ B^T D B |J| dξdηdζ
            // where B is the shape function gradient matrix
            // and D = [k 0 0; 0 k 0; 0 0 k]
            
            // Loop over Gauss points
            for (int i = 0; i < GaussPoints.Length; i++)
            {
                var xi = GaussPoints[i][0];
                var eta = GaussPoints[i][1];
                var zeta = GaussPoints[i][2];
                var weight = GaussWeights[i];
                
                // Get shape function gradients
                var gradients = CalculateShapeFunctionGradients(xi, eta, zeta);
                
                // Get Jacobian determinant
                var jacobian = CalculateJacobian(xi, eta, zeta);
                
                // Calculate B matrix
                var B = new double[3, 4];
                B[0, 0] = gradients[0];  // dN1/dx
                B[0, 1] = gradients[3];  // dN2/dx
                B[0, 2] = gradients[6];  // dN3/dx
                B[0, 3] = gradients[9];  // dN4/dx
                B[1, 0] = gradients[1];  // dN1/dy
                B[1, 1] = gradients[4];  // dN2/dy
                B[1, 2] = gradients[7];  // dN3/dy
                B[1, 3] = gradients[10]; // dN4/dy
                B[2, 0] = gradients[2];  // dN1/dz
                B[2, 1] = gradients[5];  // dN2/dz
                B[2, 2] = gradients[8];  // dN3/dz
                B[2, 3] = gradients[11]; // dN4/dz
                
                // Calculate B^T D B
                for (int j = 0; j < 4; j++)
                {
                    for (int k = 0; k < 4; k++)
                    {
                        var value = 0.0;
                        for (int l = 0; l < 3; l++)
                        {
                            value += B[l, j] * B[l, k];
                        }
                        stiffnessMatrix[j, k] += value * conductivity * weight * jacobian;
                    }
                }
            }

            return stiffnessMatrix;
        }

        public override double[] CalculateShapeFunctions(double xi, double eta, double zeta)
        {
            return new[]
            {
                1 - xi - eta - zeta,
                xi,
                eta,
                zeta
            };
        }

        public override double[] CalculateShapeFunctionGradients(double xi, double eta, double zeta)
        {
            var volume = CalculateVolume();
            var gradients = new double[12];
            
            // Calculate geometric coefficients
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
            
            var b1 = (y2 - y4) * (z3 - z4) - (y3 - y4) * (z2 - z4);
            var b2 = (y3 - y4) * (z1 - z4) - (y1 - y4) * (z3 - z4);
            var b3 = (y1 - y4) * (z2 - z4) - (y2 - y4) * (z1 - z4);
            var b4 = (y2 - y3) * (z1 - z4) - (y1 - y4) * (z2 - z3);
            
            var c1 = (x3 - x4) * (z2 - z4) - (x2 - x4) * (z3 - z4);
            var c2 = (x1 - x4) * (z3 - z4) - (x3 - x4) * (z1 - z4);
            var c3 = (x2 - x4) * (z1 - z4) - (x1 - x4) * (z2 - z4);
            var c4 = (x1 - x4) * (z2 - z3) - (x2 - x3) * (z1 - z4);
            
            var d1 = (x2 - x4) * (y3 - y4) - (x3 - x4) * (y2 - y4);
            var d2 = (x3 - x4) * (y1 - y4) - (x1 - x4) * (y3 - y4);
            var d3 = (x1 - x4) * (y2 - y4) - (x2 - x4) * (y1 - y4);
            var d4 = (x2 - x3) * (y1 - y4) - (x1 - x4) * (y2 - y3);
            
            // dN/dx = b/(6V), dN/dy = c/(6V), dN/dz = d/(6V)
            gradients[0] = b1 / (6 * volume);  // dN1/dx
            gradients[1] = c1 / (6 * volume);  // dN1/dy
            gradients[2] = d1 / (6 * volume);  // dN1/dz
            gradients[3] = b2 / (6 * volume);  // dN2/dx
            gradients[4] = c2 / (6 * volume);  // dN2/dy
            gradients[5] = d2 / (6 * volume);  // dN2/dz
            gradients[6] = b3 / (6 * volume);  // dN3/dx
            gradients[7] = c3 / (6 * volume);  // dN3/dy
            gradients[8] = d3 / (6 * volume);  // dN3/dz
            gradients[9] = b4 / (6 * volume);  // dN4/dx
            gradients[10] = c4 / (6 * volume);  // dN4/dy
            gradients[11] = d4 / (6 * volume);  // dN4/dz
            
            return gradients;
        }

        public override double CalculateJacobian(double xi, double eta, double zeta)
        {
            return CalculateVolume();
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

            var volume = CalculateVolume();
            var volume1 = CalculateVolume(x, y, z, x2, y2, z2, x3, y3, z3, x4, y4, z4);
            var volume2 = CalculateVolume(x1, y1, z1, x, y, z, x3, y3, z3, x4, y4, z4);
            var volume3 = CalculateVolume(x1, y1, z1, x2, y2, z2, x, y, z, x4, y4, z4);
            var volume4 = CalculateVolume(x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y, z);

            return new[]
            {
                volume1 / volume,
                volume2 / volume,
                volume3 / volume,
                volume4 / volume
            };
        }

        private double CalculateVolume()
        {
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

            return CalculateVolume(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);
        }

        private double CalculateVolume(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double x4, double y4, double z4)
        {
            var a = x2 - x1;
            var b = y2 - y1;
            var c = z2 - z1;
            var d = x3 - x1;
            var e = y3 - y1;
            var f = z3 - z1;
            var g = x4 - x1;
            var h = y4 - y1;
            var i = z4 - z1;

            return Math.Abs(a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)) / 6;
        }

        private double[,] CalculateJacobianMatrix(double xi, double eta, double zeta)
        {
            var J = new double[3, 3];
            var dN_dxi = new[]
            {
                -1,
                1,
                0,
                0
            };

            var dN_deta = new[]
            {
                -1,
                0,
                1,
                0
            };

            var dN_dzeta = new[]
            {
                -1,
                0,
                0,
                1
            };

            for (int i = 0; i < 4; i++)
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
            var volume = CalculateVolume();
            var B = new double[3, 4];
            
            // Calculate geometric coefficients
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
            
            var b1 = (y2 - y4) * (z3 - z4) - (y3 - y4) * (z2 - z4);
            var b2 = (y3 - y4) * (z1 - z4) - (y1 - y4) * (z3 - z4);
            var b3 = (y1 - y4) * (z2 - z4) - (y2 - y4) * (z1 - z4);
            var b4 = (y2 - y3) * (z1 - z4) - (y1 - y4) * (z2 - z3);
            
            var c1 = (x3 - x4) * (z2 - z4) - (x2 - x4) * (z3 - z4);
            var c2 = (x1 - x4) * (z3 - z4) - (x3 - x4) * (z1 - z4);
            var c3 = (x2 - x4) * (z1 - z4) - (x1 - x4) * (z2 - z4);
            var c4 = (x1 - x4) * (z2 - z3) - (x2 - x3) * (z1 - z4);
            
            var d1 = (x2 - x4) * (y3 - y4) - (x3 - x4) * (y2 - y4);
            var d2 = (x3 - x4) * (y1 - y4) - (x1 - x4) * (y3 - y4);
            var d3 = (x1 - x4) * (y2 - y4) - (x2 - x4) * (y1 - y4);
            var d4 = (x2 - x3) * (y1 - y4) - (x1 - x4) * (y2 - y3);
            
            // dN/dx = b/(6V), dN/dy = c/(6V), dN/dz = d/(6V)
            // Note: b, c, and d already contain the volume factor
            B[0, 0] = b1;  // dN1/dx
            B[0, 1] = b2;  // dN2/dx
            B[0, 2] = b3;  // dN3/dx
            B[0, 3] = b4;  // dN4/dx
            B[1, 0] = c1;  // dN1/dy
            B[1, 1] = c2;  // dN2/dy
            B[1, 2] = c3;  // dN3/dy
            B[1, 3] = c4;  // dN4/dy
            B[2, 0] = d1;  // dN1/dz
            B[2, 1] = d2;  // dN2/dz
            B[2, 2] = d3;  // dN3/dz
            B[2, 3] = d4;  // dN4/dz
            
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