using System;
using System.Linq;
using SCFEM.Core.Matrices;
using SCFEM.Core.Meshes;

namespace SCFEM.Core.Meshes.Elements
{
    public class TriangleElement : Element
    {
        // Gauss points for triangular elements
        private static readonly double[][] GaussPoints = new double[][]
        {
            new double[] { 1.0/6.0, 1.0/6.0 },
            new double[] { 2.0/3.0, 1.0/6.0 },
            new double[] { 1.0/6.0, 2.0/3.0 }
        };
        private static readonly double[] GaussWeights = { 1.0/6.0, 1.0/6.0, 1.0/6.0 };

        public TriangleElement(Node[] nodes, string physicalGroup) : base(ElementType.Triangle, nodes, physicalGroup)
        {
        }

        public override double[,] CalculateStiffnessMatrix(MaterialProperties materialProperties)
        {
            var stiffnessMatrix = new double[3, 3]; // 3 nodes * 1 DOF
            var conductivity = materialProperties.Conductivity;
            var thickness = 1.0; // Default thickness for 2D elements

            // For 2D elements, the stiffness matrix is:
            // k = ∫∫ B^T D B t |J| dξdη
            // where B is the shape function gradient matrix
            // and D = [k 0; 0 k]
            
            // Loop over Gauss points
            for (int i = 0; i < GaussPoints.Length; i++)
            {
                var xi = GaussPoints[i][0];
                var eta = GaussPoints[i][1];
                var weight = GaussWeights[i];
                
                // Get shape function gradients
                var gradients = CalculateShapeFunctionGradients(xi, eta, 0);
                
                // Get Jacobian determinant
                var jacobian = CalculateJacobian(xi, eta, 0);
                
                // Calculate B matrix
                var B = new double[2, 3];
                B[0, 0] = gradients[0]; // dN1/dx
                B[0, 1] = gradients[2]; // dN2/dx
                B[0, 2] = gradients[4]; // dN3/dx
                B[1, 0] = gradients[1]; // dN1/dy
                B[1, 1] = gradients[3]; // dN2/dy
                B[1, 2] = gradients[5]; // dN3/dy
                
                // Calculate B^T D B
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        var value = 0.0;
                        for (int l = 0; l < 2; l++)
                        {
                            value += B[l, j] * B[l, k];
                        }
                        stiffnessMatrix[j, k] += value * conductivity * thickness * weight * jacobian;
                    }
                }
            }

            return stiffnessMatrix;
        }

        public override double[] CalculateShapeFunctions(double xi, double eta, double zeta)
        {
            return new[]
            {
                1 - xi - eta,
                xi,
                eta
            };
        }

        public override double[] CalculateShapeFunctionGradients(double xi, double eta, double zeta)
        {
            var area = CalculateArea();
            var gradients = new double[6];
            
            // Calculate geometric coefficients
            var x1 = Nodes[0].Coordinates[0];
            var y1 = Nodes[0].Coordinates[1];
            var x2 = Nodes[1].Coordinates[0];
            var y2 = Nodes[1].Coordinates[1];
            var x3 = Nodes[2].Coordinates[0];
            var y3 = Nodes[2].Coordinates[1];
            
            var b1 = y2 - y3;
            var b2 = y3 - y1;
            var b3 = y1 - y2;
            var c1 = x3 - x2;
            var c2 = x1 - x3;
            var c3 = x2 - x1;
            
            // dN/dx = b/(2A), dN/dy = c/(2A)
            gradients[0] = b1 / (2 * area);  // dN1/dx
            gradients[1] = c1 / (2 * area);  // dN1/dy
            gradients[2] = b2 / (2 * area);  // dN2/dx
            gradients[3] = c2 / (2 * area);  // dN2/dy
            gradients[4] = b3 / (2 * area);  // dN3/dx
            gradients[5] = c3 / (2 * area);  // dN3/dy
            
            return gradients;
        }

        public override double CalculateJacobian(double xi, double eta, double zeta)
        {
            return CalculateArea();
        }

        public override double[] CalculateNaturalCoordinates(double[] point)
        {
            var x = point[0];
            var y = point[1];
            var x1 = Nodes[0].Coordinates[0];
            var y1 = Nodes[0].Coordinates[1];
            var x2 = Nodes[1].Coordinates[0];
            var y2 = Nodes[1].Coordinates[1];
            var x3 = Nodes[2].Coordinates[0];
            var y3 = Nodes[2].Coordinates[1];

            var area = CalculateArea();
            var area1 = CalculateArea(x, y, x2, y2, x3, y3);
            var area2 = CalculateArea(x1, y1, x, y, x3, y3);
            var area3 = CalculateArea(x1, y1, x2, y2, x, y);

            return new[]
            {
                area1 / area,
                area2 / area,
                area3 / area
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

        private double[,] CalculateJacobianMatrix(double xi, double eta)
        {
            var J = new double[2, 2];
            var dN_dxi = new[]
            {
                -1,
                1,
                0
            };

            var dN_deta = new[]
            {
                -1,
                0,
                1
            };

            for (int i = 0; i < 3; i++)
            {
                J[0, 0] += dN_dxi[i] * Nodes[i].Coordinates[0];
                J[0, 1] += dN_deta[i] * Nodes[i].Coordinates[0];
                J[1, 0] += dN_dxi[i] * Nodes[i].Coordinates[1];
                J[1, 1] += dN_deta[i] * Nodes[i].Coordinates[1];
            }

            return J;
        }

        private double[,] CalculateBMatrix(double xi, double eta)
        {
            var area = CalculateArea();
            var B = new double[2, 3];
            
            // Calculate geometric coefficients
            var x1 = Nodes[0].Coordinates[0];
            var y1 = Nodes[0].Coordinates[1];
            var x2 = Nodes[1].Coordinates[0];
            var y2 = Nodes[1].Coordinates[1];
            var x3 = Nodes[2].Coordinates[0];
            var y3 = Nodes[2].Coordinates[1];
            
            var b1 = y2 - y3;
            var b2 = y3 - y1;
            var b3 = y1 - y2;
            var c1 = x3 - x2;
            var c2 = x1 - x3;
            var c3 = x2 - x1;
            
            // dN/dx = b/(2A), dN/dy = c/(2A)
            // Note: b and c already contain the area factor
            B[0, 0] = b1;  // dN1/dx
            B[0, 1] = b2;  // dN2/dx
            B[0, 2] = b3;  // dN3/dx
            B[1, 0] = c1;  // dN1/dy
            B[1, 1] = c2;  // dN2/dy
            B[1, 2] = c3;  // dN3/dy
            
            return B;
        }

        private double[,] CalculateConductivityMatrix(MaterialProperties materialProperties)
        {
            var conductivity = materialProperties.Conductivity;
            return new double[,]
            {
                { conductivity, 0 },
                { 0, conductivity }
            };
        }
    }
} 