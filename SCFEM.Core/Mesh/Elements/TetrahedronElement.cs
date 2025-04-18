using System;
using System.Numerics;

namespace SCFEM.Core.Mesh.Elements
{
    public class TetrahedronElement : Element
    {
        public TetrahedronElement(int[] nodeIndices) : base(nodeIndices, ElementType.Tetrahedron)
        {
            if (nodeIndices.Length != 4)
                throw new ArgumentException("Tetrahedron element must have exactly 4 nodes");
        }

        public override double[,] CalculateStiffnessMatrix(Vector3[] nodes, double conductivity)
        {
            var gradients = CalculateShapeFunctionGradients(Vector3.Zero, nodes);
            var jacobian = CalculateJacobian(nodes);
            var volume = Math.Abs(jacobian) / 6.0;

            var stiffness = new double[4, 4];
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    stiffness[i, j] = conductivity * volume * 
                        Vector3.Dot(gradients[i], gradients[j]);
                }
            }
            return stiffness;
        }

        public override double[] CalculateShapeFunctions(Vector3 point, Vector3[] nodes)
        {
            var N = new double[4];
            N[0] = 1 - point.X - point.Y - point.Z;
            N[1] = point.X;
            N[2] = point.Y;
            N[3] = point.Z;
            return N;
        }

        public override Vector3[] CalculateShapeFunctionGradients(Vector3 point, Vector3[] nodes)
        {
            var gradients = new Vector3[4];
            gradients[0] = new Vector3(-1, -1, -1);
            gradients[1] = new Vector3(1, 0, 0);
            gradients[2] = new Vector3(0, 1, 0);
            gradients[3] = new Vector3(0, 0, 1);
            return gradients;
        }

        public override double CalculateJacobian(Vector3[] nodes)
        {
            var J = new double[3, 3];
            for (int i = 0; i < 3; i++)
            {
                J[0, i] = nodes[1][i] - nodes[0][i];
                J[1, i] = nodes[2][i] - nodes[0][i];
                J[2, i] = nodes[3][i] - nodes[0][i];
            }
            return J[0, 0] * (J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1]) -
                   J[0, 1] * (J[1, 0] * J[2, 2] - J[1, 2] * J[2, 0]) +
                   J[0, 2] * (J[1, 0] * J[2, 1] - J[1, 1] * J[2, 0]);
        }
    }
} 