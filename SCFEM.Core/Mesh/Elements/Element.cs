using System;
using System.Numerics;

namespace SCFEM.Core.Mesh.Elements
{
    public abstract class Element
    {
        public int[] NodeIndices { get; protected set; }
        public ElementType Type { get; protected set; }
        public string PhysicalGroup { get; set; }

        protected Element(int[] nodeIndices, ElementType type)
        {
            NodeIndices = nodeIndices;
            Type = type;
        }

        public abstract double[,] CalculateStiffnessMatrix(Vector3[] nodes, double conductivity);
        public abstract double[] CalculateShapeFunctions(Vector3 point, Vector3[] nodes);
        public abstract Vector3[] CalculateShapeFunctionGradients(Vector3 point, Vector3[] nodes);
        public abstract double CalculateJacobian(Vector3[] nodes);
    }
} 