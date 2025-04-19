using System;

namespace SCFEM.Core.Meshes
{
    public class Node
    {
        public int Id { get; set; }
        public double[] Coordinates { get; set; }

        public Node(int id, double[] coordinates)
        {
            Id = id;
            Coordinates = coordinates;
        }
    }
} 