using System;
using System.Collections.Generic;
using System.Numerics;
using SCFEM.Core.Mesh.Elements;

namespace SCFEM.Core.Mesh
{
    public class Mesh
    {
        public List<Vector3> Nodes { get; set; }
        public List<Element> Elements { get; set; }
        public Dictionary<string, List<int>> PhysicalGroups { get; set; }
        public Dictionary<string, double> MaterialProperties { get; set; }

        public Mesh()
        {
            Nodes = new List<Vector3>();
            Elements = new List<Element>();
            PhysicalGroups = new Dictionary<string, List<int>>();
            MaterialProperties = new Dictionary<string, double>();
        }

        public void LoadFromGmsh(string filePath)
        {
            var parsedMesh = GmshParser.Parse(filePath);
            Nodes = parsedMesh.Nodes;
            Elements = parsedMesh.Elements;
            PhysicalGroups = parsedMesh.PhysicalGroups;
        }
    }
} 