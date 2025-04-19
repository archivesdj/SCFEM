using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using SCFEM.Core.Meshes.Elements;
using SCFEM.Core.Solver.BoundaryConditions;

namespace SCFEM.Core.Meshes
{
    public class Mesh
    {
        public List<Node> Nodes { get; }
        public List<Element> Elements { get; }
        public Dictionary<string, HashSet<Element>> PhysicalGroups { get; }
        public List<BoundaryCondition> BoundaryConditions { get; }

        public Dictionary<int, string> PhysicalNames { get; }

        public Mesh()
        {
            Nodes = new List<Node>();
            Elements = new List<Element>();
            PhysicalGroups = new Dictionary<string, HashSet<Element>>();
            BoundaryConditions = new List<BoundaryCondition>();
            PhysicalNames = new Dictionary<int, string>();
        }

        public void LoadFromGmsh(string filePath)
        {
            var parser = new GmshParser(this);
            parser.Parse(filePath);
        }

        public void Clear()
        {
            Nodes.Clear();
            Elements.Clear();
            PhysicalGroups.Clear();
            BoundaryConditions.Clear();
            PhysicalNames.Clear();
        }

        public List<Element> GetElementsByPhysicalGroup(string name)
        {
            if (PhysicalGroups.TryGetValue(name, out HashSet<Element>? elements))
            {
                return elements.ToList();
            }
            return new List<Element>();
        }
    }
} 