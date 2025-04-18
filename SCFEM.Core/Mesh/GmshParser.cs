using System;
using System.Collections.Generic;
using System.IO;
using System.Numerics;
using SCFEM.Core.Mesh.Elements;

namespace SCFEM.Core.Mesh
{
    public class GmshParser
    {
        public static Mesh Parse(string filePath)
        {
            var mesh = new Mesh();
            var physicalGroups = new Dictionary<string, List<int>>();
            var elementTags = new Dictionary<int, string>();

            using (var reader = new StreamReader(filePath))
            {
                string line;
                while ((line = reader.ReadLine()) != null)
                {
                    line = line.Trim();
                    if (string.IsNullOrEmpty(line)) continue;

                    switch (line)
                    {
                        case "$MeshFormat":
                            ParseMeshFormat(reader);
                            break;
                        case "$PhysicalNames":
                            ParsePhysicalNames(reader, physicalGroups);
                            break;
                        case "$Nodes":
                            ParseNodes(reader, mesh);
                            break;
                        case "$Elements":
                            ParseElements(reader, mesh, elementTags, physicalGroups);
                            break;
                    }
                }
            }

            return mesh;
        }

        private static void ParseMeshFormat(StreamReader reader)
        {
            // Version number, file-type, data-size
            var formatLine = reader.ReadLine();
            var endSection = reader.ReadLine(); // $EndMeshFormat
        }

        private static void ParsePhysicalNames(StreamReader reader, Dictionary<string, List<int>> physicalGroups)
        {
            var numPhysicalNames = int.Parse(reader.ReadLine());
            for (int i = 0; i < numPhysicalNames; i++)
            {
                var parts = reader.ReadLine().Split(' ');
                var dimension = int.Parse(parts[0]);
                var tag = int.Parse(parts[1]);
                var name = parts[2].Trim('"');
                physicalGroups[name] = new List<int>();
            }
            var endSection = reader.ReadLine(); // $EndPhysicalNames
        }

        private static void ParseNodes(StreamReader reader, Mesh mesh)
        {
            var numNodes = int.Parse(reader.ReadLine());
            for (int i = 0; i < numNodes; i++)
            {
                var parts = reader.ReadLine().Split(' ');
                var x = float.Parse(parts[1]);
                var y = float.Parse(parts[2]);
                var z = float.Parse(parts[3]);
                mesh.Nodes.Add(new Vector3(x, y, z));
            }
            var endSection = reader.ReadLine(); // $EndNodes
        }

        private static void ParseElements(StreamReader reader, Mesh mesh, Dictionary<int, string> elementTags, Dictionary<string, List<int>> physicalGroups)
        {
            var numElements = int.Parse(reader.ReadLine());
            for (int i = 0; i < numElements; i++)
            {
                var parts = reader.ReadLine().Split(' ');
                var elementType = int.Parse(parts[1]);
                var numTags = int.Parse(parts[2]);
                var physicalTag = int.Parse(parts[3]);
                var elementTag = int.Parse(parts[4]);

                Element element = null;
                var nodeIndices = new int[GetNumNodesForType(elementType)];
                for (int j = 0; j < nodeIndices.Length; j++)
                {
                    nodeIndices[j] = int.Parse(parts[5 + numTags + j]) - 1; // Convert to 0-based index
                }

                switch ((ElementType)elementType)
                {
                    case ElementType.Tetrahedron:
                        element = new TetrahedronElement(nodeIndices);
                        break;
                    case ElementType.Hexahedron:
                        element = new HexahedronElement(nodeIndices);
                        break;
                    case ElementType.Prism:
                        element = new PrismElement(nodeIndices);
                        break;
                }

                if (element != null)
                {
                    mesh.Elements.Add(element);
                }
            }
            var endSection = reader.ReadLine(); // $EndElements
        }

        private static int GetNumNodesForType(int elementType)
        {
            return elementType switch
            {
                4 => 4, // Tetrahedron
                5 => 8, // Hexahedron
                6 => 6, // Prism
                _ => throw new ArgumentException($"Unsupported element type: {elementType}")
            };
        }
    }
} 