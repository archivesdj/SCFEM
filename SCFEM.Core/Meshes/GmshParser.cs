using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Diagnostics;
using SCFEM.Core.Meshes.Elements;
using SCFEM.Core.Solver.BoundaryConditions;

namespace SCFEM.Core.Meshes
{
    public class GmshParser
    {
        private readonly List<Node> _nodes;
        private readonly List<Element> _elements;
        private readonly Dictionary<string, HashSet<Element>> _physicalGroups;
        private readonly List<BoundaryCondition> _boundaryConditions;
        private readonly Dictionary<int, string> _physicalNames;

        public GmshParser(Mesh mesh)
        {
            _nodes = mesh.Nodes;
            _elements = mesh.Elements;
            _physicalGroups = mesh.PhysicalGroups;
            _boundaryConditions = mesh.BoundaryConditions;
            _physicalNames = mesh.PhysicalNames;
        }

        public void Parse(string filePath)
        {
            if (!File.Exists(filePath))
            {
                throw new FileNotFoundException($"Mesh file not found: {filePath}");
            }

            Debug.WriteLine($"Parsing mesh file: {filePath}");
            var lines = File.ReadAllLines(filePath);
            var currentSection = string.Empty;
            var lineIndex = 0;

            // Check mesh format version
            while (lineIndex < lines.Length)
            {
                var line = lines[lineIndex].Trim();
                if (line.StartsWith("$MeshFormat"))
                {
                    lineIndex++;
                    var formatLine = lines[lineIndex].Trim();
                    var formatParts = formatLine.Split(' ');
                    Debug.WriteLine($"Mesh format: {formatLine}");
                    if (formatParts[0] != "2.2")
                    {
                        throw new Exception("Only Gmsh format version 2.2 is supported");
                    }
                    break;
                }
                lineIndex++;
            }

            lineIndex = 0; // Reset to start parsing from beginning

            while (lineIndex < lines.Length)
            {
                var line = lines[lineIndex].Trim();
                if (string.IsNullOrWhiteSpace(line))
                {
                    lineIndex++;
                    continue;
                }

                if (line.StartsWith("$"))
                {
                    if (line.StartsWith("$End"))
                    {
                        Debug.WriteLine($"End of section: {currentSection}");
                        currentSection = string.Empty;
                    }
                    else
                    {
                        currentSection = line.Trim('$');
                        Debug.WriteLine($"Start of section: {currentSection}");
                    }
                    lineIndex++;
                    continue;
                }

                switch (currentSection)
                {
                    case "PhysicalNames":
                        if (int.TryParse(line, out int numPhysicalNames))
                        {
                            Debug.WriteLine($"Found {numPhysicalNames} physical names");
                            for (int i = 0; i < numPhysicalNames; i++)
                            {
                                lineIndex++;
                                var parts = lines[lineIndex].Split(' ');
                                var tag = int.Parse(parts[1]);
                                var name = parts[2].Trim('"');
                                _physicalNames[tag] = name;
                                Debug.WriteLine($"Physical name: {tag} -> {name}");
                            }
                        }
                        break;
                    case "Nodes":
                        if (int.TryParse(line, out int numNodes))
                        {
                            Debug.WriteLine($"Found {numNodes} nodes");
                            for (int i = 0; i < numNodes; i++)
                            {
                                lineIndex++;
                                var parts = lines[lineIndex].Split(' ');
                                var id = int.Parse(parts[0]);
                                var x = double.Parse(parts[1]);
                                var y = double.Parse(parts[2]);
                                var z = double.Parse(parts[3]);
                                _nodes.Add(new Node(id, new[] { x, y, z }));
                                Debug.WriteLine($"Node {id}: ({x}, {y}, {z})");
                            }
                        }
                        break;
                    case "Elements":
                        if (int.TryParse(line, out int numElements))
                        {
                            Debug.WriteLine($"Found {numElements} elements");
                            for (int i = 0; i < numElements; i++)
                            {
                                lineIndex++;
                                var parts = lines[lineIndex].Split(' ');
                                var elementType = int.Parse(parts[1]);
                                var numTags = int.Parse(parts[2]);
                                var physicalTag = int.Parse(parts[3 + numTags]);
                                
                                Debug.WriteLine($"Element {i + 1}: type={elementType}, tags={numTags}, physical={physicalTag}");
                                
                                if (!_physicalNames.ContainsKey(physicalTag))
                                {
                                    Debug.WriteLine($"Skipping element with unknown physical tag: {physicalTag}");
                                    continue;
                                }

                                var nodeIndices = parts.Skip(3 + numTags + 1).Select(int.Parse).ToArray();
                                var elementNodes = nodeIndices.Select(idx => _nodes.FirstOrDefault(n => n.Id == idx)).Where(n => n != null).ToArray();
                                
                                if (elementNodes.Length != nodeIndices.Length)
                                {
                                    Debug.WriteLine($"Skipping element with invalid node indices");
                                    continue;
                                }

                                var physicalGroup = _physicalNames[physicalTag];
                                Element element = null;

                                switch (elementType)
                                {
                                    case 1: // Line
                                        element = new LineElement(elementNodes, physicalGroup);
                                        break;
                                    case 2: // Triangle
                                        element = new TriangleElement(elementNodes, physicalGroup);
                                        break;
                                    case 4: // Tetrahedron
                                        element = new TetrahedronElement(elementNodes, physicalGroup);
                                        break;
                                    case 5: // Hexahedron
                                        element = new HexahedronElement(elementNodes, physicalGroup);
                                        break;
                                    case 6: // Prism
                                        element = new PrismElement(elementNodes, physicalGroup);
                                        break;
                                }

                                if (element != null)
                                {
                                    _elements.Add(element);
                                    if (!_physicalGroups.ContainsKey(physicalGroup))
                                    {
                                        _physicalGroups[physicalGroup] = new HashSet<Element>();
                                    }
                                    _physicalGroups[physicalGroup].Add(element);
                                    Debug.WriteLine($"Added element to physical group: {physicalGroup}");
                                }
                            }
                        }
                        break;
                }
                lineIndex++;
            }

            Debug.WriteLine($"Parsing completed: {_nodes.Count} nodes, {_elements.Count} elements, {_physicalGroups.Count} physical groups");
        }
    }
} 