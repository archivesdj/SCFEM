using System;
using System.Collections.Generic;
using System.IO;

namespace SCFEM.Core.Mesh
{
    public static class MaterialProperties
    {
        public static Dictionary<string, double> LoadFromFile(string filePath)
        {
            var properties = new Dictionary<string, double>();

            using (var reader = new StreamReader(filePath))
            {
                string line;
                while ((line = reader.ReadLine()) != null)
                {
                    line = line.Trim();
                    if (string.IsNullOrEmpty(line) || line.StartsWith("#")) continue;

                    var parts = line.Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                    if (parts.Length != 2)
                    {
                        throw new FormatException($"Invalid line format in material properties file: {line}");
                    }

                    var partName = parts[0];
                    if (!double.TryParse(parts[1], out double value))
                    {
                        throw new FormatException($"Invalid conductivity value in material properties file: {parts[1]}");
                    }

                    properties[partName] = value;
                }
            }

            return properties;
        }
    }
} 