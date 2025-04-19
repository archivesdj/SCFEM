using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace SCFEM.Core.Meshes
{
    public class MaterialProperties
    {
        public Dictionary<string, double> Properties { get; private set; } = new Dictionary<string, double>();

        public double Conductivity
        {
            get => GetProperty("conductivity");
            set => Properties["conductivity"] = value;
        }

        public MaterialProperties()
        {
        }

        public MaterialProperties(Dictionary<string, double> properties)
        {
            Properties = properties;
        }

        public static MaterialProperties LoadFromFile(string filePath)
        {
            var properties = new Dictionary<string, double>();
            var lines = File.ReadAllLines(filePath);
            
            foreach (var line in lines)
            {
                var parts = line.Split('=');
                if (parts.Length == 2)
                {
                    var key = parts[0].Trim();
                    if (double.TryParse(parts[1].Trim(), out double value))
                    {
                        properties[key] = value;
                    }
                }
            }

            return new MaterialProperties(properties);
        }

        public double GetProperty(string key)
        {
            return Properties.TryGetValue(key, out double value) ? value : 0.0;
        }
    }
} 