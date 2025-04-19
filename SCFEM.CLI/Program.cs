using System;
using System.IO;
using SCFEM.Core.Mesh;
using SCFEM.Core.Solver;

namespace SCFEM.CLI
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args.Length != 2)
            {
                Console.WriteLine("Usage: SCFEM.CLI <mesh_file> <materials_file>");
                return;
            }

            string meshFile = args[0];
            string materialsFile = args[1];

            if (!File.Exists(meshFile))
            {
                Console.WriteLine($"Mesh file not found: {meshFile}");
                return;
            }

            if (!File.Exists(materialsFile))
            {
                Console.WriteLine($"Materials file not found: {materialsFile}");
                return;
            }

            try
            {
                var parser = new GmshParser();
                var mesh = parser.Parse(meshFile);
                mesh.LoadMaterialProperties(materialsFile);

                var solver = new StationaryCurrentSolver(mesh);
                solver.Solve();

                string outputFile = Path.ChangeExtension(meshFile, ".vtk");
                solver.SaveSolution(outputFile);
                Console.WriteLine($"Solution saved to {outputFile}");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Error: {ex.Message}");
            }
        }
    }
}
