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
            if (args.Length < 2)
            {
                Console.WriteLine("Usage: SCFEM.CLI <mesh_file> <material_properties_file>");
                return;
            }

            string meshFile = args[0];
            string materialFile = args[1];

            if (!File.Exists(meshFile))
            {
                Console.WriteLine($"Error: Mesh file {meshFile} does not exist.");
                return;
            }

            if (!File.Exists(materialFile))
            {
                Console.WriteLine($"Error: Material properties file {materialFile} does not exist.");
                return;
            }

            try
            {
                var mesh = new Mesh();
                mesh.LoadFromGmsh(meshFile);
                mesh.MaterialProperties = MaterialProperties.LoadFromFile(materialFile);

                var solver = new StationaryCurrentSolver(mesh);
                solver.AssembleSystem();
                solver.ApplyBoundaryConditions();
                solver.Solve();

                string outputFile = Path.ChangeExtension(meshFile, ".vtk");
                solver.ExportToVTK(outputFile);
                Console.WriteLine($"Solution exported to {outputFile}");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Error: {ex.Message}");
            }
        }
    }
}
