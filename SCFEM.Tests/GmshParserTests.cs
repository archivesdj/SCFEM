using System;
using System.IO;
using Xunit;
using SCFEM.Core.Meshes;

namespace SCFEM.Tests
{
    public class GmshParserTests
    {
        private readonly string _meshFilePath;

        public GmshParserTests()
        {
            _meshFilePath = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "..", "..", "..", "..", "examples", "cube.msh");
        }

        [Fact]
        public void Parse_ValidMeshFile_ReturnsMeshWithNodes()
        {
            var mesh = new Mesh();
            var parser = new GmshParser(mesh);
            parser.Parse(_meshFilePath);
            Assert.NotNull(mesh.Nodes);
            Assert.NotEmpty(mesh.Nodes);
        }

        [Fact]
        public void Parse_ValidMeshFile_ReturnsMeshWithElements()
        {
            var mesh = new Mesh();
            var parser = new GmshParser(mesh);
            parser.Parse(_meshFilePath);
            Assert.NotNull(mesh.Elements);
            Assert.NotEmpty(mesh.Elements);
        }

        [Fact]
        public void Parse_ValidMeshFile_ReturnsMeshWithPhysicalGroups()
        {
            var mesh = new Mesh();
            var parser = new GmshParser(mesh);
            parser.Parse(_meshFilePath);
            Assert.NotNull(mesh.PhysicalGroups);
            Assert.NotEmpty(mesh.PhysicalGroups);
        }

        [Fact]
        public void Parse_ValidMeshFile_ReturnsMeshWithBoundaryConditions()
        {
            var mesh = new Mesh();
            var parser = new GmshParser(mesh);
            parser.Parse(_meshFilePath);
            Assert.NotNull(mesh.BoundaryConditions);
        }
    }
} 