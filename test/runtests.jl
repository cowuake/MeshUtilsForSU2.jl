using MeshUtilsForSU2
using Test

#@testset "Mesh reading" begin
#    @test readMesh( "./in/NACA0012_sharp.su2" ).nElements == 52500
#    @test readMesh( "./in/NACA0012_sharp.su2" ).nMarkers == 2
#    @test readMesh( "./in/NACA0012_sharp.su2" ).nPoints == 52850
#end

@testset "Mesh reading" begin
    @test readMesh( "./in/NACA0012_sharp.su2" ).nElements == 52500 &&
        readMesh( "./in/NACA0012_sharp.su2" ).nMarkers == 2 &&
        readMesh( "./in/NACA0012_sharp.su2" ).nPoints == 52850
end

@testset "Refinement of structured grids --- 2D" begin
    @test refineStructuredMesh(readMesh( "./in/NACA0012_sharp.su2" )).nElements == 210000
    #@test refineStructuredMesh(readMesh( "./in/NACA0012_sharp.su2" )).nMarkers == 2
    #@test refineStructuredMesh(readMesh( "./in/NACA0012_sharp.su2" )).nPoints == 315350
end

@testset "Refinement of structured grids --- 3D" begin
    @test typeof( writeMesh( refineStructuredMesh( readMesh( "./in/cube.su2" ) , 4 ) ,
                  "./out/cube_ref4.su2" ) ) == Nothing
end

@testset "Laplacian smoothing" begin
    @test typeof( writeMesh( laplacianSmoothing2D( readMesh( "./in/NACA0012_blunt.su2" ) , 5 ) ,
                       "./out/NACA0012_blunt_ref1_smooth5.su2" ) ) == Nothing
end
