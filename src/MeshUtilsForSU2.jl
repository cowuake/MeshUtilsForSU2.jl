module MeshUtilsForSU2

using ArgParse
using GraphPlot
using LightGraphs
using Parsers
using Printf

export readMesh, writeMesh, refineStructuredMesh, laplacianSmoothing2D

abstract type Point end

mutable struct Point2D <: Point
    x::Float64
    y::Float64
end

mutable struct Point3D <: Point
    x::Float64
    y::Float64
    z::Float64
end

# SU2 support the following element types:
# ========================================
# Line (ID: 3)
# Triangle (ID: 5)
# Quadrilateral (ID: 9)
# Tetrahedral (ID: 10)
# Hexahedral (ID: 12)
# Prism (ID: 13)
# Pyramid (ID: 14)
mutable struct Element
    id::Int64
    connectivity::Vector{Int64}
    boundary::Vector{Int64}
end

mutable struct MarkerElement
    id::Int64
    connectivity::Vector{Int64}
end

mutable struct Marker
    tag::String
    nElements::Int64
    elements::Vector{MarkerElement}
end

struct Mesh
    nDimensions::Int64
    nElements::Int64
    nPoints::Int64
    nMarkers::Int64
    elements::Vector{Element}
    points::Vector{Point}
    markers::Vector{Marker}
end

function parseCommandLine( )
    settings = ArgParseSettings( )
    settings.prog = "mesh.jl"
    settings.description = "PUT DESCRIPTION HERE"
    settings.suppress_warnings = true

    @add_arg_table settings begin
        "--addzone"
        help = "Add a mesh zone to the input mesh"
        arg_type = String
        default = ""
        action = :store_arg
        "--debug"
        help = "Turn on additional output for debuggin/testing purposes"
        action = :store_true
        "--refine"
        help = "Refine the mesh by the specified refinement level"
        arg_type = Int
        default = 0
        constant = 2
        action = :store_arg
        "--rewrite"
        help = "Read and write again the specified mesh"
        action = :store_true
        "mesh"
        help = "The (.su2) mesh file to be manipulated"
        required = true
        "--smooth"
        help = "Smooth the specified mesh"
        action = :store_true
    end

    # Return type is Array{String,Any}
    return parse_args( settings )
end

function numberOfEdges( elementId::Int64 , debug::Bool=false )
    dict = Dict( [ ( 3 , 1 ) , # line
                   ( 5 , 3 ) , # triangle
                   ( 9 , 4 ) , # quadrilateral
                   ( 10 , 6 ) , # tetrahedron
                   ( 12 , 12 ) , # hexahedron
                   ( 14 , 8 ) ] ) # pyramid
    return dict[ elementId ]
end

function numberOfPointsOnBoundaries( elementId:: Int64 , debug::Bool=false )
    dict = Dict( [ ( 5 , 2 ) , # triangle
                   ( 9 , 2 ) , # quadrilateral
                   ( 10 , 3 ) , # tetrahedron
                   ( 12 , 4 ) , # hexahedron
                   ( 14 , 3 ) ] ) # pyramid
    return dict[ elementId ]
end

function numberOfFaces( elementId::Int64 , debug::Bool=false )
    dict = Dict( [ ( 10 , 4 ) , # tetrahedron
                   ( 12 , 6 ) , # hexahedron
                   ( 14 , 5 ) ] ) # pyramid
    return dict[ elementId ]
end

function numberOfBoundaries( elementId::Int64 , nDimensions::Int64 , debug::Bool=false )
    if nDimensions == 2
        return numberOfEdges( elementId )
    else
        return numberOfFaces( elementId )
    end
end

function numberOfNodes( elementType::Int64 , debug::Bool=false )
    dict = Dict( [ ( 3 , 2 ) , # line
                   ( 5 , 3 ) , # triangle
                   ( 9 , 4 ) , # quadrilateral
                   ( 10 , 4 ) , # tetrahedron
                   ( 12 , 8 ) , # hexahedron
                   ( 14 , 5 ) ] ) # pyramid
    return dict[ elementType ]
end

function readMesh( meshFile::String , debug::Bool=false )
    rawText = readlines( open( meshFile ) )

    # Terrible, to be changes ASAP
    function fixLine( line::String , debug::Bool=false )
        fixed = lstrip( line )
        fixed = replace( fixed , r"\s+" => "\t" )
    end

    # Clean text from comments
    text = [ fixLine( line ) for line in rawText if !occursin( "%" , line ) ]

    # Store position of the dimension line and number of dimensions
    dimensionLine = findfirst( s -> occursin( "NDIME=" , s ) , text )
    nDimensions = parse( Int64 , split( text[ dimensionLine ] , r"[=]" )[ end ] )

    # Store position of the point block
    pointLine = findfirst( s -> occursin( "NPOIN=" , s ) , text )
    nPoints = parse( Int64 , split( text[ pointLine ] , r"[=]" )[ end ] )

    # Declare dimension-dependent arrays
    if nDimensions == 2
        points = Vector{Point2D}( undef , nPoints )
    else
        points = Vector{Point3D}( undef , nPoints )
    end

    # Retrieve data from the point block and store in a vector
    if nDimensions == 2
        @inbounds for i = 1:nPoints
            temp = split( text[ pointLine + i ] , r"[\t]" )
            points[ i ] = Point2D( parse( Float64 , temp[ 1 ] ) ,
                                   parse( Float64 , temp[ 2 ] ) )
            # Consider broadcasting (to be checked)
            #points[ i ] = Point2D( parse.( Float64 , temp ) )
        end
    else
        @inbounds for i = 1:nPoints
            temp = split( text[ pointLine + i ] , r"[\t]" )
            points[ i ] = Point3D( parse( Float64 , temp[ 1 ] ) ,
                                   parse( Float64 , temp[ 2 ] ) ,
                                   parse( Float64 , temp[ 3 ] ) )
            # Consider broadcasting (to be checked)
            #points[ i ] = Point3D( parse.( Float64 , temp ) )
        end
    end

    # Store position of the element block
    elementLine = findfirst( s -> occursin( "NELEM=" , s ) , text )
    nElements = parse( Int64 , split( text[ elementLine ] , r"[=]" )[ end ] )

    elements = Vector{Element}( undef , nElements )

    # Retrieve data from the element block and store in a vector
    @inbounds for i = 1:nElements
        temp = split( text[ elementLine + i ] , r"[\t]" )
        id = parse( Int64 , temp[ 1 ] )
        # The element index is found at the end of the line in some legacy mesh files
        connectivity = parse.( Int64 , temp[ 2:numberOfNodes( id )+1 ] )
        elements[ i ] = Element( id ,
                                 connectivity ,
                                 zeros( Int64 , numberOfBoundaries( id , nDimensions , debug ) ) )
    end

    # Store position of the marker block and respective data
    markerLine = findfirst( s -> occursin( "NMARK=" , s ) , text )
    nMarkers = parse( Int64 , split( text[ markerLine ] , r"[=]" )[ end ] )

    markers = Vector{Marker}( undef , nMarkers )

    # Retrieve data from the marker block and store in a vector
    @inbounds for m ∈ 1:nMarkers
        if m == 1
            # The index points at the beginning of the first marker sub-block
            global k = markerLine + 1
        end

        tag = lstrip( split( text[ k ] , r"[= ]" )[ end ] )
        nMarkerElements = parse( Int64 , split( text[ k + 1 ] , r"[=]" )[ end ] )

        markerElements = Vector{MarkerElement}( undef , nMarkerElements )

        @inbounds for j = 1:nMarkerElements
            # Fist number is element type, follow the nodes identifying the element
            temp = split( text[ k + 1 + j  ] , r"[\t]" )
            id = parse( Int64 , temp[ 1 ] )
            connectivity = parse.( Int64 , temp[ 2:end ] )
            markerElements[ j ] = MarkerElement( id , connectivity )
        end

        markers[ m ] = Marker( tag , nMarkerElements , markerElements )

        # The loop must advance to the next marker
        # N.B.: Better solution required, otherwise no multi-threaded outer loop
        k += 2 + nMarkerElements
    end

    return( Mesh( nDimensions , nElements , nPoints , nMarkers , elements , points , markers ) )
end

function assignMeshBoundaries!( mesh::Mesh , debug::Bool=false )
    nDimensions = mesh.nDimensions

    # Sorting avoids the need to compare all possible permutations in the following
    boundaries = [ findStructuredBoundaries( e , debug ) for e ∈ mesh.elements ]
    boundaries = [ [ sort( findPermutations( boundaries[ i ][ j ] ) )[ 1 ]
                     for j ∈ 1:length( boundaries[ 1 ] ) ]
                   for i ∈ 1:length( boundaries ) ]

    flippedBoundaries = [ [ sort( findPermutations( reverse( boundaries[ i ][ j ] ) ) )[ 1 ]
                            for j ∈ 1:length( boundaries[ 1 ] ) ]
                          for i ∈ 1:length( boundaries ) ]

    if debug
        println( "DEBUG: assignMeshBoundaries! -> boundaries[ 1:4 ] " , boundaries[ 1:4 ] )
        println( "DEBUG: assignMeshBoundaries! -> flippedBoundaries[ 1:4 ] " , flippedBoundaries[ 1:4 ] )
    end

    println( "INFO: Assigning boundary markers to boundary elements" )

    # Storing in advance indexes of the elements having external boundaries is expected to save
    # computational time ONLY?/ESPECIALLY when the number of subdivisions is greater for the parent
    # element than for its boundary elements.
    # For instance, a speedup will be surely experienced if working with hexahedra. At each
    # refinement level, 8 times the parent element and 5 times the boundary element will be
    # computed. So, for an increasing number of mesh elements, the saving in terms of time is
    # expected to be increasingly greater.
    elementIndexes = findBoundaryElements( mesh , debug )

    @inbounds for i ∈ 1:mesh.nMarkers
        for l ∈ 1:length( mesh.markers[ i ].elements )
            test = sort( findPermutations( mesh.markers[ i ].elements[ l ].connectivity ) )[ 1 ]
            Threads.@threads for j ∈ elementIndexes
                nBoundaries = length( boundaries[ j ] )
                for k ∈ 1:nBoundaries
                    # Always returns 0 if the edge/face does not belong to the boundary
                    # Returns the boundary index otherwise, with sig dependent on the normal orientation
                    if test == boundaries[ j ][ k ]
                        mesh.elements[ j ].boundary[ k ] = i
                    elseif test == flippedBoundaries[ j ][ k ]
                        mesh.elements[ j ].boundary[ k ] = -i
                        if debug
                            println( "DEBUG: assignMeshboundaries -> FOUND FLIPPED NORMAL" )
                        end
                    end

                    if debug
                        println( "DEBUG: assignMeshboundaries -> boundaries " ,
                                 boundaries[ j ][ k ] )
                        println( "DEBUG: assignMeshboundaries -> marker element connectivity " ,
                                 mesh.markers[ i ].elements[ l ].connectivity )
                        println( "DEBUG: assignMeshboundaries -> check " , check * i )
                    end
                end
            end
        end
    end
end

function outDir( )
    dir = abspath( pwd( ) , "out" )
    if !isdir( dir )
        mkdir( dir )
    end
    return dir
end

function writeMesh( mesh::Mesh , file::String , append::Bool=false , debug::Bool=false )
    path = abspath( file )

    println( "INFO: Writing output to " , path )

    if append
        fileOption = "a"
    else
        fileOption = "w"
    end

    f = open( path , fileOption )

    println( f , "NDIME= " , mesh.nDimensions )

    println( f , "NPOIN= " , mesh.nPoints )

    if mesh.nDimensions == 2
        @inbounds for i ∈ 1:mesh.nPoints
            @printf( f , "%.14f\t%.14f\t%d\n" ,
                     mesh.points[ i ].x ,
                     mesh.points[ i ].y ,
                     i - 1 ) # The index is retained for backwark compatibility (PointWise)
        end
    else
        @inbounds for i ∈ 1:mesh.nPoints
            @printf( f , "%.14f\t%.14f\t%.14f\t%d\n" ,
                     mesh.points[ i ].x ,
                     mesh.points[ i ].y ,
                     mesh.points[ i ].z ,
                     i - 1 ) # The index is retained for backwark compatibility (PointWise)
        end
    end

    println( f , "NELEM= " , mesh.nElements )

    @inbounds for i ∈ 1:mesh.nElements
        @printf( f , "%d" , mesh.elements[ i ].id )

        @inbounds for j ∈ 1:length( mesh.elements[ i ].connectivity )
            @printf( f , "\t%d" , mesh.elements[ i ].connectivity[ j ] )
        end

        @printf( f , "\n" )
    end

    println( f , "NMARK= " , mesh.nMarkers )

    @inbounds for i ∈ 1:mesh.nMarkers
        println( f , "MARKER_TAG= " , mesh.markers[ i ].tag )
        println( f , "MARKER_ELEMS= " , mesh.markers[ i ].nElements )

        @inbounds for j ∈ 1:length( mesh.markers[ i ].elements )
            @printf( f , "%d" , mesh.markers[ i ].elements[ j ].id )
            @inbounds for k ∈ 1:length( mesh.markers[ i ].elements[ j ].connectivity )
                @printf( f , "\t%d" , mesh.markers[ i ].elements[ j ].connectivity[ k ] )
            end
            @printf( f , "\n" )
        end
    end

    close( f )
end

function findQuadEdges( element::Element , debug::Bool=false )
    nEdges = numberOfEdges( element.id , debug )

    # Vector of vectors instead of matrix for direct comparison of outer elements
    edges = Vector{Vector{Int64}}( undef , nEdges )

    @inbounds @simd for i ∈ 1:nEdges-1
        edges[ i ] = [ element.connectivity[ i ] ,
                       element.connectivity[ i + 1 ] ]
    end
    edges[ end ] = [ element.connectivity[ end ] , element.connectivity[ 1 ] ]

    if debug
        println( "DEBUG: findQuadEdges -> edges " , edges )
    end

    return edges
end

function findHexaFaces( element::Element , debug::Bool=false )
    c = element.connectivity

    faces = [ [ c[ 1 ] , c[ 4 ] , c[ 3 ] , c[ 2 ] ] , # XY-
              [ c[ 5 ] , c[ 6 ] , c[ 7 ] , c[ 8 ] ] , # XY+
              [ c[ 2 ] , c[ 3 ] , c[ 7 ] , c[ 6 ] ] , # XZ-
              [ c[ 1 ] , c[ 5 ] , c[ 8 ] , c[ 4 ] ] , # XZ+
              [ c[ 1 ] , c[ 2 ] , c[ 6 ] , c[ 5 ] ] , # YZ-
              [ c[ 3 ] , c[ 4 ] , c[ 8 ] , c[ 7 ] ] ] # YZ+

    if debug
        println( "DEBUG: findHexaFaces -> faces " , faces )
    end

    return faces
end

function findStructuredBoundaries( element::Element , debug::Bool=false )
    if element.id == 9
        return findQuadEdges( element , debug )
    elseif element.id == 12
        return findHexaFaces( element , debug )
    end
end

function isNeighbor( mesh::Mesh , elementIndex1::Int64 , elementIndex2::Int64 , debug::Bool=false )
    if mesh.nDimensions == 2
        if debug
            println( "DEBUG: isNeighbor -> findEdge( mesh , elementIndex1 ) " ,
                     findQuadEdges( mesh , elementIndex1 , debug ) )
            println( "DEBUG: isNeighbor -> findEdge( mesh , elementIndex2 ) " ,
                     findQuadEdges( mesh , elementIndex2 , debug ) )
        end

        edges1 = findQuadEdges( mesh , elementIndex1 , debug )
        edges2 = findQuadEdges( mesh , elementIndex2 , debug )
        edges2bis = [ [ edges2[ i ][ 2 ] , edges2[ i ][ 1 ] ] for i ∈ 1:length( edges2 ) ]

        if length( intersect( edges1 , edges2 ) ) + length( intersect( edges1 , edges2bis ) ) > 0
            return true
        else
            return false
        end
    end
end

function findNeighbors( mesh::Mesh , elementIndex::Int64 , debug::Bool=false )
    # Number of neighbor elements not known in advance, same for their indexes
    indexes = Vector{Int64}( )

    @inbounds for i ∈ filter!( x -> x ≠ elementIndex , collect( 1:mesh.nElements ) )
        if isNeighbor( mesh , elementIndex , i , debug ) == true
            push!( indexes , i )
        end
    end

    if debug
        println( "DEBUG: findNeighbors -> indexes " , indexes )
    end

    return indexes
end

function midPoint( points::Vector{Point2D} , debug::Bool=false )
    nPoints = length( points )

    x = [ point.x for point ∈ points ]
    y = [ point.y for point ∈ points ]

    if debug
        println( "DEBUG: centroid -> x " , x )
        println( "DEBUG: centroid -> y " , y )
    end

    return Point2D( sum( x ) / nPoints , sum( y ) / nPoints )
end

function midPoint( points::Vector{Point3D} , debug::Bool=false )
    nPoints = length( points )

    x = [ point.x for point ∈ points ]
    y = [ point.y for point ∈ points ]
    z = [ point.z for point ∈ points ]

    if debug
        println( "DEBUG: centroid -> x " , x )
        println( "DEBUG: centroid -> y " , y )
        println( "DEBUG: centroid -> z " , z )
    end

    return Point3D( sum( x ) / nPoints , sum( y ) / nPoints , sum( z ) / nPoints )
end

function polygonCentroid( p::Vector{Point2D} , debug::Bool=false )
    pp = vcat( p , p[ 1 ] )

    A = 0.5*sum( [ ( pp[ i ].x * pp[ i+1 ].y - pp[ i+1 ].x * pp[ i ].y )
                   for i ∈ 1:length( pp ) - 1 ] )

    Cx = sum( [ ( pp[ i ].x + pp[ i+1 ].x ) * ( pp[ i ].x * pp[ i+1 ].y - pp[ i+1 ].x * pp[ i ].y )
                for i ∈ 1:length( pp )  - 1 ] ) / A / 6

    Cy = sum( [ ( pp[ i ].y + pp[ i+1 ].y ) * ( pp[ i ].x * pp[ i+1 ].y - pp[ i+1 ].x * pp[ i ].y )
                for i ∈ 1:length( pp ) - 1 ] ) / A / 6

    C = Point2D( Cx , Cy )

    return C
end

function polygonCentroid( p::Vector{Point3D} , debug::Bool=false )
    # Group i-th coordinates
    X = [ point.x for point ∈ p ]
    Y = [ point.y for point ∈ p ]
    Z = [ point.z for point ∈ p ]

    # Check if polygon is normal to any axis
    XY = all( c -> c == Z[ 1 ] , Z )
    YZ = all( c -> c == X[ 1 ] , X )
    XZ = all( c -> c == Y[ 1 ] , Y )

    if XY
        C = polygonCentroid( [ Point2D( p[ i ].x , p[ i ].y ) for i in 1:length( p ) ] )
        C = Point3D( C.x , C.y , Z[ 1 ] )
    elseif YZ
        C = polygonCentroid( [ Point2D( p[ i ].y , p[ i ].z ) for i in 1:length( p ) ] )
        C = Point3D( X[ 1 ] , C.x , C.y )
    elseif XZ
        C = polygonCentroid( [ Point2D( p[ i ].x , p[ i ].z ) for i in 1:length( p ) ] )
        C = Point3D( C.x , Y[ 1 ] , C.y )
    else
        pp = vcat( p , p[ 1 ] )

        A = 0.5*sum( [ ( pp[ i ].x * pp[ i+1 ].y - pp[ i+1 ].x * pp[ i ].y )
                       for i ∈ 1:length( pp ) - 1 ] )
        A2 = 0.5*sum( [ ( pp[ i ].x * pp[ i+1 ].z - pp[ i+1 ].x * pp[ i ].z )
                        for i ∈ 1:length( pp ) - 1 ] )
        Cx = sum( [ ( pp[ i ].x + pp[ i+1 ].x ) * ( pp[ i ].x * pp[ i+1 ].y - pp[ i+1 ].x * pp[ i ].y )
                    for i ∈ 1:length( pp )  - 1 ] ) / A / 6
        Cx2 = sum( [ ( pp[ i ].x + pp[ i+1 ].x ) * ( pp[ i ].x * pp[ i+1 ].z - pp[ i+1 ].x * pp[ i ].z )
                     for i ∈ 1:length( pp ) - 1 ] ) / A
        Cy = sum( [ ( pp[ i ].y + pp[ i+1 ].y ) * ( pp[ i ].x * pp[ i+1 ].y - pp[ i+1 ].x * pp[ i ].y )
                    for i ∈ 1:length( pp ) - 1 ] ) / A / 6
        Cz = sum( [ ( pp[ i ].z + pp[ i+1 ].z ) * ( pp[ i ].x * pp[ i+1 ].z - pp[ i+1 ].x * pp[ i ].z )
                    for i ∈ 1:length( pp ) - 1 ] ) / A2 / 6

        C = Point3D( Cx , Cy , Cz )
    end

    return C
end

function centroid2( points::Vector{Point2D} , debug::Bool=false )
    let nPoints = length( points )
        return Point2D( sum( [ point.x for point ∈ points ] ) / nPoints ,
                        sum( [ point.y for point ∈ points ] ) / nPoints )
    end
end

function centroid2( points::Vector{Point3D} , debug::Bool=false )
    let nPoints = length( points )
        return Point3D( sum( [ point.x for point ∈ points ] ) / nPoints ,
                        sum( [ point.y for point ∈ points ] ) / nPoints ,
                        sum( [ point.z for point ∈ points ] ) / nPoints )
    end
end

function hexa2facesIndexes( i::Vector{Int64} , debug::Bool=false )
    return [ [ i[ 1 ] , i[ 2 ] , i[ 3 ] , i[ 4 ] ] , # XZ-
             [ i[ 8 ] , i[ 7 ] , i[ 6 ] , i[ 5 ] ] , # XZ+
             [ i[ 3 ] , i[ 2 ] , i[ 6 ] , i[ 8 ] ] , # XY-
             [ i[ 1 ] , i[ 4 ] , i[ 7 ] , i[ 5 ] ] , # XY+
             [ i[ 4 ] , i[ 3 ] , i[ 8 ] , i[ 7 ] ] , # YZ-
             [ i[ 2 ] , i[ 1 ] , i[ 5 ] , i[ 6 ] ] ] # YZ+
end

function permute( vector , debug::Bool=false )
    if length( vector ) < length( vector[ 1 ] )
        return hcat( [ vector[ end ] ] , vector[ 1:end-1 ] )
    else
        return vcat( [ vector[ end ] ] , vector[ 1:end-1 ] )
    end
end

function findPermutations( vector , debug::Bool=false )
    permutations = Vector{}( undef , length( vector ) )
    permutations[ 1 ] = vector
    @inbounds for i ∈ 2:length( permutations )
        permutations[ i ] = permute( permutations[ i - 1 ] )
    end

    return permutations
end

function nullElements( n , d , debug::Bool=false )
    # N.B.: Only works for quadrilaterals and tetrahedra!
    return [ Element( 0 , zeros( 2 * ( d - 1 ) ) ) for i ∈ 1:n ]
end

function disordered( a , b , debug::Bool=false )
    return sort( a ) == sort( b )
end

function findBoundaryPoints( mesh::Mesh , debug::Bool=false )
    indexes = Vector{Int64}( )

    @inbounds for m ∈ 1:mesh.nMarkers
        for e ∈ 1:length( mesh.markers[ m ].elements )
            for p ∈ mesh.markers[ m ].elements[ e ].connectivity
                push!( indexes , p )
            end
        end
    end

    indexes = sort( unique( indexes ) )

    if debug
        println( "DEBUG: findBoundaryPoints -> indexes " , indexes )
    end

    return indexes
end

function findInnerPoints( mesh::Mesh , debug::Bool=false )
    boundaryPointIndexes = findBoundaryPoints( mesh , debug )

    # Note that points are natively 0-base indexed
    allIndexes = collect( 0:mesh.nPoints-1 )

    indexes = findall( x -> x ∉ boundaryPointIndexes , allIndexes )

    # Note that points are natively 0-base indexed
    indexes .-= 1

    if debug
        println( "DEBUG: findInnerPoints -> indexes " , indexes )
    end

    return indexes
end

function isBoundaryElement( element::MarkerElement , mesh::Mesh , debug::Bool=false )
    boundaryPoints = findBoundaryPoints( mesh , debug )
    return all( x -> x ∈ boundaryPoints , element.connectivity )
end

function findBoundaryElements( mesh::Mesh , debug::Bool=false )
    indexes = Vector{Int64}( )
    points = findBoundaryPoints( mesh , debug )

    @inbounds for i ∈ 1:mesh.nElements
        nPoints = numberOfPointsOnBoundaries( mesh.elements[ i ].id , debug )
        if length( findall( x -> x ∈ points , mesh.elements[ i ].connectivity ) ) ≥
            numberOfPointsOnBoundaries( mesh.elements[ i ].id )
            push!( indexes , i )
        end
    end

    return indexes
end

function refineQuad( element::Element , points::Vector{Point2D} ,
                     newElements::Vector{Element} , newPoints::Vector{Point2D} , newMarkers::Vector{Marker} ,
                     eCounter::Int64 , mCounter::Vector{Int64} , debug::Bool=false )
    # G --- H --- I
    # |     |     |
    # D --- E --- F
    # |     |     |
    # A --- B --- C

    ( A , C , I , G ) = ( points[ 1 ] , points[ 2 ] , points[ 3 ] , points[ 4 ] )

    # Compute old element centroid
    #E = polygonCentroid( [ A , C , I , G ] , debug)
    E = midPoint( [ midPoint( [ A , I ] ) , midPoint( [ C , G ] ) ] )

    # Compute edge mid-points
    ( B , F , H , D ) = midPoint.( [ [ A , C ] , [ C , I ] , [ G , I ] , [ G , A ] ] , debug )

    # Group all element points under a common name
    newPointSet = [ A , B , C , D , E , F , G , H , I ]

    if debug
        println( "DEBUG: refineQuad -> newPointSet " , newPointSet )
    end

    # Assign index to all element points
    newIndexes = zeros( Int64 , 9 )
    @inbounds for i ∈ 1:9
        temp = findfirst( x -> x == newPointSet[ i ] , newPoints )

        if typeof( temp ) == Nothing
            # New mesh point!
            newIndexes[ i ] = length( newPoints )

            # Store new mesh point
            push!( newPoints , newPointSet[ i ] )
        else
            newIndexes[ i ] = temp - 1
        end
    end

    if debug
        println( "DEBUG: refineQuad -> newIndexes " , newIndexes )
    end

    idxA = newIndexes[ 1 ] ; idxB = newIndexes[ 2 ] ; idxC = newIndexes[ 3 ]
    idxD = newIndexes[ 4 ] ; idxE = newIndexes[ 5 ] ; idxF = newIndexes[ 6 ]
    idxG = newIndexes[ 7 ] ; idxH = newIndexes[ 8 ] ; idxI = newIndexes[ 9 ]

    # Instantiate new elements
    quad = ( Element( 9 , [ idxA , idxB , idxE , idxD ] , zeros( Int64 , 4 ) ) ,
             Element( 9 , [ idxB , idxC , idxF , idxE ] , zeros( Int64 , 4 ) ) ,
             Element( 9 , [ idxD , idxE , idxH , idxG ] , zeros( Int64 , 4 ) ) ,
             Element( 9 , [ idxE , idxF , idxI , idxH ] , zeros( Int64 , 4 ) ) )

    if debug
        println( "DEBUG: refineQuad -> quad " , quad )
        println( "DEBUG: refineQuad -> eCounter " , eCounter )
    end

    # Store new mesh elements
    @inbounds @simd for i ∈ 1:4
        newElements[ eCounter ] = quad[ i ]
        eCounter += 1
    end

    if debug
        println( "DEBUG: refineQuad -> newElements " , newElements )
    end

    positiveIndexes = [ [ [ idxA , idxB ] , [ idxB , idxC ] ] ,
                        [ [ idxC , idxF ] , [ idxF , idxI ] ] ,
                        [ [ idxG , idxH ] , [ idxH , idxI ] ] ,
                        [ [ idxA , idxD ] , [ idxD , idxG ] ] ]

    negativeIndexes = [ [ reverse( positiveIndexes[ i ][ j ] ) for j ∈ 1:2 ] for i ∈ 1:4 ]

    # Transfer marker (boundary) information
    @inbounds @simd for i ∈ 1:4
        if element.boundary[ i ] > 0
            markerIndex = element.boundary[ i ]
            k = mCounter[ markerIndex ]
            newMarkers[ markerIndex ].elements[ k ].connectivity = positiveIndexes[ i ][ 1 ]
            newMarkers[ markerIndex ].elements[ k + 1 ].connectivity = positiveIndexes[ i ][ 2 ]
            mCounter[ markerIndex ] += 2
        elseif element.boundary[ i ] < 0
            markerIndex = - element.boundary[ i ]
            k = mCounter[ markerIndex ]
            newMarkers[ markerIndex ].elements[ k ].connectivity = negativeIndexes[ i ][ 1 ]
            newMarkers[ markerIndex ].elements[ k + 1 ].connectivity = negativeIndexes[ i ][ 2 ]
            mCounter[ markerIndex ] += 2
        end
    end

    return eCounter , mCounter
end

function refineHexa( element::Element , points::Vector{Point3D} ,
                     newElements::Vector{Element} , newPoints::Vector{Point3D} , newMarkers::Vector{Marker} ,
                     eCounter::Int64 , mCounter::Vector{Int64} , debug::Bool=false )
    #                            G3 --- H3 --- I3
    #                           /|      |     /|
    #                          / D3 --- E3 --/ F3
    #                         /  |      |   /  |
    #                        /   A3 --- B3 /-- C3
    #                       /             /    /
    #                      G2 --- H2 --- I2   /
    #                     /|      |     /|   /
    #                    / D2 --- E2 --/ F2 /
    #                   /  |      |   /  | /
    #                  /   A2 --- B2 /-- C2
    #                 /             /    /
    #                 G1 --- H1 --- I1  /
    #   Z             |      |      |  /
    #   |  Y          D1 --- E1 --- F1/
    #   | /           |      |      |/
    #   |/_ _ X       A1 --- B1 --- C1

    ( A3 , A1 , C1 , C3 ) = ( points[ 1 ] , points[ 2 ] , points[ 3 ] , points[ 4 ] )
    ( G3 , G1 , I1 , I3 ) = ( points[ 5 ] , points[ 6 ] , points[ 7 ] , points[ 8 ] )

    # Compute edge mid-points
    ( A2 , C2 , I2 , G2 ) = midPoint.( [ [ A1 , A3 ] , [ C1 , C3 ] ,
                                         [ I1 , I3 ] , [ G1 , G3 ] ] , debug )

    ( B1 , F1 , H1 , D1 ) = midPoint.( [ [ A1 , C1 ] , [ C1 , I1 ] ,
                                         [ G1 , I1 ] , [ G1 , A1 ] ] , debug )

    ( B3 , F3 , H3 , D3 ) = midPoint.( [ [ A3 , C3 ] , [ C3 , I3 ] ,
                                         [ G3 , I3 ] , [ G3 , A3 ] ] , debug )

    #( B2 , F2 , H2 , D2 ) = midPoint.( [ [ B1 , B3 ] , [ F1 , F3 ] ,
    #                                     [ H1 , H3 ] , [ D1 , D3 ] ] , debug )

    B2 = midPoint( [ midPoint( [ A1 , C3 ] ) , midPoint( [ C1 , A3 ] ) ] )
    H2 = midPoint( [ midPoint( [ G1 , I3 ] ) , midPoint( [ I1 , G3 ] ) ] )

    E1 = midPoint( [ midPoint( [ A1 , I1 ] ) , midPoint( [ C1 , G1 ] ) ] )
    E3 = midPoint( [ midPoint( [ A3 , I3 ] ) , midPoint( [ C3 , G3 ] ) ] )

    D2 = midPoint( [ midPoint( [ A1 , G3 ] ) , midPoint( [ G1 , A3 ] ) ] )
    F2 = midPoint( [ midPoint( [ C1 , I3 ] ) , midPoint( [ I1 , C3 ] ) ] )

    E2 = midPoint( [ midPoint( [ C1 , G3 ] ) ,
                     midPoint( [ A1 , I3 ] ) ,
                     midPoint( [ I1 , A3 ] ) ,
                     midPoint( [ G1 , C3 ] ) ] )

    # Compute face centroids
    #E2 = centroid2( [ A1 , C1 , I1 , G1 , A3 , C3 , I3 , G3 ] , debug )
    #( E1 , E2 , E3 ) = polygonCentroid.( [ [ A1 , C1 , I1 , G1 ] ,
    #                                       [ A2 , C2 , I2 , G2 ] ,
    #                                       [ A3 , C3 , I3 , G3 ] ] , debug )
    #
    #( B2 , H2 , D2 , F2 ) = polygonCentroid.( [ [ A1 , A3 , C3 , C1 ] ,
    #                                            [ G1 , I1 , I3 , G3 ] ,
    #                                            [ A1 , A3 , G3 , G1 ] ,
    #                                            [ C1 , C3 , I3 , I1 ] ] , debug )

    # Group all element points under a common name
    newPointSet = [ A1 , B1 , C1 , D1 , E1 , F1 , G1 , H1 , I1 ,
                    A2 , B2 , C2 , D2 , E2 , F2 , G2 , H2 , I2 ,
                    A3 , B3 , C3 , D3 , E3 , F3 , G3 , H3 , I3 ]

    if debug
        println( "DEBUG: refineHexa -> newPointSet " , newPointSet )
    end

    # Assign index to all element points
    newIndexes = zeros( Int64 , 27 )
    @inbounds for i ∈ 1:27
        temp = findfirst( x -> x == newPointSet[ i ] , newPoints )

        if typeof( temp ) == Nothing
            # New mesh point!
            newIndexes[ i ] = length( newPoints )

            # Store new mesh point
            push!( newPoints , newPointSet[ i ] )
        else
            newIndexes[ i ] = temp - 1
        end
    end

    if debug
        println( "DEBUG: refineHexa -> newIndexes " , newIndexes )
    end

    idxA1 = newIndexes[ 1 ] ; idxB1 = newIndexes[ 2 ] ; idxC1 = newIndexes[ 3 ]
    idxD1 = newIndexes[ 4 ] ; idxE1 = newIndexes[ 5 ] ; idxF1 = newIndexes[ 6 ]
    idxG1 = newIndexes[ 7 ] ; idxH1 = newIndexes[ 8 ] ; idxI1 = newIndexes[ 9 ]

    idxA2 = newIndexes[ 10 ] ; idxB2 = newIndexes[ 11 ] ; idxC2 = newIndexes[ 12 ]
    idxD2 = newIndexes[ 13 ] ; idxE2 = newIndexes[ 14 ] ; idxF2 = newIndexes[ 15 ]
    idxG2 = newIndexes[ 16 ] ; idxH2 = newIndexes[ 17 ] ; idxI2 = newIndexes[ 18 ]

    idxA3 = newIndexes[ 19 ] ; idxB3 = newIndexes[ 20 ] ; idxC3 = newIndexes[ 21 ]
    idxD3 = newIndexes[ 22 ] ; idxE3 = newIndexes[ 23 ] ; idxF3 = newIndexes[ 24 ]
    idxG3 = newIndexes[ 25 ] ; idxH3 = newIndexes[ 26 ] ; idxI3 = newIndexes[ 27 ]

    # Instantiate and store new elements
    newElements[ eCounter ] =
        Element( 12 , [ idxA2 , idxA1 , idxB1 , idxB2 , idxD2 , idxD1 , idxE1 , idxE2 ] ,
                 zeros( Int64 , 6 ) )
    newElements[ eCounter + 1 ] =
        Element( 12 , [ idxB2 , idxB1 , idxC1 , idxC2 , idxE2 , idxE1 , idxF1 , idxF2 ] ,
                 zeros( Int64 , 6 ) )
    newElements[ eCounter + 2 ] =
        Element( 12 , [ idxD2 , idxD1 , idxE1 , idxE2 , idxG2 , idxG1 , idxH1 , idxH2 ] ,
                 zeros( Int64 , 6 ) )
    newElements[ eCounter + 3 ] =
        Element( 12 , [ idxE2 , idxE1 , idxF1 , idxF2 , idxH2 , idxH1 , idxI1 , idxI2 ] ,
                 zeros( Int64 , 6 ) )
    newElements[ eCounter + 4 ] =
        Element( 12 , [ idxA3 , idxA2 , idxB2 , idxB3 , idxD3 , idxD2 , idxE2 , idxE3 ] ,
                 zeros( Int64 , 6 ) )
    newElements[ eCounter + 5 ] =
        Element( 12 , [ idxB3 , idxB2 , idxC2 , idxC3 , idxE3 , idxE2 , idxF2 , idxF3 ] ,
                 zeros( Int64 , 6 ) )
    newElements[ eCounter + 6 ] =
        Element( 12 , [ idxD3 , idxD2 , idxE2 , idxE3 , idxG3 , idxG2 , idxH2 , idxH3 ] ,
                 zeros( Int64 , 6 ) )
    newElements[ eCounter + 7 ] =
        Element( 12 , [ idxE3 , idxE2 , idxF2 , idxF3 , idxH3 , idxH2 , idxI2 , idxI3 ] ,
                 zeros( Int64 , 6 ) )

    eCounter += 8

    if debug
        println( "DEBUG: refineHexa -> eCounter ", eCounter )
        println( "DEBUG: refineHexa -> newElements ", newElements )
    end

    positiveIndexes = [ [ [ idxA1 , idxA2 , idxB2 , idxB1 ] , [ idxB1 , idxB2 , idxC2 , idxC1 ] ,
                          [ idxA2 , idxA3 , idxB3 , idxB2 ] , [ idxB2 , idxB3 , idxC3 , idxC2 ] ] ,
                        [ [ idxG2 , idxG1 , idxH1 , idxH2 ] , [ idxH2 , idxH1 , idxI1 , idxI2 ] ,
                          [ idxG3 , idxG2 , idxH2 , idxH3 ] , [ idxH3 , idxH2 , idxI2 , idxI3 ] ] ,
                        [ [ idxA1 , idxB1 , idxE1 , idxD1 ] , [ idxB1 , idxC1 , idxF1 , idxE1 ] ,
                          [ idxD1 , idxE1 , idxH1 , idxG1 ] , [ idxE1 , idxF1 , idxI1 , idxH1 ] ] ,
                        [ [ idxA3 , idxD3 , idxE3 , idxB3 ] , [ idxB3 , idxE3 , idxF3 , idxC3 ] ,
                          [ idxD3 , idxG3 , idxH3 , idxE3 ] , [ idxE3 , idxH3 , idxI3 , idxF3 ] ] ,
                        [ [ idxA2 , idxA1 , idxD1 , idxD2 ] , [ idxD2 , idxD1 , idxG1 , idxG2 ] ,
                          [ idxA3 , idxA2 , idxD2 , idxD3 ] , [ idxD3 , idxD2 , idxG2 , idxG3 ] ] ,
                        [ [ idxC1 , idxC2 , idxF2 , idxF1 ] , [ idxF1 , idxF2 , idxI2 , idxI1 ] ,
                          [ idxC2 , idxC3 , idxF3 , idxF2 ] , [ idxF2 , idxF3 , idxI3 , idxI2 ] ] ]

    negativeIndexes = [ [ reverse( positiveIndexes[ i ][ j ] ) for j ∈ 1:4 ] for i ∈ 1:6 ]

    # Transfer marker (boundary) information
    @inbounds for i ∈ 1:6
        if debug
            println( "DEBUG: refineHexa -> element.boundary[ i ] " , element.boundary[ i ] )
        end
        if element.boundary[ i ] > 0
            markerIndex = element.boundary[ i ]
            k = mCounter[ markerIndex ]
            newMarkers[ markerIndex ].elements[ k ].connectivity = positiveIndexes[ i ][ 1 ]
            newMarkers[ markerIndex ].elements[ k + 1 ].connectivity = positiveIndexes[ i ][ 2 ]
            newMarkers[ markerIndex ].elements[ k + 2 ].connectivity = positiveIndexes[ i ][ 3 ]
            newMarkers[ markerIndex ].elements[ k + 3 ].connectivity = positiveIndexes[ i ][ 4 ]
            mCounter[ markerIndex ] += 4
        elseif element.boundary[ i ] < 0
            markerIndex =  - element.boundary[ i ]
            k = mCounter[ markerIndex ]
            newMarkers[ markerIndex ].elements[ k ].connectivity = negativeIndexes[ i ][ 1 ]
            newMarkers[ markerIndex ].elements[ k + 1 ].connectivity = negativeIndexes[ i ][ 2 ]
            newMarkers[ markerIndex ].elements[ k + 2 ].connectivity = negativeIndexes[ i ][ 3 ]
            newMarkers[ markerIndex ].elements[ k + 3 ].connectivity = negativeIndexes[ i ][ 4 ]
            mCounter[ markerIndex ] += 4
        end
    end

    return eCounter , mCounter
end

function refineStructured( element::Element , points::Vector{Point2D} ,
                           newElements::Vector{Element} , newPoints::Vector{Point2D} , newMarkers::Vector{Marker} ,
                           eCounter::Int64 , mCounter::Vector{Int64} , debug::Bool=false )
    return refineQuad( element , points , newElements , newPoints , newMarkers , eCounter , mCounter, debug )
end

function refineStructured( element::Element , points::Vector{Point3D} ,
                           newElements::Vector{Element} , newPoints::Vector{Point3D} , newMarkers::Vector{Marker} ,
                           eCounter::Int64 , mCounter::Vector{Int64} , debug::Bool=false )
    return refineHexa( element , points , newElements , newPoints , newMarkers , eCounter , mCounter, debug )
end

function refineStructuredMesh( mesh::Mesh , refinementLevel=1 , debug::Bool=false )
    # Assign marker indexes to elements before any further operation
    assignMeshBoundaries!( mesh , debug )

    # Copy marker structured data (mutable)
    newMarkers = deepcopy( mesh.markers )

    if mesh.nDimensions == 2
        newPoints = Vector{Point2D}( )
        newElements = Vector{Element}( undef , mesh.nElements * 4 )

        @inbounds for i ∈ 1:mesh.nMarkers
            newMarkers[ i ].nElements = mesh.markers[ i ].nElements * 2
            baseN = mesh.markers[ i ].nElements
            # Every edge will be divided into 2
            resize!( newMarkers[ i ].elements , baseN * 2 )
            # Initialize new marker elements
            @inbounds for j ∈ ( baseN + 1 ):baseN * 2
                newMarkers[ i ].elements[ j ] = MarkerElement( 3 , zeros( Int64 , 2 ) )
            end
        end
    else
        newPoints = Vector{Point3D}( )
        newElements = Vector{Element}( undef , mesh.nElements * 8 )

        @inbounds for i ∈ 1:mesh.nMarkers
            newMarkers[ i ].nElements = mesh.markers[ i ].nElements * 4
            baseN = mesh.markers[ i ].nElements
            # Every face will be divided into 4
            resize!( newMarkers[ i ].elements , baseN * 4 )
            # Initialize new marker elements
            @inbounds for j ∈ ( baseN + 1 ):baseN * 4
                newMarkers[ i ].elements[ j ] = MarkerElement( 9 , zeros( Int64 , 4 ) )
            end
        end
    end

    # Support counters for helping array preallocation strategy
    eCounter = 1
    mCounter = ones( Int64 , mesh.nMarkers )

    println( "INFO: Subdividing elements" )

    @inbounds for i ∈ 1:mesh.nElements
        element = mesh.elements[ i ]
        points = [ mesh.points[ i + 1 ] for i ∈ element.connectivity ]

        eCounter , mCounter = refineStructured( element , points ,
                                                newElements , newPoints , newMarkers ,
                                                eCounter , mCounter , debug )

        if debug
            println( "DEBUG: refineStructuredMesh -> newPoints " , newPoints )
            println( "DEBUG: refineStructuredMesh -> newElements " , newElements )
            println( "DEBUG: refineStructuredMesh -> newMarkers " , newMarkers )
        end
    end

    newMesh = Mesh( mesh.nDimensions ,
                    length( newElements ) , length( newPoints ) , length( newMarkers ) ,
                    newElements , newPoints , newMarkers )

    if refinementLevel == 1
        return newMesh
    else
        return refineStructuredMesh( newMesh , refinementLevel - 1 , debug )
    end
end

function couples( series::Vector{Int64} , debug::Bool=false )
    augmented = vcat( series , series[ 1 ] )
    new = [ [ augmented[ i ] , augmented[ i + 1 ] ] for i ∈ 1:length( augmented ) - 1 ]
    return new
end

function buildNetwork2D( mesh::Mesh , debug::Bool=false )
    associations = Vector{Vector{Int64}}( )

    for e ∈ mesh.elements
        for c ∈ couples( e.connectivity .+ 1 )
            push!( associations , c )
        end
    end

    sort!.( associations )
    sort!( associations )
    unique!( associations )

    if debug
        println( "DEBUG: buildNetwork2D -> associations " , associations )
    end

    network = Graph( length( associations ) )
    for a ∈ associations
        add_edge!( network , a[ 1 ]  , a[ 2 ] )
    end

    return network
end

function laplacianSmoothing2D( mesh::Mesh , nIterations::Int64 , debug::Bool=false )
    newMesh = deepcopy( mesh )

    # Build the network based on mesh points
    network = buildNetwork2D( newMesh , debug )

    # Only the inner points must be (eventually) moved
    innerPointIndexes = findInnerPoints( newMesh , debug )

    # Pass to 1-based indexing
    innerPointIndexes .+= 1

    # Update coordinates
    for k ∈ 1:nIterations
        for i ∈ innerPointIndexes
            point = newMesh.points[ i ]
            pointCoordinates = [ point.x , point.y ]

            neighborIndexes = neighbors( network , i )
            neighborPoints = [ newMesh.points[ j ] for j ∈ neighborIndexes ]

            N = length( neighborPoints )

            newMesh.points[ i ].x = sum( [ p.x for p in neighborPoints ] ) / N
            newMesh.points[ i ].y = sum( [ p.y for p in neighborPoints ] ) / N
        end
    end

    return newMesh
end

function desplit( str::Array{SubString{String},1} , debug::Bool=false )
    s = str[ 1 ]
    for i ∈ 2:length( str )
        s *= "/"
        s *= str[ i ]
    end
    return s
end

function mergeZones( zones::Vector{Mesh} , debug::Bool=false )
    nZones = length( zones )
    out = outDir( )
    fileName = "merged.su2"
    path = abspath( out , fileName )

    f = open( path , "w" )
    println( f , "NZONE= " , nZones )
    close( f )

    for i ∈ 1:nZones
        f = open( path , "a" )
        println( f )
        println( f , "IZONE= " , i )
        close( f )
        writeMesh( zones[ i ] , fileName , true , debug )
    end
end

end # module
