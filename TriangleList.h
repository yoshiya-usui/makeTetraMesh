//--------------------------------------------------------------------------
// This file is part of makeTetraMesh.
//
// makeTetraMesh is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// makeTetraMesh is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with makeTetraMesh. If not, see <http://www.gnu.org/licenses/>.
//--------------------------------------------------------------------------
#ifndef DBLDEF_TRIANGLE_LIST
#define DBLDEF_TRIANGLE_LIST

#include "Triangle.h"
#include "NodeList.h"
#include "CommonParameters.h"
#include "BoundaryCurveList.h"
#include "BoundaryCurve.h"
#include "TopographyData.h"
#include "AnalysisDomain.h"

#include <vector>
#include <map>

// Class of the list of triangles
class TriangleList{

public:

	// Constructer
	explicit TriangleList();

	// Destructer
	~TriangleList();

	// Create surface triangles
	void createSurfaceTriangles();

	// Create surface triangles on additional plane
	void createSurfaceTrianglesOnAdditionalPlane();

	// Write triangles without the ones belonging to supert triangle to intermediate file after step 2
	void writeTrianglesToIntermediateFileStep2( const std::string& fileName ) const;

	// Write triangles without the ones belonging to supert triangle to intermediate file after step 3
	void writeTrianglesToIntermediateFileStep3( const std::string& fileName ) const;

	// Write nodes of triangles without the ones of supert triangle to intermediate file after step 2
	void writeNode2DListToIntermediateFileStep2( const std::string& fileName ) const;

	// Write nodes of triangles without the ones of supert triangle to intermediate file after step 3
	void writeNodeListToIntermediateFileStep3( const std::string& fileName ) const;

	// Read triangles without the ones belonging to supert triangle from intermediate file
	void readTrianglesFromIntermediateFile( const std::string& fileName );

	// Read nodes of triangles without the ones of supert triangle from intermediate file
	void readNode2DListFromIntermediateFile( const std::string& fileName );

	// Read nodes of triangles without the ones of supert triangle from intermediate file after step 3
	void readNodeListFromIntermediateFileStep3( const std::string& fileName );

	// Relate nodes to triangles
	void relateNodesToTriangles();

	// Assign location flag to nodes
	void assignLocationToNodes();

	// Interpolate altitude to nodes
	void interpolateAltitudeToNodes();

	// Write array convert node ID of boundary curve list to the one of triangle list 
	void writeNodeIDBoundCurve2Triangle( const std::string& fileName ) const;

	// Write array convert node ID of boundary curve list to the one of triangle list after step3
	void writeNodeIDBoundCurve2TriangleStep3( const std::string& fileName ) const;

	// Read array convert node ID of boundary curve list to the one of triangle list 
	void readNodeIDBoundCurve2Triangle( const std::string& fileName );

	// Write PLC for TETGEN
	void writePLCs();

	// Get pointer to the class of boundary curve list
	BoundaryCurveList* getPointerToBoundaryCurveList();

	// Set data of one outer boundary curve
	void setDataOfOneOuterBoundaryCurve( const int numNodes, const int* nodeIDs, const CommonParameters::XY* coord2D, const int geologicalType );

	// Set location of the plane on which 2D mesh is created
	void setPlaneLocation( const CommonParameters::Boundary& planeLocation, const int iLayer = -1 );

private:

	//typedef std::vector< std::pair< CommonParameters::XY , std::pair<int,int> > > ArrayType;
	typedef std::pair<int,int> NodePair;

	//typedef std::pair<int,NodePair> ElemAndNodePair;
	typedef std::map< double, NodePair > ArrayType;

	typedef std::vector<int> PolygonNodes;

	struct Facet{
		std::vector<PolygonNodes> polygons;
		std::vector<CommonParameters::XYZ> holes;
	};

	enum CornerPoint{
		UNDEFINED_CORNER = -1,
		XPLUS_YMINUS_TOP = 0,
		XPLUS_YPLUS_TOP,
		XMINUS_YPLUS_TOP,
		XMINUS_YMINUS_TOP,
		XPLUS_YMINUS_SURFACE,
		XPLUS_YPLUS_SURFACE,
		XMINUS_YPLUS_SURFACE,
		XMINUS_YMINUS_SURFACE,
		XPLUS_YMINUS_BOTTOM,
		XPLUS_YPLUS_BOTTOM,
		XMINUS_YPLUS_BOTTOM,
		XMINUS_YMINUS_BOTTOM,
	};

	enum Side{
		UNDEFINED_SIDE = -1,
		XMINUS_TOP = 0,
		YMINUS_TOP,
		XPLUS_TOP,
		YPLUS_TOP,
		XMINUS_BOTTOM,
		YMINUS_BOTTOM,
		XPLUS_BOTTOM,
		YPLUS_BOTTOM,
		XPLUS_YMINUS_VERTICAL_AIR,
		XPLUS_YPLUS_VERTICAL_AIR,
		XMINUS_YPLUS_VERTICAL_AIR,
		XMINUS_YMINUS_VERTICAL_AIR,
		XPLUS_YMINUS_VERTICAL_SEA,
		XPLUS_YPLUS_VERTICAL_SEA,
		XMINUS_YPLUS_VERTICAL_SEA,
		XMINUS_YMINUS_VERTICAL_SEA,
		XPLUS_YMINUS_VERTICAL_LAND,
		XPLUS_YPLUS_VERTICAL_LAND,
		XMINUS_YPLUS_VERTICAL_LAND,
		XMINUS_YMINUS_VERTICAL_LAND,
	};

	// Copy constructer
	TriangleList(const TriangleList& rhs);

	// Assignment operator
	TriangleList& operator=(const TriangleList& rhs);

	std::vector< Triangle > m_triangles;

	int m_lastSearchedTriangle;

	// List of nodes
	NodeList m_nodeList;

	// List of boundary curve
	BoundaryCurveList m_boundaryCurveList;

	// Location of the plane on which 2D mesh is created
	CommonParameters::Boundary m_locationOfPlane;

	// Index of layer interface
	int m_indexLayerInterface;

	// Additional node IDs on the corners of the analysis domain
	int m_additionalCornerNodes[12];

	// Additional node IDs on the sides of the analysis domain
	std::vector<int> m_additionalSideNodes[20];

	// Array convert node ID of triangle list to the one of boundary curve list
	std::map<int,int> m_nodeIDTriangle2BoundCurve;

	// Array convert node ID of boundary curve list to the one of triangle list 
	std::map<int,int> m_nodeIDBoundCurve2Triangle;

	// Array convert node ID of triangle list to the one of sediment layer
	std::map<int, int> m_nodeLandSurfToSediment;

	// Add super triangle
	void addSuperTriangle();

	// Edge flipping
	void edgeFlipping( const CommonParameters::XY& coord, std::vector<int>& stack );

	// Search triangle enclosing the specified point
	int searchTriangleEnclosingPoint( const CommonParameters::XY& coord, int& onSide, int& locType );

	// Insert new node and flip
	int insertNewNodeAndFlip( const CommonParameters::XY& coord, int& locType, const int nodeIDBoundCurveList );

	// Search for triangles between specified two nodes
	void searchElemBetweenNodes( const int node0, const int node1, std::vector<int>& triangleIDs ) const;

	// Subdivide polygon constructed from specified triangles 
	void subDividePolygon( const int node0, const int node1, const std::vector<int>& triangleIDs );

	// Write triangles to VTK
	void writeTrianglesToVTK( const std::string& fileName );

	// Write PLCs to VTK
	void writePLCsToVTK( const std::string& fileName, const std::vector<CommonParameters::XYZ>& nodeList, const std::vector<TriangleList::Facet>& facetList );

	// Make rough triangles from nodes of boundaries
	void makeRoughTriangles();

	//----- Do not delete for future use >>>>>
	//// Remove out-of-domain triangle 
	//void removeOutOfDomainTriangles();
	//----- Do not delete for future use <<<<<

	// Assign domain type to triangles
	void assignDomainTypeToTriangles();

	// Assign domain type to triangles within the specified boundary curve
	void assignDomainTypeToTrianglesWithinBoundaryCurve( const BoundaryCurve* const ptrBoundaryCurve, const bool clockwise );

	// Search triangles within the specified boundary curve
	void searchTrianglesWithinBoundaryCurve( const BoundaryCurve* const ptrBoundaryCurve, const bool clockwise, std::vector<int>& triangles );

	//// Search boundary curve containing the segment containing specified two nodes as ends
	//const BoundaryCurve* searchBoundaryCurve( const int nodeID0, const int nodeID1, int& nodeIDBoundCurve0, int& nodeIDBoundCurve1 );

	// Get triangle sharing the specified two nodes
	//void getTriangleShareTwoNodes( const int node0, const int node1, int& triLeft, int& triRight, int& edgeLeft, int& edgeRight );
	void getTriangleShareTwoNodes( const int node0, const int node1, int& triLeft, int& triRight );

	// Refine triangles having large edge
	void refineTrianglesHavingLongEdges();

	// Calculate maximum edge length of the specified node
	double calcMaximumEdgeLength( const int iNode );

	// Insert node pair of element having large angle to container
	bool insertNodePairToContainer( const int triangleID, ArrayType& container );

	// Refine triangles having large angle
	void refineTrianglesHavingLargeAngle( const bool notInsertNodesOnBoundary = false );

	// Insert node to the triangles more than two nodes of which locate on the coast line
	void insertNodeToTriangleWithAllNodesOnCoast();

	// Perform laplacian method
	void laplacian();

	// Get node ID of boundary curve list from node ID of triangle list
	int getNodeIDBoundaryCurveList( const int nodeIDTriangleList ) const;

	// Make PLCs
	void makePLCs( std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList );

	// Make surface mesh of the earth
	void makeSurfMeshOfEarth( std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList );

	// Make PLCs for the inner surface inside of land
	void makeFacetOnInnnerSurface( std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList );

	// Make surface mesh of the sea
	void makeSurfMeshOfSea( const std::map<int,int>* const outerBounNode2MeshNode, const std::map<int,int>* const innerBounNode2MeshNode, const std::map<int,int>* const subInnerBounNode2MeshNode,
		std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList );

	// Make surface mesh for sediment layer
	void makeSurfMeshForSedimentLayer( std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList );

	// Make PLCs for the top and bottom
	void makeFacetOnTopAndBot( std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList );

	// Make PLCs for the sides of the sea
	void makeFacetOnSideOfSea( const std::map<int,int>* const outerBounNode2MeshNode, std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList );

	// Make boundary curve of the sides of the sea
	void makeBoundaryCurveOfSideOfSea( std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList );

	// Make PLCs for the surface of the sea
	void makeFacetOnSeaSurface( std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList );

	// Make PLCs for the sides of the air and land
	void makeFacetOnSideOfAirAndLand( const std::map<int,int>* const outerBounNode2MeshNode, const CommonParameters::Boundary& boundaryType, std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList );

	// Make PLCs for lower layers
	void makeFacetOfLowerLayers( std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList );

	// Make PLCs for extended region
	void makeFacetOfExtendedRegion( std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList );

	// Write PLCs to poly file
	void writePLCsToPolyFile( const std::string& fileName, const std::vector<CommonParameters::XYZ>& nodeList, const std::vector<TriangleList::Facet>& facetList );

	// Convert node ID of boundary list to the one of triangle list
	int convertNodeIDBoundCurve2Triangle( const int nodeIDBoundCurve ) const;

	// Convert node ID of triangle to the one of boundary list
	int convertNodeIDTriangle2BoundCurve( const int nodeIDTriangle ) const;

	// Get ID of the triangle which has the two specified nodes.
	// Return -1 if the two specified nodes does not share any triangle.
	int shareSameTriangle( const int nodeID1, const int nodeID2 ) const;

	// Get ID of the triangle which has the two specified nodes except domain boundary.
	// Return -1 if the two specified nodes does not share any triangle or edge connecting the two nodes corresponds to a domann boundary.
	int shareSameTriangleExceptDomainBoundary( const int nodeID1, const int nodeID2 ) const;

	// Subdivide polygon to keep boundary
	void subdividePolygon( std::vector<int>& triangles, const int nodeIDStart , const int nodeIDEnd );

	// Edge flipping
	void flipEdgeOfTwoTriangles( const std::vector<int>& triangles, const int nodeIDStart , const int nodeIDEnd );

	// Delete the data relating to the polygon
	void deleteDataRelatingToPolygon( const std::vector<int>& triangles );

	// Search nodes of polygon
	void searchNodesOfPolygon( const std::vector<int>& triangles, const int nodeIDStart , const int nodeIDEnd, 
		std::vector<int>& nodeGroup1, std::vector<int>& nodeGroup2, std::vector< std::pair<int,int> >& surroundingElemEdges );

	// Fill triangles in polygon
	void fillTrianglesInPolygon( const std::vector<int>& triangles, std::vector<int>& nodeGroup, const int nodeIDStart , const int nodeIDEnd, int& iElem );

	// Reconstruct adjacency relationship
	void reconstructAdjacencyRelationship( const std::vector<int>& triangles, const std::vector< std::pair<int,int> >& surroundingElemEdges );

	// Get the nodes surroudding the analysis domain for setting the extended region
	void getNodesSurroundingAnalysisDomain( const int* additionalCornerNodes, std::vector<CommonParameters::XYZ>& nodeList, std::vector<int>& nodes ) const;

	// Edit height of nodes deviating from the ones of neighbors
	void editDeviatingHeight();

	// Calculate height from the ones of neighbors
	double calcHeightFromTheOnesOfNeighbors( const int iNode ) const;

	// Calculate deviation of height
	void calcDeviationOfHeight( double* deviation ) const;

	// Make boundary curve for the surface meshes of the sea and lakes
	void makeBoundaryCurveForSeaAndLakeMeshes( std::vector<CommonParameters::XYZ>& nodeList, std::map<int,int>*& outerBounNode2MeshNode, std::map<int,int>*& innerBounNode2MeshNode, std::map<int,int>*& subInnerBounNode2MeshNode );

	// Make facet from the surrounding nodes
	void makeFacetFromSurroundingNodes( const std::vector<int>& nodes, const bool reverseOrder, const CommonParameters::Boundary& planeType, 
		const CommonParameters::DomainType& domainType, std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList,
		const int iLayer = -1 );

	// Search the type of the side boundary where specified nodes locate on
	CommonParameters::Boundary searchTypeOfSideBoundary( const std::vector<CommonParameters::XYZ>& nodeList, const std::vector<int>& nodes ) const;

	// Calculate corner coordinates of the land surface
	void calcCornerCoordsOfLandSurface( CommonParameters::XYZ coords[4] );

	//// Add boundary curvers from facet list
	//void addBoundaryCurvesFromFacetList( const std::vector<CommonParameters::XYZ>& nodeList, const std::vector<TriangleList::Facet>& facetList );

	//// Calculate coordinate of 3D model from the coordinate of the plane on which 2D mesh is created
	//CommonParameters::XYZ calcCoordOf3DModel( const CommonParameters::XY& coord2D ) const;

	//// Coordinate transform to the one on XY plane from the one of the plane on which 2D mesh is created
	//CommonParameters::XY coordinateTransform( const CommonParameters::XY& coord2D ) const;

};

#endif
