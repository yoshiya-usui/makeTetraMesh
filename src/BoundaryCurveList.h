//--------------------------------------------------------------------------
// MIT License
//
// Copyright (c) 2021 Yoshiya Usui
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//--------------------------------------------------------------------------
#ifndef DBLDEF_BOUNDARY_CURVE_LIST
#define DBLDEF_BOUNDARY_CURVE_LIST

#include "BoundaryCurveOuter.h"
#include "BoundaryCurveInner.h"
#include "BoundaryCurveSubInner.h"
#include "NodeList.h"
#include <vector>
#include <map>
#include <string>

// Class of the list of boundary 
class BoundaryCurveList{

public:

	// Constructer
	BoundaryCurveList();

	// Destructer
	~BoundaryCurveList();

	// Write boundary curve of 2D mesh to VTK file
	void writeBoundaryCurveListVTK( const std::string& fileName );

	// Make closed boudary curves
	void makeBoundaryCurvesCoastLine();

#ifdef _MOD_FOR_NMT
	// Add NMT dipoles to boudary curve list
	void addNMTDipolesToBoundaryCurveList();
#endif

	// Relate inner boundaries to outer boundaries containing them
	// [note] This function does not support sub-inner boundary
	void relateInnerBoundToOuterBound();

	// Relate sb-inner boundaries to inner boundaries containing them
	// [note] Implementation has not been completed
	void relateSubInnerBoundToInnerBound();

	//----- Don't delete for future use >>>>>
	//// Add boundary curves of anomalies to list
	//// [note] Change future to the function adding boundary curves of lakes
	//void addBoudaryCurveOfAnomalies();
	//----- Don't delete for future use <<<<<

	//// Relate node to boundary curve
	//void relateNodeToBoundaryCurve();

	// Write boundary curve list
	void writeBoundaryCurveList( const std::string& fileName ) const;

	// Write containment relationship of boundary curves
	void writeBoundaryCurveRelatios( const std::string& fileName ) const;

	// Read boundary curve list
	void readBoundaryCurveList( const std::string& fileName );

	// Read containment relationship of boundary curves 
	void readBoundaryCurveRelatios( const std::string& fileName );

	// Set one outer boundary curve
	void setOneOuterBoundaryCurve( const int numNodes, const int* nodeIDs, const int geologicalType );

	// Make all node fixed
	void makeAllNodeFixed();

	// Make node data to VTK file
	void writeNode2DListToVTK( const std::string& fileName ) const;

	// Write node 2D data to the Intermediate file of step 1
	void writeNode2DListToIntermediateFileStep1( const std::string& fileName ) const;

	// Read node list of 2D mesh from intermediate file after step 1
	void readNode2DListFromIntermediateFile( const std::string& fileName );

	// Set node list of 2D mesh
	void setNode2DList( const int numNodes, const CommonParameters::XY* coord2D );

	// Get pointer to the list of the node 2D
	const NodeList* getPointerToNodeList() const;

	// Get total number of the outer boundaries
	int getTotalNumberOuterBoundaries() const;

	// Get total number of the inner boundaries
	int getTotalNumberInnerBoundaries() const;

	// Get total number of the sub-inner boundaries
	int getTotalNumberSubInnerBoundaries() const;

	// Get pointer to the outer boundaries
	const BoundaryCurveOuter* getPointerToOuterBoundary( const int id ) const;

	// Get pointer to the inner boundaries
	const BoundaryCurveInner* getPointerToInnerBoundary( const int id ) const;

	// Get pointer to the sub-inner boundaries
	const BoundaryCurveSubInner* getPointerToSubInnerBoundary( const int id ) const;

	// Get total number of nodes belong to the outer boundaries
	int getNumNodeOuterBoundary( const int iBoun ) const;

	// Get total number of nodes belong to the inner boundaries
	int getNumNodeInnerBoundary( const int iBoun ) const;

	// Get total number of nodes belong to the sub-inner boundaries
	int getNumNodeSubInnerBoundary( const int iBoun ) const;

	// Get coordinate of nodes belong to the outer boundaries
	CommonParameters::XY getPointCoordOuterBoundary( const int iBoun, const int iNode ) const;

	// Get coordinate of nodes belong to the inner boundaries
	CommonParameters::XY getPointCoordInnerBoundary( const int iBoun, const int iNode ) const;

	// Get coordinate of nodes belong to the sub-inner boundaries
	CommonParameters::XY getPointCoordSubInnerBoundary( const int iBoun, const int iNode ) const;

	// Get number of inner boundaries within the specified outer boundary
	int getNumInnerBoundaryIncluded( const int iBounOuter ) const;

	// Get number of sub-inner boundaries within the specified inner boundary
	int getNumSubInnerBoundaryIncluded( const int iBounInner ) const;

	// Get ID of an inner boundary within the specified outer boundary
	int getInnerBoundaryIDIncluded( const int iBounOuter, const int id ) const;

	// Get ID of a sub-inner boundary within the specified inner boundary
	int getSubInnerBoundaryIDIncluded( const int iBounInner, const int id ) const;

	// Get flag specifing whether the specified nodes are the two end of a segment
	bool nodePairOfBoundaryCurve( const int nodeID0, const int nodeID1 );

	// Add new node between segment
	//int addNewNode( const CommonParameters::XY coord );
	int addNewNodeBetweenSegment( const int node0, const int node1 );

	// Remove specified node
	void removeNode( const int nodeID );

	// Search the point of boundary having specified coordinates and return whether the point belong to the boundary curve
	bool hasPointOfGivenCoordInOuterBoundary( const CommonParameters::XY& coord, const int& iBounOuter, int& pointID ) const;

	// Search point of outer boundaries having specified coordinates
	void searchPointOfOuterBoundaryByCoord( const CommonParameters::XY& coord, int& boundCurveID, int& pointID ) const;

	// Search point of outer boundaries having specified coordinates except the specified curve
	void searchPointOfOuterBoundaryByCoordExceptACurve( const CommonParameters::XY& coord, const int exceptCurveID, int& boundCurveID, int& pointID ) const;

	// Exchange the types of boundary curves (Sea => Land, Land => Sea)
	void exchangeTypesOfBoundaryCurves();

	// Setup relation between nodes on boundary and lake data
	void setupRelationBetweenNodesOnBoundaryAndLakeData();

	// Get lake height from coordinate
	int getLakeIndexFromCoordinate( const CommonParameters::XY& coord ) const;

	// Get lake index from node ID on boundary curves
	int getLakeIndexFromNodeID( const int nodeID );

private:

	// Copy constructer
	BoundaryCurveList(const BoundaryCurveList& rhs);

	// Assignment operator
	BoundaryCurveList& operator=(BoundaryCurveList& rhs);

	// Outer boundaries
	std::vector<BoundaryCurveOuter> m_outerBoundaries;

	// Inner boundaries
	std::vector<BoundaryCurveInner> m_innerBoundaries;

	// Sub-inner boundaries
	std::vector<BoundaryCurveSubInner> m_subInnerBoundaries;

	const static double m_eps;

	// List of nodes
	NodeList m_nodeList; 

	// Map relating outer boundary to inner boundary
	std::vector<int>* m_outer2InnerBoundaries;

	// Map relating inner boundary to sub-inner boundary
	std::vector<int>* m_inner2SubInnerBoundaries;

	// Array convert node ID of boundary curve list to lake index
	std::map<int,int> m_nodeIDBoundCurveToLakeIndex;

	// Add node between the specified two points
	void addNodesBetweenTwoPoins( const CommonParameters::XY& startCoord, const CommonParameters::XY& endCoord, const double maxLength, std::vector<int>& nodeIDs );

	// Add node between the specified two points ( maximum edge is calculated automatically )
	void addNodesBetweenTwoPoins( const CommonParameters::XY& startCoord, const CommonParameters::XY& endCoord, std::vector<int>& nodeIDs );

	// Calculate maximum edge length from specified two points
	double calculateMaxEdgeLength( const CommonParameters::XY& startCoord, const CommonParameters::XY& endCoord ) const;

	// Check whether coastline is clockwise
	bool isClockWise(  std::vector<int>& nodeIDVec ) const;

};

#endif
