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
#ifndef DBLDEF_BOUNDARY_CURVE
#define DBLDEF_BOUNDARY_CURVE

#include "Node.h"
#include "NodeList.h"
#include <vector>
#include <fstream>
#include <iostream>
#include "CommonParameters.h"

// Class of boundary of two dimensional field
class BoundaryCurve{

public:

	const static int UNKNOWN = -1;

	//enum GeologicalType{
	//	SEA = 0,
	//	LAND,
	//	LAKE,
	//};

	enum BoundaryType{
		OUTER = 0,
		INNER,
		SUB_INNER,
	};

	// Constructer
	explicit BoundaryCurve( const std::vector<int>& coords, const int itype );

	// Constructer
	explicit BoundaryCurve( const int itype );

	// Default constructer
	explicit BoundaryCurve();

	// Destructer
	virtual ~BoundaryCurve();

	// Copy constructer
	explicit BoundaryCurve(const BoundaryCurve& rhs);
	
	// Get total number of points
	int getNumOfPoints() const;

	// Get coordinates of point
	CommonParameters::XY getCoordXYOfPoints( const int id, const NodeList* const ptrNode2DList ) const;

	// Get pointer to the node of the boundary
	Node* getPointerToNode( const int id, NodeList* ptrNode2DList ) const;

	// Get node ID
	int getNodeID( const int id ) const;

	// Get geological type of boundary
	int getGeologicalType() const;

	// Get flag this boundary include the specified boudary within it
	//virtual bool include( const BoundaryCurve& rhs, const NodeList* const ptrNode2DList ) const;
	bool include( const BoundaryCurve& rhs, const NodeList* const ptrNode2DList ) const;

	// Get flag this boundary include the specified point within it
	//virtual bool include( const CommonParameters::XY& coord, const NodeList* const ptrNode2DList ) const;
	bool include( const CommonParameters::XY& coord, const NodeList* const ptrNode2DList ) const;

	// Add number of segment on which the specified point locate left hand side or right hand side
	void addNumLhsAndRhs( const CommonParameters::XY& coord, const NodeList* const ptrNode2DList, int& iLeft, int& iRight ) const;

	// Get flag specifing whether object locate within the boundary from 
	// number of segment on which the object locate left hand side or right hand side
	virtual bool include( const int iLeft, const int iRight ) const = 0;

	// Write boundary curve to intermediate file
	void writeBoudaryCurve( std::ofstream& ofs ) const;

	// Read  boundary curve from intermediate file
	void readBoudaryCurve( std::ifstream& ifs );

	// Set one boundary curve
	void setOneBoundaryCurve( const int numNodes, const int* nodeIDs );

	// Get flag specifing whether the specified nodes are the two end of a segment
	bool nodePairOfBoundaryCurve( const int nodeID0, const int nodeID1 ) const;

	// Add new node to the boundary curve
	void addNewNode( const int nodeID0, const int nodeID1, const int nodeIDNew );
	
	// Remove specified node
	void removeNode( const int nodeID );

	// Search the point of boundary having specified coordinates and return whether the point belong to the boundary curve
	bool hasPointOfGivenCoord( const CommonParameters::XY& coord, const NodeList* const ptrNode2DList, int& pointID ) const;

	// Exchange the types of boundary curves (Sea => Land, Land => Sea)
	void exchangeTypesOfBoundaryCurves();

	// Check the specified points is included in the boudary curve
	bool included( const CommonParameters::XY& targetCoord, const NodeList* const ptrNode2DList ) const;

	// Check whether boudary curve is clockwise
	bool isClockWise( const NodeList* const ptrNode2DList ) const;

protected:

	// Nodes of boundary
	std::vector<int> m_nodeIDs;

	// Geological type of boundary
	int m_geologicalType;

private:

	// Assignment operator
	BoundaryCurve& operator=(BoundaryCurve& rhs);

};

#endif
