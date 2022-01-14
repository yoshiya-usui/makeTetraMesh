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
#ifndef DBLDEF_NODE_2D
#define DBLDEF_NODE_2D

#include "CommonParameters.h"
#include <vector>
#include <set>

// Class of node of two dimensinal field
class Node{

public:
	enum Location{
		UNKNOWN = -1,
		COAST_LINE = 0,
		SEA,
		LAND,
		LAKE,
		LAKE_LINE,
		OUT_OF_DOMAIN = 999,
	};

	// Default constructer
	Node();

	// Constructer
	Node( const CommonParameters::XY& coord, const bool bfix );

	// Constructer
	Node( const CommonParameters::XYZ& coord, const bool bfix, const int location );

	// Destructer
	~Node();

	// Copy constructer
	Node(const Node& rhs);

	// Assignment operator
	Node& operator=(const Node& rhs);

	// Set X Y coordinate
	void setCoordXY( const CommonParameters::XY& coord );

	// Set Z coordinate
	void setCoordZ( const double coordZ );

	// Get X Y coordinate
	CommonParameters::XY getCoordXY() const;

	// Get X Y Z coordinate
	CommonParameters::XYZ getCoordXYZ() const;

	// Get X coordinate
	double getCoordX() const;

	// Get Y coordinate
	double getCoordY() const;

	// Get Z coordinate
	double getCoordZ() const;

	// Get flag specifing this node is fixed or not
	bool isFixed() const;

	// Calculate maximum edge length of this points assming this is one of the point of the coast line
	double calcMaximumEdgeLengthForCoastLine() const;

	// Calculate distance of two node in horizontal plane
	double calcHorizontalDistance(const Node& rhs) const;

	// Make this node to be fixed
	void fixThisNode();

	//// Write nodes of 2D mesh to intermediate file
	//void writeNode2DToIntermediateFile( std::ofstream& ofs ) const;

	//// Add IDs of boundary curve sharing the node
	//void addBoundaryCurveID( const int itype, const int ID );

	// Add IDs of triangles sharing the node
	void addTriangleIDs( const int triangleID );

	// Erase IDs of triangles sharing the node
	void eraseTriangleIDs( const int triangleID );

	// Get IDs of triangles sharing the node
	std::vector<int> getTriangleIDs() const;

	// Set location of the node
	void setLocation( const int loc );

	// Get location of the node
	int getLocation() const;

private:
	// Coordinates of the node
	CommonParameters::XYZ m_coord;

	// Flag specifing the node is fixed or not
	bool m_fixed;

	// Location of the node
	int m_location;

	//// IDs of outer boundary curves sharing the node
	//std::vector<int> m_outerBoundaryCurves;

	//// IDs of inner boundary curves sharing the node
	//std::vector<int> m_innerBoundaryCurves;

	//// IDs of sub-inner boundary curves sharing the node
	//std::vector<int> m_subInnerBoundaryCurves;

	// IDs of triangles sharing the node
	std::vector<int> m_triangleIDs;

};

#endif
