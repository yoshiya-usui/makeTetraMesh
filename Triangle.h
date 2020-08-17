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
#ifndef DBLDEF_TRIANGLE
#define DBLDEF_TRIANGLE

#include "Node.h"
#include "NodeList.h"
#include "CommonParameters.h"
#include <vector>

// Class of triangle
class Triangle{

public:

	//enum DomainType{
	//	UNKNOWN = -1,
	//	SEA = 0,
	//	LAND,
	//	LAKE,
	//	OUTSIDE_OF_DOMAIN
	//};

	enum DecisionFlag{
		OUTSIDE_TRIANGLE = -1,
		ON_NODE = 0,
		ON_EDGE,
		INSIDE_TRIANGLE,
	};

	// Default constructer
	explicit Triangle();

	// Constructer
	explicit Triangle( const int* nodes, const int* adjTriangles, const int* edgeOfAdjTriangles, const bool satisfy, const int domainType );

	// Destructer
	~Triangle();

	// Copy constructer
	Triangle( const Triangle& rhs );

	// Assignment operator
	Triangle& operator=(const Triangle& rhs);

	// Get flag specifing whether the specified point locate in the specified triagle
	int doesLocateInTriangle( const CommonParameters::XY& coord, int& onSide, const NodeList* const ptrNode2DList ) const;

	// Get global node ID belongint to the triangle
	int getNodeID( const int num ) const;

	// Get local node ID from global node ID
	int getLocalNodeID( const int globalNodeID ) const;

	// Get adjacent triangle ID
	int getAdjacentTriangleID( const int num ) const;

	// Get edge IDs of adjacent triangle ID
	int getEdgeIDOfAdjacentTriangle( const int num ) const;

	// Set node ID belongint to the triangle
	void setNodeID( const int num, const int id );

	// Set adjacent triangle ID
	void setAdjacentTriangleID( const int num, const int id );

	// Set edge IDs of adjacent triangle ID
	void setEdgeIDOfAdjacentTriangle( const int num, const int id );

	// Get coordinate of gravity center
	CommonParameters::XY getCoordGravCenter( const NodeList* const ptrNode2DList ) const;

	// Set flag specifing whether the triangle satisfies criteria
	void setSatisfyCriteria( const bool flag );

	// Set domain type
	void setDomainType( const int domainType );

	// Get flag specifing whether the triangle satisfies criteria
	bool getSatisfyCriteria() const;

	// Get domain type
	int getDomainType() const;

	// Caluculate determinant
	double calcDeterminant( const NodeList* const ptrNode2DList ) const;

private:

	// Node consisting of triangle
	int m_nodes[3];

	// Adjacent triangles
	int m_adjacentTriangles[3];

	// Adjacent triangles
	int m_edgeOfAdjacentTriangles[3];

	// Flag specifing whether the triangle satisfies criteria
	bool m_satisfyCriteria;

	//// IDs of anomaly to which this triangle belong to.
	//// The value is -1 if this triangle does not belong to any anomaly,
	//int m_anomalyID;

	// Domain type of the triangle
	int m_domainType;

};

#endif
