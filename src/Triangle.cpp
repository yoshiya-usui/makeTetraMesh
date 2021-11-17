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
#include "Triangle.h"
#include "Util.h"
#include "OutputFiles.h"
#include "NodeList.h"
#include "math.h"
#include <iostream>
#include <iomanip>

// Default constructer
Triangle::Triangle():
	m_satisfyCriteria(false),
	m_domainType(CommonParameters::UNKNOWN)
{

	for( int i = 0; i < 3; ++i ){
		m_nodes[i] = -1;
		m_adjacentTriangles[i] = -1;
		m_edgeOfAdjacentTriangles[i] = -1;
	}
	
}

// Constructer
Triangle::Triangle( const int* nodes, const int* adjTriangles, const int* edgeOfAdjTriangles, const bool satisfy, const int domainType ):
	m_satisfyCriteria(satisfy),
	m_domainType(domainType)
{
	for( int i = 0; i < 3; ++i ){
		m_nodes[i] = nodes[i];
		m_adjacentTriangles[i] = adjTriangles[i];
		m_edgeOfAdjacentTriangles[i] = edgeOfAdjTriangles[i];
	}

//#ifdef _DEBUG_WRITE
//	std::cout << "m_edgeOfAdjacentTriangles : ";
//	for( int iNode = 0; iNode < 3; ++iNode ){
//		std::cout << m_edgeOfAdjacentTriangles[iNode] << " ";
//	}
//	std::cout << std::endl;
//#endif

}

// Destructer
Triangle::~Triangle(){

}

// Copy constructer
Triangle::Triangle( const Triangle& rhs ):
	m_satisfyCriteria( rhs.m_satisfyCriteria ),
	m_domainType( rhs.m_domainType )
{

	for( int i = 0; i < 3; ++i ){
		m_nodes[i] = rhs.m_nodes[i];
		m_adjacentTriangles[i] = rhs.m_adjacentTriangles[i];
		m_edgeOfAdjacentTriangles[i] = rhs.m_edgeOfAdjacentTriangles[i];
	}

}

// Assignment operator
Triangle& Triangle::operator=(const Triangle& rhs){

	m_satisfyCriteria = rhs.m_satisfyCriteria;
	m_domainType = rhs.m_domainType;

	for( int i = 0; i < 3; ++i ){
		m_nodes[i] = rhs.m_nodes[i];
		m_adjacentTriangles[i] = rhs.m_adjacentTriangles[i];
		m_edgeOfAdjacentTriangles[i] = rhs.m_edgeOfAdjacentTriangles[i];
	}

	return *this;
}

// Get flag specifing whether the specified point locate in the specified triagle
int Triangle::doesLocateInTriangle( const CommonParameters::XY& coord, int& onSide, const NodeList* const ptrNode2DList ) const{

	const double EPS = 1.0e-9;

	onSide = -1;

	//Node2DList* ptrNode2DList = Node2DList::getInstance();

	const CommonParameters::XY coordTri[3] = {
		ptrNode2DList->getCoordXYOfPoints( m_nodes[0] ),
		ptrNode2DList->getCoordXYOfPoints( m_nodes[1] ),
		ptrNode2DList->getCoordXYOfPoints( m_nodes[2] )
	};

//#ifdef _DEBUG_WRITE
//	std::cout << "coordTri[0] = " << coordTri[0].X << " " << coordTri[0].Y << std::endl;
//	std::cout << "coordTri[1] = " << coordTri[1].X << " " << coordTri[1].Y << std::endl;
//	std::cout << "coordTri[2] = " << coordTri[2].X << " " << coordTri[2].Y << std::endl;
//#endif

	const double det0 = Util::calcOuterProduct( coordTri[0], coordTri[1], coord  );
	const double det1 = Util::calcOuterProduct( coordTri[1], coordTri[2], coord  );
	const double det2 = Util::calcOuterProduct( coordTri[2], coordTri[0], coord  );
	//const double detT = Util::calcOuterProduct( coordTri[0], coordTri[1], coordTri[2] );

//#ifdef _DEBUG_WRITE
//	std::cout << "det0 = " << det0 << std::endl;
//	std::cout << "det1 = " << det1 << std::endl;
//	std::cout << "det2 = " << det2 << std::endl;
//	//std::cout << "detT = " << detT << std::endl;
//	//std::cout << "diff = " << fabs( detT - det0 - det1 - det2 ) << std::endl;
//#endif

	//if( fabs( detT - det0 - det1 - det2 ) < EPS ){// Locate in triangle
	if( det0 > -EPS && det1 > -EPS && det2 > -EPS ){// Locate in triangle

//#ifdef _DEBUG_WRITE
//		std::cout << "Locate in triangle" << std::endl;
//#endif

		// Specified point is a duplicate of one of the corner nodes
		if( fabs(det0) < EPS && fabs(det1) < EPS ){
			onSide = 1;
			return Triangle::ON_NODE;
		}
		else if( fabs(det1) < EPS && fabs(det2) < EPS ){
			onSide = 2;
			return Triangle::ON_NODE;
		}
		else if( fabs(det2) < EPS && fabs(det0) < EPS ){
			onSide = 0;
			return Triangle::ON_NODE;
		}

		// Specified point locate one of the edge of the triangle
		if( fabs(det0) < EPS ){
			onSide = 0;
			return Triangle::ON_EDGE;
		}
		else if( fabs(det1) < EPS ){
			onSide = 1;
			return Triangle::ON_EDGE;
		}
		else if( fabs(det2) < EPS ){
			onSide = 2;
			return Triangle::ON_EDGE;
		}

		// Specified point locate completely inside of the triangle
		onSide = -1;
		return Triangle::INSIDE_TRIANGLE;

	}
	
	// Specified point locate out of the triangle
//#ifdef _DEBUG_WRITE
//	std::cout << "Locate out of triangle" << std::endl;
//#endif

	const CommonParameters::XY centerCoord = Util::calcCenterOfGravity( coordTri[0], coordTri[1], coordTri[2] );

//#ifdef _DEBUG_WRITE
//	std::cout << "centerCoord : " << centerCoord.X << " " << centerCoord.Y << std::endl;
//#endif

	for( int i = 0; i < 3; ++i ){
		if( Util::intersectTwoSegments( coordTri[i], coordTri[(i+1)%3], centerCoord, coord ) ){
			onSide = m_adjacentTriangles[i];
//#ifdef _DEBUG_WRITE
//			std::cout << "i onSide : " << i << " " <<  onSide << std::endl;
//#endif
			if( onSide >= 0 ){
				break;
			}
		}
	}
//
//#ifdef _DEBUG_WRITE
//	std::cout << "onSide : " << onSide << std::endl;
//#endif

	return Triangle::OUTSIDE_TRIANGLE;


}

// Get node ID belongint to the triangle
int Triangle::getNodeID( const int num ) const{

	if( num < 0 && num > 2 ){
		OutputFiles::m_logFile << " Error : Specified number is out of the range. num = " << num << std::endl;
		exit(1);
	}

	return m_nodes[num];

}

// Get local node ID from global node ID
int Triangle::getLocalNodeID( const int globalNodeID ) const{

	for( int iNode = 0; iNode < 3; ++iNode ){
		if( m_nodes[iNode] == globalNodeID ){
			return iNode;
		}
	}

	return -1;

}

// Get adjacent triangle ID
int Triangle::getAdjacentTriangleID( const int num ) const{

	if( num < 0 && num > 2 ){
		OutputFiles::m_logFile << " Error : Specified number is out of the range. num = " << num << std::endl;
		exit(1);
	}

	return m_adjacentTriangles[num];

}

// Get edge IDs of adjacent triangle ID
int Triangle::getEdgeIDOfAdjacentTriangle( const int num ) const{

	if( num < 0 && num > 2 ){
		OutputFiles::m_logFile << " Error : Specified number is out of the range. num = " << num << std::endl;
		exit(1);
	}

	return m_edgeOfAdjacentTriangles[num];

}

// Set node ID belongint to the triangle
void Triangle::setNodeID( const int num, const int id ){

	if( num < 0 && num > 2 ){
		OutputFiles::m_logFile << " Error : Specified number is out of the range. num = " << num << std::endl;
		exit(1);
	}

	m_nodes[num] = id;

}

// Set adjacent triangle ID
void Triangle::setAdjacentTriangleID( const int num, const int id ){

	if( num < 0 && num > 2 ){
		OutputFiles::m_logFile << " Error : Specified number is out of the range. num = " << num << std::endl;
		exit(1);
	}

	m_adjacentTriangles[num] = id;

}

// Set edge IDs of adjacent triangle ID
void Triangle::setEdgeIDOfAdjacentTriangle( const int num, const int id ){

	if( num < 0 && num > 2 ){
		OutputFiles::m_logFile << " Error : Specified number is out of the range. num = " << num << std::endl;
		exit(1);
	}

	m_edgeOfAdjacentTriangles[num] = id;

}

// Get coordinate of gravity center
CommonParameters::XY Triangle::getCoordGravCenter( const NodeList* const ptrNode2DList ) const{
	
	return Util::calcCenterOfGravity( ptrNode2DList->getCoordXYOfPoints( m_nodes[0] ), ptrNode2DList->getCoordXYOfPoints( m_nodes[1] ), ptrNode2DList->getCoordXYOfPoints( m_nodes[2] ) );

}

// Set flag specifing whether the triangle satisfies criteria
void Triangle::setSatisfyCriteria( const bool flag ){
	m_satisfyCriteria = flag;
}

// Set domain type
void Triangle::setDomainType( const int domainType ){

	if( domainType < CommonParameters::UNKNOWN || domainType > CommonParameters::OUTSIDE_OF_DOMAIN ){
		OutputFiles::m_logFile << " Error : Wrong domaint type. domainType = " << domainType << std::endl;
		exit(1);
	}

	m_domainType = domainType;

}

// Get domain type
bool Triangle::getSatisfyCriteria() const{
	return m_satisfyCriteria;
}

// Get domain type
int Triangle::getDomainType() const{
	return m_domainType;
}

// Caluculate jacobian
double Triangle::calcDeterminant( const NodeList* const ptrNode2DList ) const{

	return Util::calcOuterProduct( ptrNode2DList->getCoordXYOfPoints( m_nodes[0] ), ptrNode2DList->getCoordXYOfPoints( m_nodes[1] ), ptrNode2DList->getCoordXYOfPoints( m_nodes[2] ) );

}
