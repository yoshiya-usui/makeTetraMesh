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
#include "BoundaryCurve.h"
#include "NodeList.h"
#include "OutputFiles.h"
#include "Control.h"
#include "Util.h"
#include "math.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iomanip>

// Constructer
BoundaryCurve::BoundaryCurve( const std::vector<int>& nodes, const int itype ):
	m_nodeIDs( nodes ),
	m_geologicalType( itype )
{

}

// Constructer
BoundaryCurve::BoundaryCurve( const int itype ):
	m_geologicalType( itype )
{

}

// Default constructer
BoundaryCurve::BoundaryCurve():
	m_geologicalType( CommonParameters::UNKNOWN )
{

}

// Destructer
BoundaryCurve::~BoundaryCurve(){

}

// Copy constructer
BoundaryCurve::BoundaryCurve(const BoundaryCurve& rhs){

	m_nodeIDs = rhs.m_nodeIDs;

	m_geologicalType = rhs.m_geologicalType;

}

// Get total number of points
int BoundaryCurve::getNumOfPoints() const{
	
	return static_cast<int>( m_nodeIDs.size() );

}

// Get coordinates of node
CommonParameters::XY BoundaryCurve::getCoordXYOfPoints( const int id, const NodeList* const ptrNode2DList ) const{

	if( id < 0 || id >= getNumOfPoints() ){
		OutputFiles::m_logFile << "Error : ID of point is out of range at the function " << __FUNCTION__  << " !! ID = " << id << std::endl;
		exit(1);
	}

	return ptrNode2DList->getCoordXYOfPoints( m_nodeIDs[id] );
	
}

// Get pointer to the node of the boundary
Node* BoundaryCurve::getPointerToNode( const int id, NodeList* ptrNode2DList ) const{

	if( id < 0 || id >= getNumOfPoints() ){
		OutputFiles::m_logFile << "Error : ID of point is out of range at the function " << __FUNCTION__  << " !! ID = " << id << std::endl;
		exit(1);
	}

	return ptrNode2DList->getPointerToNode( m_nodeIDs[id] );

}

// Get node ID
int BoundaryCurve::getNodeID( const int id ) const{

	if( id < 0 || id >= getNumOfPoints() ){
		OutputFiles::m_logFile << "Error : ID of point is out of range at the function " << __FUNCTION__  << " !! ID = " << id << std::endl;
		exit(1);
	}

	return m_nodeIDs[id];

}

// Get geological type of boundary
int BoundaryCurve::getGeologicalType() const{
	return m_geologicalType;
}

// Get flag this boundary include specified boudary within it
bool BoundaryCurve::include( const BoundaryCurve& rhs, const NodeList* const ptrNode2DList ) const{

	int iLeft(0);
	int iRight(0);
	for( std::vector<int>::const_iterator itrIn = rhs.m_nodeIDs.begin(); itrIn != rhs.m_nodeIDs.end(); ++itrIn ){
		const CommonParameters::XY coordInnerBound = ptrNode2DList->getCoordXYOfPoints(*itrIn);
		addNumLhsAndRhs( coordInnerBound, ptrNode2DList, iLeft, iRight );
	}

#ifdef _DEBUG_WRITE
	std::cout << "iLeft  : " << iLeft << std::endl;
	std::cout << "iRight : " << iRight << std::endl;
#endif

	return include( iLeft, iRight );

}

// Get flag this boundary include the specified point within it
bool BoundaryCurve::include( const CommonParameters::XY& coord, const NodeList* const ptrNode2DList ) const{

	//NodeList* ptrNode2DList = NodeList::getInstance();

	int iLeft(0);
	int iRight(0);

	addNumLhsAndRhs( coord, ptrNode2DList, iLeft, iRight );
		
#ifdef _DEBUG_WRITE
	std::cout << "iLeft  : " << iLeft << std::endl;
	std::cout << "iRight : " << iRight << std::endl;
#endif

	return include( iLeft, iRight );

}

// Add number of segment on which the specified point locate left hand side or right hand side
void BoundaryCurve::addNumLhsAndRhs( const CommonParameters::XY& coord, const NodeList* const ptrNode2DList, int& iLeft, int& iRight ) const{

	//NodeList* ptrNode2DList = NodeList::getInstance();

	int nodeIDPre = m_nodeIDs.back();
	for( std::vector<int>::const_iterator itrOut = m_nodeIDs.begin(); itrOut != m_nodeIDs.end(); ++itrOut ){

		const int nodeIDCur = *itrOut;
		const CommonParameters::XY coordOuterBound = ptrNode2DList->getCoordXYOfPoints(nodeIDCur);

		int icount(0);
		int nodeIDPre2 = m_nodeIDs.back();
		for( std::vector<int>::const_iterator itrOut2 = m_nodeIDs.begin(); itrOut2 != m_nodeIDs.end(); ++itrOut2 ){
			const int nodeIDCur2 = *itrOut2;

			const CommonParameters::XY segmentStart2 = ptrNode2DList->getCoordXYOfPoints(nodeIDPre2);
			const CommonParameters::XY segmentEnd2 = ptrNode2DList->getCoordXYOfPoints(nodeIDCur2);
			if( nodeIDCur != nodeIDPre2 && nodeIDCur != nodeIDCur2 && // Exclude segment including evaluation point
				Util::intersectTwoSegments( coordOuterBound, coord, segmentStart2, segmentEnd2 ) ){
				++icount;
			} 

			nodeIDPre2 = nodeIDCur2;

		}

		if( icount == 0 ){// No intersect point with the other segments
			if( Util::locateRightHandSide( ptrNode2DList->getCoordXYOfPoints(nodeIDPre), coordOuterBound , coord ) ){// Locate left hand side of the segmeht
				++iRight;
			}
			else{
				++iLeft;
			}
		}

		nodeIDPre = nodeIDCur;

	}

}

// Write boundary curve to intermediate file
void BoundaryCurve::writeBoudaryCurve( std::ofstream& ofs ) const{

	if( !ofs.is_open() ){
		std::cerr << " Error : Output file is not open." << std::endl;
		exit(1);
	}

	//NodeList* ptrNode2DList = NodeList::getInstance();

	const int numNode = static_cast<int>( m_nodeIDs.size() );
	ofs << std::setw(10) << numNode << std::setw(10) << m_geologicalType << std::endl;

	for( std::vector<int>::const_iterator itr = m_nodeIDs.begin(); itr != m_nodeIDs.end(); ++itr ){
		const int nodeID = *itr;
		ofs << std::setw(10) << nodeID << std::endl;
	}
	
}

// Read boundary curve from intermediate file
void BoundaryCurve::readBoudaryCurve( std::ifstream& ifs ){

	if( !ifs.is_open() ){
		std::cerr << " Error : Input file has not been opened." << std::endl;
		exit(1);
	}

	int numNode(0);
	ifs >> numNode >> m_geologicalType;
	m_nodeIDs.reserve(numNode);

#ifdef _DEBUG_WRITE
	std::cout << "numNode  : " << numNode << std::endl;
	std::cout << "m_geologicalType : " << m_geologicalType << std::endl;
#endif


	int ibuf(0);
	for( int iNode = 0; iNode < numNode; ++iNode ){
		ifs >> ibuf;
		m_nodeIDs.push_back( ibuf );
	}

}

// Set one boundary curve
void BoundaryCurve::setOneBoundaryCurve( const int numNodes, const int* nodeIDs ){

	for( int iNode = 0; iNode < numNodes; ++iNode ){
		m_nodeIDs.push_back( nodeIDs[iNode] );
	}

}

// Get flag specifing whether the specified nodes are the two end of a segment
bool BoundaryCurve::nodePairOfBoundaryCurve( const int nodeID0, const int nodeID1 ) const{

	std::vector<int>::const_iterator itr0 = find( m_nodeIDs.begin(), m_nodeIDs.end(), nodeID0 );
	std::vector<int>::const_iterator itr1 = find( m_nodeIDs.begin(), m_nodeIDs.end(), nodeID1 );

	if( itr0 != m_nodeIDs.end() && itr1 != m_nodeIDs.end() ){
		const int distance0 = static_cast<int>( distance( m_nodeIDs.begin(), itr0 ) );
		const int distance1 = static_cast<int>( distance( m_nodeIDs.begin(), itr1 ) );
		if( abs( distance0 - distance1 ) == 1 || abs( distance0 - distance1 ) == static_cast<int>( m_nodeIDs.size() ) - 1 ){
			return true;
		}
	}

	return false;

}

// Add new node to the boundary curve
void BoundaryCurve::addNewNode( const int nodeID0, const int nodeID1, const int nodeIDNew ){

	std::vector<int>::iterator itr0 = find( m_nodeIDs.begin(), m_nodeIDs.end(), nodeID0 );
	std::vector<int>::iterator itr1 = find( m_nodeIDs.begin(), m_nodeIDs.end(), nodeID1 );

	if( itr0 != m_nodeIDs.end() && itr1 != m_nodeIDs.end() ){

		const int distance0 = static_cast<int>( distance( m_nodeIDs.begin(), itr0 ) );
		const int distance1 = static_cast<int>( distance( m_nodeIDs.begin(), itr1 ) );
		if( abs( distance0 - distance1 ) == 1 ){

			if( distance1 > distance0 ){
				m_nodeIDs.insert( itr1, nodeIDNew );
			}
			else{
				m_nodeIDs.insert( itr0, nodeIDNew );
			}
			return;

		}
		else if( abs( distance0 - distance1 ) == static_cast<int>( m_nodeIDs.size() ) - 1 ){

			m_nodeIDs.push_back( nodeIDNew );
			return;

		}

	}

	std::cerr << " Error : Specified nodes ( " << nodeID0 << ", " << nodeID1 << " ) are not two end of a segment." << std::endl;
	exit(1);
	
}

// Remove specified node
void BoundaryCurve::removeNode( const int nodeID ){

	std::vector<int>::iterator itr = find( m_nodeIDs.begin(), m_nodeIDs.end(), nodeID );

	if( itr != m_nodeIDs.end() ){
		m_nodeIDs.erase( itr );
	}

}

// Search the point of boundary having specified coordinates and return whether the point belong to the boundary curve
bool BoundaryCurve::hasPointOfGivenCoord( const CommonParameters::XY& coord, const NodeList* const ptrNode2DList, int& pointID ) const{

	const double EPS = 1.0e-12;

	int icount(0);
	for( std::vector<int>::const_iterator itr = m_nodeIDs.begin(); itr != m_nodeIDs.end(); ++itr ){
		const CommonParameters::XY coordTmp = ptrNode2DList->getCoordXYOfPoints( *itr );
	
		if( hypot( coord.X - coordTmp.X, coord.Y - coordTmp.Y ) < EPS ){
			pointID = icount;
			return true;
		}

		++icount;
	}

	pointID = -1;
	return false;

}

// Exchange the types of boundary curves (Sea => Land, Land => Sea)
void BoundaryCurve::exchangeTypesOfBoundaryCurves(){

	if( m_geologicalType == CommonParameters::SEA ){
		m_geologicalType = CommonParameters::LAND;
	}
	else if( m_geologicalType == CommonParameters::LAND ){
		m_geologicalType = CommonParameters::SEA;
	}	

}

// Check the specified points is included in the boudary curve
bool BoundaryCurve::included( const CommonParameters::XY& targetCoord, const NodeList* const ptrNode2DList ) const{

	double angleSum(0.0);
	int nodeIDPre = m_nodeIDs.back();
	for( std::vector<int>::const_iterator itr = m_nodeIDs.begin(); itr != m_nodeIDs.end(); ++itr ){
		const int nodeIDCur = *itr;
		const CommonParameters::XY coord1 = ptrNode2DList->getCoordXYOfPoints(nodeIDPre);
		const CommonParameters::XY coord2 = ptrNode2DList->getCoordXYOfPoints(nodeIDCur);
		const double angle = Util::calcAngle( coord1, targetCoord, coord2 );
		if( Util::calcOuterProduct( targetCoord, coord1, coord2 ) > 0.0 ){
			// Clockwise
			angleSum += angle;
		}else{
			// Counter clockwise
			angleSum -= angle;
		}
		nodeIDPre = nodeIDCur;
	}
	if( !isClockWise( ptrNode2DList) ){
		// Boundary curve is counter clockwise
		angleSum *= -1.0;
	}

	if( fabs( angleSum - CommonParameters::PI * 2.0 ) < 1.0e-6 ){
		return true;
	}else{
		return false;
	}

}

// Check whether boudary curve is clockwise
bool BoundaryCurve::isClockWise( const NodeList* const ptrNode2DList ) const{
	
	double areaSum(0.0);
	int nodeIDPre = m_nodeIDs.back();
	for( std::vector<int>::const_iterator itr = m_nodeIDs.begin(); itr != m_nodeIDs.end(); ++itr ){
		const int nodeIDCur = *itr;
		const CommonParameters::XY coord0 = { 0.0, 0.0 };
		const CommonParameters::XY coord1 = ptrNode2DList->getCoordXYOfPoints(nodeIDPre);
		const CommonParameters::XY coord2 = ptrNode2DList->getCoordXYOfPoints(nodeIDCur);
		areaSum += Util::calcOuterProduct( coord0, coord1, coord2 );
		nodeIDPre = nodeIDCur;
	}
	if( areaSum > 0.0 ){
		return true;
	}
	else{
		return false;
	}
	
}


