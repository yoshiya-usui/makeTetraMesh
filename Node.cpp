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
#include "Node.h"
#include "Control.h"
#include "ObservingSiteList.h"
#include "OutputFiles.h"
#include "math.h"
#include "Util.h"
#include "AnalysisDomain.h"
#include "CoastLineList.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

// Default constructer
Node::Node():
	m_fixed( false ),
	m_location( Node::UNKNOWN )
{

	m_coord.X = 0.0;
	m_coord.Y = 0.0;
	m_coord.Z = 0.0;

}

// Constructer
Node::Node( const CommonParameters::XY& coord, const bool bfix ):
	m_fixed( bfix ),
	m_location( Node::UNKNOWN )
{
	m_coord.X = coord.X;
	m_coord.Y = coord.Y;
	m_coord.Z = 0.0;
}

// Constructer
Node::Node( const CommonParameters::XYZ& coord, const bool bfix, const int location ):
	m_fixed( bfix ),
	m_location( location )
{
	m_coord.X = coord.X;
	m_coord.Y = coord.Y;
	m_coord.Z = coord.Z;
}

// Destructer
Node::~Node(){

}

// Copy constructer
Node::Node(const Node& rhs){

	m_coord = rhs.m_coord;

	m_fixed = rhs.m_fixed;

	m_location = rhs.m_location;

	if( !m_triangleIDs.empty() ){
		m_triangleIDs.clear();
	}

	m_triangleIDs = rhs.m_triangleIDs;
	
}

// Assignment operator
Node& Node::operator=(const Node& rhs){

    if(this == &rhs){
		return *this;
	}

	m_coord = rhs.m_coord;

	m_fixed = rhs.m_fixed;

	m_location = rhs.m_location;

	if( !m_triangleIDs.empty() ){
		m_triangleIDs.clear();
	}

	m_triangleIDs = rhs.m_triangleIDs;

	return *this;

}

// Set X Y coordinate
void Node::setCoordXY( const CommonParameters::XY& coord ){
	m_coord.X = coord.X;
	m_coord.Y = coord.Y;
}

// Set Z coordinate
void Node::setCoordZ( const double coordZ ){
	m_coord.Z = coordZ;
}

// Get X Y coordinate
CommonParameters::XY Node::getCoordXY() const{
	CommonParameters::XY coord2D = { m_coord.X, m_coord.Y };
	return coord2D;
}

// Get X Y Z coordinate
CommonParameters::XYZ Node::getCoordXYZ() const{
	CommonParameters::XYZ coord3D = { m_coord.X, m_coord.Y, m_coord.Z };
	return coord3D;
}

// Get X coordinate
double Node::getCoordX() const{
	return m_coord.X;
}

// Get X coordinate
double Node::getCoordY() const{
	return m_coord.Y;
}

// Get Z coordinate
double Node::getCoordZ() const{
	return m_coord.Z;
}

// Get flag specifing this node is fixed or not
bool Node::isFixed() const{
	return m_fixed;
}

// Calculate maximum edge length of this points assming this is one of the point of the coast line
double Node::calcMaximumEdgeLengthForCoastLine() const{

	return ( CoastLineList::getInstance() )->calcMaximumEdgeLengthForCoastLine( getCoordXY() );

}

// Calculate distance of two node in horizontal plane
double Node::calcHorizontalDistance(const Node& rhs) const{

	return hypot( m_coord.X - (rhs.m_coord).X, m_coord.Y - (rhs.m_coord).Y );
	
}

// Make this node to be fixed
void Node::fixThisNode(){

	m_fixed = true;

}

//// Write nodes of 2D mesh to intermediate file
//void Node::writeNode2DToIntermediateFile( std::ofstream& ofs ) const{
//
//	if( !ofs.is_open() ){
//		std::cerr << " Error : Output file is not open." << std::endl;
//		exit(1);
//	}
//
//	ofs.precision(9);
//
//	ofs << std::setw(20) << std::scientific << m_coord.X
//		<< std::setw(20) << std::scientific << m_coord.Y
//		<< std::setw(20) << std::scientific << 0.0 << std::endl;
//
//	const int ifix = m_fixed ? 1 : 0;
//
//	ofs << std::setw(10) << ifix << std::endl;
//
//	//ofs << std::setw(10) << static_cast<int>( m_outerBoundaryCurves.size() );
//	//for( std::vector<int>::const_iterator itr = m_outerBoundaryCurves.begin(); itr != m_outerBoundaryCurves.end(); ++itr ){
//	//	ofs << std::setw(10) << *itr;
//	//}
//	//ofs << std::endl;
//
//	//ofs << std::setw(10) << static_cast<int>( m_innerBoundaryCurves.size() );
//	//for( std::vector<int>::const_iterator itr = m_innerBoundaryCurves.begin(); itr != m_innerBoundaryCurves.end(); ++itr ){
//	//	ofs << std::setw(10) << *itr;
//	//}
//	//ofs << std::endl;
//
//	//ofs << std::setw(10) << static_cast<int>( m_subInnerBoundaryCurves.size() );
//	//for( std::vector<int>::const_iterator itr = m_subInnerBoundaryCurves.begin(); itr != m_subInnerBoundaryCurves.end(); ++itr ){
//	//	ofs << std::setw(10) << *itr;
//	//}
//	//ofs << std::endl;
//
//	//ofs << std::setw(10) << static_cast<int>( m_triangleIDs.size() );
//	//for( std::vector<int>::const_iterator itr = m_triangleIDs.begin(); itr != m_triangleIDs.end(); ++itr ){
//	//	ofs << std::setw(10) << *itr;
//	//}
//	//ofs << std::endl;
//
//}

//// Add IDs of boundary curve sharing the node
//void Node::addBoundaryCurveID( const int itype, const int ID ){
//
//	if( ID < 0 ){
//		std::cerr << " Error : Boundary curve ID is negative !! ID = " << ID << std::endl;
//		exit(1);
//	}
//
//	switch( itype ){
//		case BoundaryCurve::OUTER:
//			if( std::find( m_outerBoundaryCurves.begin(), m_outerBoundaryCurves.end(), ID ) == m_outerBoundaryCurves.end() ){
//				// Specified ID has not been inserted yet
//				m_outerBoundaryCurves.push_back(ID);
//			}
//			break;
//		case BoundaryCurve::INNER:
//			if( std::find( m_innerBoundaryCurves.begin(), m_innerBoundaryCurves.end(), ID ) == m_innerBoundaryCurves.end() ){
//				// Specified ID has not been inserted yet
//				m_innerBoundaryCurves.push_back(ID);
//			}
//			break;
//		case BoundaryCurve::SUB_INNER:
//			if( std::find( m_subInnerBoundaryCurves.begin(), m_subInnerBoundaryCurves.end(), ID ) == m_subInnerBoundaryCurves.end() ){
//				// Specified ID has not been inserted yet
//				m_subInnerBoundaryCurves.push_back(ID);
//			}
//			break;
//		default:
//			std::cerr << "Type of boundary is wrong !! itype = " << itype << std::endl;
//			exit(1);
//			break;
//	}
//
//}

// Add IDs of triangles sharing the node
void Node::addTriangleIDs( const int triangleID ){
	if( find( m_triangleIDs.begin(), m_triangleIDs.end(), triangleID ) == m_triangleIDs.end() ){
		m_triangleIDs.push_back( triangleID );
	}
}

// Erase IDs of triangles sharing the node
void Node::eraseTriangleIDs( const int triangleID ){

	std::vector<int>::iterator itr = find( m_triangleIDs.begin(), m_triangleIDs.end(), triangleID );
	if( itr != m_triangleIDs.end() ){// Found
		m_triangleIDs.erase( itr );
	}

}

// Get IDs of triangles sharing the node
std::vector<int> Node::getTriangleIDs() const{

	return m_triangleIDs;

}

// Set location of the node
void Node::setLocation( const int loc ){
	m_location = loc;
}

// Get location of the node
int Node::getLocation() const{
	return m_location;
}
