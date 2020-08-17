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
#include "AnalysisDomain.h"
#include "OutputFiles.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <math.h>

const double AnalysisDomain::m_eps = 1.0e-9;

// Return the instance of the class
AnalysisDomain* AnalysisDomain::getInstance(){
   	static AnalysisDomain instance;// The only instance
  	return &instance;
}

// Default constructer
AnalysisDomain::AnalysisDomain():
	m_maxCoordX(-1.0),
	m_maxCoordY(-1.0),
	m_maxCoordZ(-1.0),
	m_minCoordX(-1.0),
	m_minCoordY(-1.0),
	m_minCoordZ(-1.0)
{

}

// Destructer
AnalysisDomain::~AnalysisDomain(){

}

// Read data of analysis domain from input file
void AnalysisDomain::readAnalysisDomainData(){

	std::ifstream ifs( "analysis_domain.dat", std::ios::in );

	if( ifs.fail() ){
		OutputFiles::m_logFile << "Error : File open error : analysis_domain.dat !!" << std::endl;
		exit(1);
	}

	OutputFiles::m_logFile << "# Read data of analysis domain from input file." << std::endl;

	ifs >> m_minCoordX >> m_maxCoordX;
	ifs >> m_minCoordY >> m_maxCoordY;
	ifs >> m_minCoordZ >> m_maxCoordZ;

	if( m_minCoordX >= m_maxCoordX ){
		OutputFiles::m_logFile << "Error : Minimum value of X coordinate ( " << m_minCoordX << " ) is greater than or equal to maximum of X coordinate ( " << m_maxCoordX << " )" << std::endl;
		exit(1);
	}

	if( m_minCoordY >= m_maxCoordY ){
		OutputFiles::m_logFile << "Error : Minimum value of Y coordinate ( " << m_minCoordY << " ) is greater than or equal to maximum of Y coordinate ( " << m_maxCoordY << " )" << std::endl;
		exit(1);
	}

	if( m_minCoordZ >= m_maxCoordZ ){
		OutputFiles::m_logFile << "Error : Minimum value of Z coordinate ( " << m_minCoordZ << " ) is greater than or equal to maximum of Z coordinate ( " << m_maxCoordZ << " )" << std::endl;
		exit(1);
	}

	OutputFiles::m_logFile << "# Minimum value of X coordinate [km] : " << m_minCoordX << std::endl;
	OutputFiles::m_logFile << "# Maximum value of X coordinate [km] : " << m_maxCoordX << std::endl;
	OutputFiles::m_logFile << "# Minimum value of Y coordinate [km] : " << m_minCoordY << std::endl;
	OutputFiles::m_logFile << "# Maximum value of Y coordinate [km] : " << m_maxCoordY << std::endl;
	OutputFiles::m_logFile << "# Minimum value of Z coordinate [km] : " << m_minCoordZ << std::endl;
	OutputFiles::m_logFile << "# Maximum value of Z coordinate [km] : " << m_maxCoordZ << std::endl;

	ifs.close();

}

// Get minimum value of X coordinate
double AnalysisDomain::getMinCoordX() const{
	return m_minCoordX;
}

// Get maximum value of X coordinate
double AnalysisDomain::getMaxCoordX() const{
	return m_maxCoordX;
}

// Get minimum value of Y coordinate
double AnalysisDomain::getMinCoordY() const{
	return m_minCoordY;
}

// Get maximum value of Y coordinate
double AnalysisDomain::getMaxCoordY() const{
	return m_maxCoordY;
}

// Get minimum value of Z coordinate
double AnalysisDomain::getMinCoordZ() const{
	return m_minCoordZ;
}

// Get maximum value of Z coordinate
double AnalysisDomain::getMaxCoordZ() const{
	return m_maxCoordZ;
}

// Get edge ID of the analysis domain where specified point locate
int AnalysisDomain::getEdgeID( const CommonParameters::XY& coord ) const{

	const double eps = 1.0e-9;

	if( fabs( coord.X - m_maxCoordX ) < eps ){
		return AnalysisDomain::PLUS_X;
	}
	else if( fabs( coord.X - m_minCoordX ) < eps ){
		return AnalysisDomain::MINUS_X;
	}
	else if( fabs( coord.Y - m_maxCoordY ) < eps ){
		return AnalysisDomain::PLUS_Y;
	}
	else if( fabs( coord.Y - m_minCoordY ) < eps ){
		return AnalysisDomain::MINUS_Y;
	}

	OutputFiles::m_logFile << "Error : Specified point ( " << coord.X << " , " << coord.Y << " ) does't locate on the boundary of the analysis domain." << std::endl;
	exit(1);

	return AnalysisDomain::BAD_DATA;

}

// Get flag specifing whether the direction from start point to end point is colock wise on the boundary of analysis domain
bool AnalysisDomain::isCWROnBoundary( const CommonParameters::XY& startCoord, const CommonParameters::XY& endCoord ) const{

	const int startEdge = getEdgeID( startCoord ); 
	const int endEdge = getEdgeID( endCoord ); 

	if( startEdge != endEdge ){
		OutputFiles::m_logFile << "Error : Start point locate on edge " << startEdge << " of analysis domain although end point locate on edge " << endEdge << " !!" << std::endl;
		exit(1);
	}

	switch (startEdge)
	{
		case AnalysisDomain::PLUS_X:
			if( endCoord.Y > startCoord.Y ){
				return true;
			}
			break;
		case AnalysisDomain::MINUS_X:
			if( endCoord.Y < startCoord.Y ){
				return true;
			}
			break;
		case AnalysisDomain::PLUS_Y:
			if( endCoord.X < startCoord.X ){
				return true;
			}
			break;
		case AnalysisDomain::MINUS_Y:
			if( endCoord.X > startCoord.X ){
				return true;
			}
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong edge ID !! ID = " << startEdge << " !!" << std::endl;
			exit(1);
			break;
	}

	return false;

}

// Get start coordinate of the edge of the analysis domain
CommonParameters::XY AnalysisDomain::getStartCoord( const int edgeID ) const{

	CommonParameters::XY coord;

	switch(edgeID)
	{
		case AnalysisDomain::PLUS_X:
			coord.X = m_maxCoordX;
			coord.Y = m_minCoordY;
			break;
		case AnalysisDomain::PLUS_Y:
			coord.X = m_maxCoordX;
			coord.Y = m_maxCoordY;
			break;
		case AnalysisDomain::MINUS_X:
			coord.X = m_minCoordX;
			coord.Y = m_maxCoordY;
			break;
		case AnalysisDomain::MINUS_Y:
			coord.X = m_minCoordX;
			coord.Y = m_minCoordY;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong edge ID !! ID = " << edgeID << " !!" << std::endl;
			exit(1);
			break;
	}

	return coord;

}

// Get end coordinate of the edge of the analysis domain
CommonParameters::XY AnalysisDomain::getEndCoord( const int edgeID ) const{

	switch(edgeID)
	{
		case AnalysisDomain::PLUS_X:
			// Go through next
		case AnalysisDomain::PLUS_Y:
			// Go through next
		case AnalysisDomain::MINUS_X:
			// Go through next
		case AnalysisDomain::MINUS_Y:
			return getStartCoord( getNextEdgeIDCWR( edgeID ) );
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong edge ID !! ID = " << edgeID << " !!" << std::endl;
			exit(1);
			break;
	}

	CommonParameters::XY coord = { 0. , 0. };
	return coord;

}

// Get ID of next edge in colock wise order
int AnalysisDomain::getNextEdgeIDCWR( const int edgeID ) const{

	switch(edgeID)
	{
		case AnalysisDomain::PLUS_X:
			return PLUS_Y;
			break;
		case AnalysisDomain::PLUS_Y:
			return MINUS_X;
			break;
		case AnalysisDomain::MINUS_X:
			return MINUS_Y;
			break;
		case AnalysisDomain::MINUS_Y:
			return PLUS_X;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong edge ID !! ID = " << edgeID << " !!" << std::endl;
			exit(1);
			break;
	}

	return -1;

}

// Get ID of next edge in anticolock wise order
int AnalysisDomain::getNextEdgeIDAntiCWR( const int edgeID ) const{

	switch(edgeID)
	{
		case AnalysisDomain::PLUS_X:
			return MINUS_Y;
			break;
		case AnalysisDomain::PLUS_Y:
			return PLUS_X;
			break;
		case AnalysisDomain::MINUS_X:
			return PLUS_Y;
			break;
		case AnalysisDomain::MINUS_Y:
			return MINUS_X;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong edge ID !! ID = " << edgeID << " !!" << std::endl;
			exit(1);
			break;
	}

	return -1;

}

// Calculate distance on the boundary of analysis domain from upper left point
double AnalysisDomain::calcDistanceOnBoundaryFromUpperLeft( const CommonParameters::XY& coord ) const{

	const double edgeLengthX =  m_maxCoordX - m_minCoordX;
	const double edgeLengthY =  m_maxCoordY - m_minCoordY;

	const int edgeID = getEdgeID( coord );
	switch( edgeID ){
		case AnalysisDomain::PLUS_X:
			return coord.Y - m_minCoordY;
			break;
		case AnalysisDomain::PLUS_Y:
			return edgeLengthY + m_maxCoordX - coord.X;
			break;
		case AnalysisDomain::MINUS_X:
			return edgeLengthY + edgeLengthX + m_maxCoordY - coord.Y;
			break;
		case AnalysisDomain::MINUS_Y:
			return 2.0 * edgeLengthY + edgeLengthX + coord.X - m_minCoordX;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong edge ID !! ID = " << edgeID << " !!" << std::endl;
			exit(1);
			break;
	}

	return 0.0;

}

// Get the coordinate of the next corner when traveling in the clockwise direction
CommonParameters::XY AnalysisDomain::getNextCoordCWRDirection( const CommonParameters::XY& coord ) const{

	return getEndCoord( getEdgeID( coord ) );

}

// Get the coordinate of the next corner when traveling in the anticlockwise direction
CommonParameters::XY AnalysisDomain::getNextCoordACWRDirection( const CommonParameters::XY& coord ) const{

	return getStartCoord( getEdgeID( coord ) );

}

// Get flag whether specified point intersects with the boundary of the analysis domain
bool AnalysisDomain::doesIntersectWithBoundary( const CommonParameters::XY& coord ) const{

	if( doesIntersectWithPlusXEdgeOfAnalysisDomain(coord) || doesIntersectWithMinusXEdgeOfAnalysisDomain(coord) ||
		doesIntersectWithPlusYEdgeOfAnalysisDomain(coord) || doesIntersectWithMinusYEdgeOfAnalysisDomain(coord) ){
		return true;
	}

	return false;

}

// Get flag whether specified point intersects with the edge of + X side of the analysis domain
bool AnalysisDomain::doesIntersectWithPlusXEdgeOfAnalysisDomain( const CommonParameters::XY& coord ) const{
	return ( fabs( coord.X - m_maxCoordX ) < m_eps && m_minCoordY <= coord.Y && coord.Y <= m_maxCoordY );
}

// Get flag whether specified point intersects with the edge of - X side of the analysis domain
bool AnalysisDomain::doesIntersectWithMinusXEdgeOfAnalysisDomain( const CommonParameters::XY& coord ) const{
	return ( fabs( coord.X - m_minCoordX ) < m_eps && m_minCoordY <= coord.Y && coord.Y <= m_maxCoordY );
}

// Get flag whether specified point intersects with the edge of + Y side of the analysis domain
bool AnalysisDomain::doesIntersectWithPlusYEdgeOfAnalysisDomain( const CommonParameters::XY& coord ) const{
	return ( fabs( coord.Y - m_maxCoordY ) < m_eps && m_minCoordX <= coord.X && coord.X <= m_maxCoordX );
}

// Get flag whether specified point intersects with the edge of - Y side of the analysis domain
bool AnalysisDomain::doesIntersectWithMinusYEdgeOfAnalysisDomain( const CommonParameters::XY& coord ) const{
	return ( fabs( coord.Y - m_minCoordY ) < m_eps && m_minCoordX <= coord.X && coord.X <= m_maxCoordX );
}

// Get flag whether specified point locates on the specified corner of analysis domain
bool AnalysisDomain::doesLocateOnGivenCorner( const CommonParameters::XY& coord, const AnalysisDomain::FourCorner& icorner ) const{

	switch( icorner ){
		case AnalysisDomain::XPLUS_YMINUS:
			return ( fabs( coord.X - m_maxCoordX ) < m_eps && fabs( coord.Y - m_minCoordY ) < m_eps );
			break;
		case AnalysisDomain::XPLUS_YPLUS:
			return ( fabs( coord.X - m_maxCoordX ) < m_eps && fabs( coord.Y - m_maxCoordY ) < m_eps );
			break;
		case AnalysisDomain::XMINUS_YPLUS:
			return ( fabs( coord.X - m_minCoordX ) < m_eps && fabs( coord.Y - m_maxCoordY ) < m_eps );
			break;
		case AnalysisDomain::XMINUS_YMINUS:
			return ( fabs( coord.X - m_minCoordX ) < m_eps && fabs( coord.Y - m_minCoordY ) < m_eps );
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong corner ID !! ID = " << icorner << " !!" << std::endl;
			exit(1);
			break;
	}

}

// Get flag whether specified point locates on a corner of analysis domain
bool AnalysisDomain::doesLocateOnCorner( const CommonParameters::XY& coord ) const{

	//if( ( fabs( coord.X - m_maxCoordX ) < m_eps || fabs( coord.X - m_minCoordX ) < m_eps ) &&
	//	( fabs( coord.Y - m_maxCoordY ) < m_eps || fabs( coord.Y - m_minCoordY ) < m_eps ) ){
	//	return true;
	//}

	//return false;

	if( doesLocateOnGivenCorner( coord, AnalysisDomain::XPLUS_YMINUS ) ||
		doesLocateOnGivenCorner( coord, AnalysisDomain::XPLUS_YPLUS ) ||
		doesLocateOnGivenCorner( coord, AnalysisDomain::XMINUS_YPLUS ) ||
		doesLocateOnGivenCorner( coord, AnalysisDomain::XMINUS_YMINUS ) ){
		return true;
	}

	return false;

}

// Get flag whether specified point locate within the analysis domain
bool AnalysisDomain::doesLocateWithinAnalysisDomain( const CommonParameters::XY& coord ) const{
	
	if( m_minCoordX <= coord.X && coord.X <= m_maxCoordX && m_minCoordY <= coord.Y && coord.Y <= m_maxCoordY ){
		return true;
	}

	return false;

}

// Calculate length of X direction
double AnalysisDomain::calcXLength() const{
	return fabs( m_maxCoordX - m_minCoordX );
}

// Calculate length of Y direction
double AnalysisDomain::calcYLength() const{
	return fabs( m_maxCoordY - m_minCoordY );
}

// Calculate length of Z direction
double AnalysisDomain::calcZLength() const{
	return fabs( m_maxCoordZ - m_minCoordZ );
}

// Calculate X coordinate of center
double AnalysisDomain::calcCenterCoordX() const{
	return 0.5 * ( m_maxCoordX + m_minCoordX );
}

// Calculate Y coordinate of center
double AnalysisDomain::calcCenterCoordY() const{
	return 0.5 * ( m_maxCoordY + m_minCoordY );
}

// Calculate Z coordinate of center
double AnalysisDomain::calcCenterCoordZ() const{
	return 0.5 * ( m_maxCoordZ + m_minCoordZ );
}

//
//// Get coordinate of the intersection point with the boundary of the analysis domain
//CommonParameters::XY AnalysisDomain::getCoordOfIntersectionPoint( const CommonParameters::XY& coord1, const CommonParameters::XY& coord2 ) const{
//
//	const CommonParameters::XY coordUpperLeft  = { m_maxCoordX, m_minCoordY };
//	const CommonParameters::XY coordUpperRight = { m_maxCoordX, m_maxCoordY };
//	const CommonParameters::XY coordLowerRight = { m_minCoordX, m_maxCoordY };
//	const CommonParameters::XY coordLowerLeft  = { m_minCoordX, m_minCoordY };
//
//	if( intersectTwoSegments( coord1, coord2, coordUpperLeft, coordUpperRight ) ){
//		return calcCoordOfIntersectionPointOfTwoLines( coord1, coord2, coordUpperLeft,  coordUpperRight );
//	}
//	else if( intersectTwoSegments( coord1, coord2, coordUpperRight, coordLowerRight ) ){
//		return calcCoordOfIntersectionPointOfTwoLines( coord1, coord2, coordUpperRight, coordLowerRight );
//	}
//	else if( intersectTwoSegments( coord1, coord2, coordLowerRight, coordLowerLeft ) ){
//		return calcCoordOfIntersectionPointOfTwoLines( coord1, coord2, coordLowerRight, coordLowerLeft  );
//	}
//	else if( intersectTwoSegments( coord1, coord2, coordLowerLeft, coordUpperLeft ) ){
//		return calcCoordOfIntersectionPointOfTwoLines( coord1, coord2, coordLowerLeft,  coordUpperLeft  );
//	}
//
//	OutputFiles::m_logFile << "Error : Specified vector ( " << coord1.X << " , " << coord1.Y << " ) -> ( " << coord2.X << " , " << coord2.Y << " ) does't intersect with boundary of the analysis domain." << std::endl;
//	exit(1);
//
//	CommonParameters::XY val = { 0.0, 0.0 };
//	return val;
//
//}
//
//// Function determine if two segments intersect or not
//bool AnalysisDomain::intersectTwoSegments( const CommonParameters::XY& startPointOf1stSegment, const CommonParameters::XY& endPointOf1stSegment,
//	const CommonParameters::XY& startPointOf2ndSegment, const CommonParameters::XY& endPointOf2ndSegment ) const{
//
//	const double val1 = ( endPointOf1stSegment.Y - startPointOf1stSegment.Y ) * ( startPointOf2ndSegment.X - startPointOf1stSegment.X ) + ( endPointOf1stSegment.X - startPointOf1stSegment.X ) * ( startPointOf1stSegment.Y - startPointOf2ndSegment.Y );
//	const double val2 = ( endPointOf1stSegment.Y - startPointOf1stSegment.Y ) * (   endPointOf2ndSegment.X - startPointOf1stSegment.X ) + ( endPointOf1stSegment.X - startPointOf1stSegment.X ) * ( startPointOf1stSegment.Y -   endPointOf2ndSegment.Y );
//
//	if( val1*val2 <= 0.0 ){
//
//		const double val3 = ( endPointOf2ndSegment.Y - startPointOf2ndSegment.Y ) * ( startPointOf1stSegment.X - startPointOf2ndSegment.X ) + ( endPointOf2ndSegment.X - startPointOf2ndSegment.X ) * ( startPointOf2ndSegment.Y - startPointOf1stSegment.Y );
//		const double val4 = ( endPointOf2ndSegment.Y - startPointOf2ndSegment.Y ) * (   endPointOf1stSegment.X - startPointOf2ndSegment.X ) + ( endPointOf2ndSegment.X - startPointOf2ndSegment.X ) * ( startPointOf2ndSegment.Y -   endPointOf1stSegment.Y );
//
//		if( fabs(val1*val2) < m_eps && fabs(val3*val4) < m_eps ){
//			return false;
//		}else if( val3*val4 <= 0.0 ){
//			return true;
//		}
//
//	}
//
//	return false;
//
//}
//
//// Calculate coordinates of intersection point of two lines
//CommonParameters::XY AnalysisDomain::calcCoordOfIntersectionPointOfTwoLines( const CommonParameters::XY& coord1stLine1, const CommonParameters::XY& coord1stLine2,
//	const CommonParameters::XY& coord2ndLine1, const CommonParameters::XY& coord2ndLine2 ) const{
//
//	const double temp1 = ( coord1stLine2.Y - coord1stLine1.Y )*( coord2ndLine2.X - coord2ndLine1.X ) - ( coord1stLine2.X - coord1stLine1.X )*( coord2ndLine2.Y - coord2ndLine1.Y );
//
//	if( fabs( temp1 ) < m_eps ){
//		OutputFiles::m_logFile << " Error : Divide by zero in calculating X coordinate of intersection point of two lines !!" << std::endl;
//		exit(1);
//	}
//
//	CommonParameters::XY coordIntersectionPoint = { 0.0, 0.0 };
//
//	coordIntersectionPoint.X = ( coord1stLine2.Y - coord1stLine1.Y )*( coord2ndLine2.X - coord2ndLine1.X )*coord1stLine1.X
//							 - ( coord1stLine2.X - coord1stLine1.X )*( coord2ndLine2.Y - coord2ndLine1.Y )*coord2ndLine1.X  
//							 + ( coord1stLine2.X - coord1stLine1.X )*( coord2ndLine2.X - coord2ndLine1.X )*( coord2ndLine1.Y - coord1stLine1.Y );
//	coordIntersectionPoint.X /= temp1;
//
//	const double temp2 = coord1stLine2.X - coord1stLine1.X;
//	const double temp3 = coord2ndLine2.X - coord2ndLine1.X;
//
//	if( fabs( temp2 ) < m_eps && fabs( temp3 ) < m_eps ){
//		OutputFiles::m_logFile << " Error : Divide by zero in calculating Y coordinate of intersection point of two lines !!" << std::endl;
//		exit(1);
//	}
//
//	if( fabs( temp2 ) > fabs( temp3 ) ){
//		coordIntersectionPoint.Y = ( coord1stLine2.Y - coord1stLine1.Y )/temp2*( coordIntersectionPoint.X - coord1stLine1.X ) + coord1stLine1.Y;
//	}else{
//		coordIntersectionPoint.Y = ( coord2ndLine2.Y - coord2ndLine1.Y )/temp3*( coordIntersectionPoint.X - coord2ndLine1.X ) + coord2ndLine1.Y;
//	}
//
//	return coordIntersectionPoint;
//
//}
