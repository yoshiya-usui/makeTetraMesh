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
#include "math.h"
#include "Util.h"
#include "OutputFiles.h"
#include "AnalysisDomain.h"


// Get coordinate of the intersection point with the boundary of the analysis domain
CommonParameters::XY Util::getCoordOfIntersectionPoint( const CommonParameters::XY& coord1, const CommonParameters::XY& coord2 ) {

	const AnalysisDomain* const ptrAnalysisDomain = AnalysisDomain::getInstance();
	const double xMin = ptrAnalysisDomain->getMinCoordX(); 
	const double xMax = ptrAnalysisDomain->getMaxCoordX(); 
	const double yMin = ptrAnalysisDomain->getMinCoordY(); 
	const double yMax = ptrAnalysisDomain->getMaxCoordY(); 

	const CommonParameters::XY coordUpperLeft  = { xMax, yMin };
	const CommonParameters::XY coordUpperRight = { xMax, yMax };
	const CommonParameters::XY coordLowerRight = { xMin, yMax };
	const CommonParameters::XY coordLowerLeft  = { xMin, yMin };

	if( Util::intersectTwoSegments( coord1, coord2, coordUpperLeft, coordUpperRight ) ){
		return Util::calcCoordOfIntersectionPointOfTwoLines( coord1, coord2, coordUpperLeft,  coordUpperRight );
	}
	else if( Util::intersectTwoSegments( coord1, coord2, coordUpperRight, coordLowerRight ) ){
		return Util::calcCoordOfIntersectionPointOfTwoLines( coord1, coord2, coordUpperRight, coordLowerRight );
	}
	else if( Util::intersectTwoSegments( coord1, coord2, coordLowerRight, coordLowerLeft ) ){
		return Util::calcCoordOfIntersectionPointOfTwoLines( coord1, coord2, coordLowerRight, coordLowerLeft  );
	}
	else if( Util::intersectTwoSegments( coord1, coord2, coordLowerLeft, coordUpperLeft ) ){
		return Util::calcCoordOfIntersectionPointOfTwoLines( coord1, coord2, coordLowerLeft,  coordUpperLeft  );
	}

	OutputFiles::m_logFile << "Error : Specified vector ( " << coord1.X << " , " << coord1.Y << " ) -> ( " << coord2.X << " , " << coord2.Y << " ) does't intersect with boundary of the analysis domain." << std::endl;
	exit(1);

	CommonParameters::XY val = { 0.0, 0.0 };
	return val;

}

// Function determine if two segments intersect or not
bool Util::intersectTwoSegments( const CommonParameters::XY& startPointOf1stSegment, const CommonParameters::XY& endPointOf1stSegment,
	const CommonParameters::XY& startPointOf2ndSegment, const CommonParameters::XY& endPointOf2ndSegment ){

	const double val1 = ( endPointOf1stSegment.Y - startPointOf1stSegment.Y ) * ( startPointOf2ndSegment.X - startPointOf1stSegment.X ) + ( endPointOf1stSegment.X - startPointOf1stSegment.X ) * ( startPointOf1stSegment.Y - startPointOf2ndSegment.Y );
	const double val2 = ( endPointOf1stSegment.Y - startPointOf1stSegment.Y ) * (   endPointOf2ndSegment.X - startPointOf1stSegment.X ) + ( endPointOf1stSegment.X - startPointOf1stSegment.X ) * ( startPointOf1stSegment.Y -   endPointOf2ndSegment.Y );

	const double EPS = 1.0e-12;

	if( val1*val2 <= 0.0 ){

		const double val3 = ( endPointOf2ndSegment.Y - startPointOf2ndSegment.Y ) * ( startPointOf1stSegment.X - startPointOf2ndSegment.X ) + ( endPointOf2ndSegment.X - startPointOf2ndSegment.X ) * ( startPointOf2ndSegment.Y - startPointOf1stSegment.Y );
		const double val4 = ( endPointOf2ndSegment.Y - startPointOf2ndSegment.Y ) * (   endPointOf1stSegment.X - startPointOf2ndSegment.X ) + ( endPointOf2ndSegment.X - startPointOf2ndSegment.X ) * ( startPointOf2ndSegment.Y -   endPointOf1stSegment.Y );

		if( fabs(val1*val2) < EPS && fabs(val3*val4) < EPS ){
			return false;
		}else if( val3*val4 <= 0.0 ){
			return true;
		}

	}

	return false;

}

// Calculate coordinates of intersection point of two lines
CommonParameters::XY Util::calcCoordOfIntersectionPointOfTwoLines( const CommonParameters::XY& coord1stLine1, const CommonParameters::XY& coord1stLine2,
	const CommonParameters::XY& coord2ndLine1, const CommonParameters::XY& coord2ndLine2 ){

	const double temp1 = ( coord1stLine2.Y - coord1stLine1.Y )*( coord2ndLine2.X - coord2ndLine1.X ) - ( coord1stLine2.X - coord1stLine1.X )*( coord2ndLine2.Y - coord2ndLine1.Y );

	const double EPS = 1.0e-12;

	if( fabs( temp1 ) < EPS ){
		OutputFiles::m_logFile << " Error : Divide by zero in calculating X coordinate of intersection point of two lines !!" << std::endl;
		exit(1);
	}

	CommonParameters::XY coordIntersectionPoint = { 0.0, 0.0 };

	coordIntersectionPoint.X = ( coord1stLine2.Y - coord1stLine1.Y )*( coord2ndLine2.X - coord2ndLine1.X )*coord1stLine1.X
							 - ( coord1stLine2.X - coord1stLine1.X )*( coord2ndLine2.Y - coord2ndLine1.Y )*coord2ndLine1.X  
							 + ( coord1stLine2.X - coord1stLine1.X )*( coord2ndLine2.X - coord2ndLine1.X )*( coord2ndLine1.Y - coord1stLine1.Y );
	coordIntersectionPoint.X /= temp1;

	const double temp2 = coord1stLine2.X - coord1stLine1.X;
	const double temp3 = coord2ndLine2.X - coord2ndLine1.X;

	if( fabs( temp2 ) < EPS && fabs( temp3 ) < EPS ){
		OutputFiles::m_logFile << " Error : Divide by zero in calculating Y coordinate of intersection point of two lines !!" << std::endl;
		exit(1);
	}

	if( fabs( temp2 ) > fabs( temp3 ) ){
		coordIntersectionPoint.Y = ( coord1stLine2.Y - coord1stLine1.Y )/temp2*( coordIntersectionPoint.X - coord1stLine1.X ) + coord1stLine1.Y;
	}else{
		coordIntersectionPoint.Y = ( coord2ndLine2.Y - coord2ndLine1.Y )/temp3*( coordIntersectionPoint.X - coord2ndLine1.X ) + coord2ndLine1.Y;
	}

	return coordIntersectionPoint;

}

// Add coordinates between the specified two points
void Util::addCoordssBetweenTwoPoins( const CommonParameters::XY& startCoord, const CommonParameters::XY& endCoord,
							   const double maxLength, std::vector<CommonParameters::XY>& coords ){

	const double distance = hypot( endCoord.X - startCoord.X, endCoord.Y - startCoord.Y );

	if( distance <= maxLength ){
		return;
	}

	const int numDiv = static_cast<int>(distance/maxLength) + 1;

	const double incX = ( endCoord.X - startCoord.X ) / static_cast<double>( numDiv );
	const double incY = ( endCoord.Y - startCoord.Y ) / static_cast<double>( numDiv );

	for( int i = 1; i < numDiv; ++i ){
		CommonParameters::XY coord = { startCoord.X + incX * static_cast<double>(i), startCoord.Y + incY * static_cast<double>(i) };
		coords.push_back( coord );
	}

}

// Coordinate transform in a clockwise manner
CommonParameters::XY Util::rotateCoordCWR( const CommonParameters::XY& coord, const CommonParameters::XY& centerCoord, const double rotationAngle ){

	const double vecXOrg = coord.X - centerCoord.X;
	const double vecYOrg = coord.Y - centerCoord.Y;

	// Coordinate transform
	const double vecX = vecXOrg * cos( rotationAngle ) - vecYOrg * sin( rotationAngle );
	const double vecY = vecXOrg * sin( rotationAngle ) + vecYOrg * cos( rotationAngle );

	CommonParameters::XY coordRotated = { vecX + centerCoord.X, vecY + centerCoord.Y };

	return coordRotated;

}

// Calculate outer product
double Util::calcOuterProduct( const CommonParameters::XY& startCoord, const CommonParameters::XY& endCoord, const CommonParameters::XY& coord ){

	return ( endCoord.X - startCoord.X )*( coord.Y - startCoord.Y ) - ( endCoord.Y - startCoord.Y )*( coord.X - startCoord.X );

}

// Get flag whether specified locate right hand side of segment
double Util::locateRightHandSide( const CommonParameters::XY& startCoord, const CommonParameters::XY& endCoord, const CommonParameters::XY& coord ){

	if( Util::calcOuterProduct( startCoord, endCoord, coord ) > 0 ){
		return true;
	}

	return false;

}

// Calculate center of gravity
CommonParameters::XY Util::calcCenterOfGravity( const CommonParameters::XY& coord0, const CommonParameters::XY& coord1, const CommonParameters::XY& coord2 ){

	CommonParameters::XY centerCoord = {
		( coord0.X + coord1.X + coord2.X ) / 3.0,
		( coord0.Y + coord1.Y + coord2.Y ) / 3.0
	};

	return centerCoord;

}

// Calculate length of edge
double Util::calcEdgeLength( const CommonParameters::XY& coord0, const CommonParameters::XY& coord1 ){

	return hypot( coord0.X - coord1.X,  coord0.Y - coord1.Y );

}

// Calculate angle of vertex
double Util::calcAngle( const CommonParameters::XY& coord0, const CommonParameters::XY& coord1, const CommonParameters::XY& coord2 ){

	//const double leng0 = Util::calcEdgeLength( coord0, coord1 );
	//const double leng1 = Util::calcEdgeLength( coord1, coord2 );
	//const double leng2 = Util::calcEdgeLength( coord2, coord0 );

	//const double val = ( leng0*leng0 + leng1*leng1 - leng2*leng2 ) / ( 2.0 * leng0 * leng1 );

	const double leng0 = Util::calcEdgeLength( coord0, coord1 );
	const double leng1 = Util::calcEdgeLength( coord1, coord2 );
	const std::pair<double,double> vec0( coord0.X - coord1.X, coord0.Y - coord1.Y );
	const std::pair<double,double> vec1( coord2.X - coord1.X, coord2.Y - coord1.Y );

	const double val = ( vec0.first * vec1.first + vec0.second * vec1.second ) / ( leng0 * leng1 );

	return acos( val );

}

// Get flag whether node d belongs to the circumcircule of triangle abc
bool Util::inCircle( const CommonParameters::XY& a, const CommonParameters::XY& b, const CommonParameters::XY& c, const CommonParameters::XY& d ){

	const double row0[3] = {
		a.X - d.X,
		a.Y - d.Y,
		pow( a.X - d.X, 2.0 ) + pow( a.Y - d.Y, 2.0 )
	};

	const double row1[3] = {
		b.X - d.X,
		b.Y - d.Y,
		pow( b.X - d.X, 2.0 ) + pow( b.Y - d.Y, 2.0 )
	};

	const double row2[3] = {
		c.X - d.X,
		c.Y - d.Y,
		pow( c.X - d.X, 2.0 ) + pow( c.Y - d.Y, 2.0 )
	};

	const double val = row0[0] * row1[1] * row2[2] + row0[1] * row1[2] * row2[0] + row0[2] * row1[0] * row2[1]
					 - row0[2] * row1[1] * row2[0] - row0[1] * row1[0] * row2[2] - row0[0] * row1[2] * row2[1];

	if( val > 0.0 ){
		return true;
	}

	return false;

}

// Calculate coordinate of 3D model from the coordinate of the plane on which 2D mesh is created
CommonParameters::XYZ Util::calcCoordOf3DModel( const CommonParameters::XY& coord2D, const CommonParameters::Boundary& planeType, const double coordZ ){

	CommonParameters::XYZ coord3D = { coord2D.X, coord2D.Y, 0.0 };

	const AnalysisDomain* const ptrAnalysisDomain = AnalysisDomain::getInstance();

	switch(planeType){
		case CommonParameters::SURFACE:
			break;
		case CommonParameters::TOP:
			coord3D.Z = ptrAnalysisDomain->getMinCoordZ();
			break;
		case CommonParameters::BOT:
			coord3D.Z = ptrAnalysisDomain->getMaxCoordZ();
			break;
		case CommonParameters::LAYER:
			coord3D.Z = coordZ;
			break;
		case CommonParameters::YZ_MINUS:
			coord3D.X = ptrAnalysisDomain->getMinCoordX();
			coord3D.Y = coord2D.X;
			coord3D.Z = coord2D.Y;
			break;
		case CommonParameters::YZ_PLUS:
			coord3D.X = ptrAnalysisDomain->getMaxCoordX();
			coord3D.Y = coord2D.X;
			coord3D.Z = coord2D.Y;
			break;
		case CommonParameters::ZX_MINUS:
			coord3D.Z = coord2D.X;
			coord3D.X = coord2D.Y;
			coord3D.Y = ptrAnalysisDomain->getMinCoordY();
			break;
		case CommonParameters::ZX_PLUS:
			coord3D.Z = coord2D.X;
			coord3D.X = coord2D.Y;
			coord3D.Y = ptrAnalysisDomain->getMaxCoordY();
			break;
		default:
			OutputFiles::m_logFile << " Error : Unknown plane type : " << planeType << std::endl;
			exit(1);
			break;
	}

	return coord3D;
}

// Transform XYZ to XY
CommonParameters::XY Util::transformXYZToXY( const CommonParameters::XYZ& coord ){
	const CommonParameters::XY coordXY = { coord.X, coord.Y };
	return coordXY;
}

// Calculate terms of geometric progression
void Util::calculateTermsOfGeometricProgression( const double firstTerm, const double lastTerm, const double sumOfTerms, std::vector<double>& terms ){

	const double EPS = 1.0e-3;

	if( fabs(firstTerm - lastTerm) < EPS ){
		const double minVal = std::min( firstTerm, lastTerm );
		const int numOfTerms = static_cast<int>( sumOfTerms / minVal );
		if( numOfTerms < 1 ){
			terms.push_back( sumOfTerms );
			return;
		}
		double sum(0.0);
		for( int i = 0; i < numOfTerms; ++i ){
			sum += minVal;
			terms.push_back( sum );
		}
	}
	else{
		const double ratio = (sumOfTerms - firstTerm) / (sumOfTerms - lastTerm);
		const int numOfTerms = 1 + static_cast<int>( log(lastTerm/firstTerm) / log(ratio) );
		if( numOfTerms < 1 ){
			terms.push_back( sumOfTerms );
			return;
		}
		double sum(0.0);
		for( int i = 0; i < numOfTerms; ++i ){
			sum += firstTerm * pow( ratio, i );
			terms.push_back( sum );
		}
	}


};

// Calculate terms of geometric progression except the last one
void Util::calculateTermsOfGeometricProgressionExceptLast( const double firstTerm, const double lastTerm, const double sumOfTerms, std::vector<double>& terms ){
	Util::calculateTermsOfGeometricProgression( firstTerm, lastTerm, sumOfTerms, terms );
	if( static_cast<int>( terms.size() ) < 1 ){
		OutputFiles::m_logFile << " Error : Number of terms are less than 1 : " << terms.size() << std::endl;
		exit(1);
	}
	terms.pop_back();
}

//// Coordinate transform to the one on XY plane from the one of the plane on which 2D mesh is created
//CommonParameters::XY Util::transform2DCoord( const CommonParameters::XY& coord2D, const CommonParameters::Boundary& planeType ){
//
//	const CommonParameters::XYZ coord3D = Util::calcCoordOf3DModel(coord2D, planeType);
//	const CommonParameters::XY coordTranformed = { coord3D.X, coord3D.Y };
//	return coordTranformed;
//
//}
