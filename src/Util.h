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
#ifndef DBLDEF_UTIL
#define DBLDEF_UTIL

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

#include "CommonParameters.h"

namespace Util{

// Get coordinate of the intersection point with the boundary of the analysis domain
CommonParameters::XY getCoordOfIntersectionPoint( const CommonParameters::XY& coord1, const CommonParameters::XY& coord2 ); 

// Function determine if two segments intersect or not
bool intersectTwoSegments( const CommonParameters::XY& startPointOf1stSegment, const CommonParameters::XY& endPointOf1stSegment,
	const CommonParameters::XY& startPointOf2ndSegment, const CommonParameters::XY& endPointOf2ndSegment );

// Calculate coordinates of intersection point of two lines
CommonParameters::XY calcCoordOfIntersectionPointOfTwoLines( const CommonParameters::XY& coord1stLine1, const CommonParameters::XY& coord1stLine2,
	const CommonParameters::XY& coord2ndLine1, const CommonParameters::XY& coord2ndLine2 );

// Add coordinates between the specified two points
void addCoordssBetweenTwoPoins( const CommonParameters::XY& startCoord, const CommonParameters::XY& endCoord, const double maxLength, std::vector<CommonParameters::XY>& coords );

// Coordinate transform
CommonParameters::XY rotateCoordCWR( const CommonParameters::XY& coord, const CommonParameters::XY& centerCoord, const double rotationAngle );

// Calculate outer product
double calcOuterProduct( const CommonParameters::XY& startCoord, const CommonParameters::XY& endCoord, const CommonParameters::XY& coord );

// Get flag whether specified locate right hand side of segment
double locateRightHandSide( const CommonParameters::XY& startCoord, const CommonParameters::XY& endCoord, const CommonParameters::XY& coord );

// Calculate center of gravity
CommonParameters::XY calcCenterOfGravity( const CommonParameters::XY& coord0, const CommonParameters::XY& coord1, const CommonParameters::XY& coord2 );

// Calculate length of edge
double calcEdgeLength( const CommonParameters::XY& coord0, const CommonParameters::XY& coord1 );

// Calculate angle of vertex
double calcAngle( const CommonParameters::XY& coord0, const CommonParameters::XY& coord1, const CommonParameters::XY& coord2 );

// Get flag whether node d belongs to the circumcircule of triangle abc
bool inCircle( const CommonParameters::XY& a, const CommonParameters::XY& b, const CommonParameters::XY& c, const CommonParameters::XY& d );

// Calculate coordinate of 3D model from the coordinate of the plane on which 2D mesh is created
CommonParameters::XYZ calcCoordOf3DModel( const CommonParameters::XY& coord2D, const CommonParameters::Boundary& planeType, const double coordZ = 0.0 );

// Transform XYZ to XY
CommonParameters::XY transformXYZToXY( const CommonParameters::XYZ& coord ); 

// Calculate terms of geometric progression
void calculateTermsOfGeometricProgression( const double firstTerm, const double lastTerm, const double sumOfTerms, std::vector<double>& terms ); 

// Calculate terms of geometric progression except the last one
void calculateTermsOfGeometricProgressionExceptLast( const double firstTerm, const double lastTerm, const double sumOfTerms, std::vector<double>& terms ); 

//// Coordinate transform to the one on XY plane from the one of the plane on which 2D mesh is created
//CommonParameters::XY transform2DCoord( const CommonParameters::XY& coord2D, const CommonParameters::Boundary& planeType );

}

#endif
