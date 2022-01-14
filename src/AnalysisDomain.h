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
#ifndef DBLDEF_ANALYSIS_DOMAIN
#define DBLDEF_ANALYSIS_DOMAIN

#include "CommonParameters.h"

// Class of the analysis domain
class AnalysisDomain{

public:

	enum EdgeID{
		BAD_DATA = -1,
		PLUS_X = 0, 
		PLUS_Y, 
		MINUS_X, 
		MINUS_Y, 
	};

	enum FourCorner{
		UNDEFINED_CORNER = -1,
		XPLUS_YMINUS = 0,
		XPLUS_YPLUS,
		XMINUS_YPLUS,
		XMINUS_YMINUS,
	};

	// Return the the instance of the class
    static AnalysisDomain* getInstance();

	// Read data of analysis domain from input file
	void readAnalysisDomainData();

	// Get minimum value of X coordinate
	double getMinCoordX() const;

	// Get maximum value of X coordinate
	double getMaxCoordX() const;

	// Get minimum value of Y coordinate
	double getMinCoordY() const;

	// Get maximum value of Y coordinate
	double getMaxCoordY() const;

	// Get minimum value of Z coordinate
	double getMinCoordZ() const;

	// Get maximum value of Z coordinate
	double getMaxCoordZ() const;

	// Get edge ID of the analysis domain where specified point locate
	int getEdgeID( const CommonParameters::XY& coord ) const;

	// Get flag specifing whether the direction from start point to end point is colock wise on the boundary of analysis domain
	bool isCWROnBoundary( const CommonParameters::XY& startCoord, const CommonParameters::XY& endCoord ) const;

	// Get start coordinate of the edge of the analysis domain
	CommonParameters::XY getStartCoord( const int edgeID ) const;

	// Get end coordinate of the edge of the analysis domain
	CommonParameters::XY getEndCoord( const int edgeID ) const;

	// Get ID of next edge in colock wise order
	int getNextEdgeIDCWR( const int edgeID ) const;

	// Get ID of next edge in anticolock wise order
	int getNextEdgeIDAntiCWR( const int edgeID ) const;

	// Calculate distance on the boundary of analysis domain from upper left point
	double calcDistanceOnBoundaryFromUpperLeft( const CommonParameters::XY& coord ) const;

	// Get the coordinate of the next corner when traveling in the clockwise direction
	CommonParameters::XY getNextCoordCWRDirection( const CommonParameters::XY& coord ) const;

	// Get the coordinate of the next corner when traveling in the anticlockwise direction
	CommonParameters::XY getNextCoordACWRDirection( const CommonParameters::XY& coord ) const;

	//// Get coordinate of the intersection point with the boundary of the analysis domain
	//CommonParameters::XY getCoordOfIntersectionPoint( const CommonParameters::XY& coord1, const CommonParameters::XY& coord2 ) const; 

	// Get flag whether specified point intersects with the boundary of the analysis domain
	bool doesIntersectWithBoundary( const CommonParameters::XY& coord ) const;

	// Get flag whether specified point intersects with the edge of + X side of the analysis domain
	bool doesIntersectWithPlusXEdgeOfAnalysisDomain( const CommonParameters::XY& coord ) const;

	// Get flag whether specified point intersects with the edge of - X side of the analysis domain
	bool doesIntersectWithMinusXEdgeOfAnalysisDomain( const CommonParameters::XY& coord ) const;

	// Get flag whether specified point intersects with the edge of + Y side of the analysis domain
	bool doesIntersectWithPlusYEdgeOfAnalysisDomain( const CommonParameters::XY& coord ) const;

	// Get flag whether specified point intersects with the edge of - Y side of the analysis domain
	bool doesIntersectWithMinusYEdgeOfAnalysisDomain( const CommonParameters::XY& coord ) const;

	// Get flag whether specified point locates on the specified corner of analysis domain
	bool doesLocateOnGivenCorner( const CommonParameters::XY& coord, const AnalysisDomain::FourCorner& icorner ) const;

	// Get flag whether specified point locates on a corner of analysis domain
	bool doesLocateOnCorner( const CommonParameters::XY& coord ) const;

	// Get flag whether specified point locate within the analysis domain
	bool doesLocateWithinAnalysisDomain( const CommonParameters::XY& coord ) const;

	// Calculate length of X direction
	double calcXLength() const;

	// Calculate length of Y direction
	double calcYLength() const;

	// Calculate length of Z direction
	double calcZLength() const;

	// Calculate X coordinate of center
	double calcCenterCoordX() const;

	// Calculate Y coordinate of center
	double calcCenterCoordY() const;

	// Calculate Z coordinate of center
	double calcCenterCoordZ() const;

private:
	// Constructer
	AnalysisDomain();

	// Destructer
	~AnalysisDomain();

	// Copy constructer
	AnalysisDomain(const AnalysisDomain& rhs);

	// Assignment operator
	AnalysisDomain& operator=(const AnalysisDomain& rhs);

	// Minimum value of X coordinate
	double m_minCoordX;

	// Maximum value of X coordinate
	double m_maxCoordX;

	// Minimum value of Y coordinate
	double m_minCoordY;

	// Maximum value of Y coordinate
	double m_maxCoordY;

	// Minimum value of Z coordinate
	double m_minCoordZ;

	// Maximum value of Z coordinate
	double m_maxCoordZ;

	// Threshold value
	static const double m_eps;

	//// Function determine if two segments intersect or not
	//bool intersectTwoSegments( const CommonParameters::XY& startPointOf1stSegment, const CommonParameters::XY& endPointOf1stSegment,
	//	const CommonParameters::XY& startPointOf2ndSegment, const CommonParameters::XY& endPointOf2ndSegment ) const;

	//// Calculate coordinates of intersection point of two lines
	//CommonParameters::XY calcCoordOfIntersectionPointOfTwoLines( const CommonParameters::XY& coord1stLine1, const CommonParameters::XY& coord1stLine2,
	//	const CommonParameters::XY& coord2ndLine1, const CommonParameters::XY& coord2ndLine2 ) const;


};

#endif
