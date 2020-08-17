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
