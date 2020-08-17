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
#ifndef DBLDEF_OBSERVING_SITE_LIST
#define DBLDEF_OBSERVING_SITE_LIST

#include "ObservationPoint.h"
#include "ObservationLine.h"
#ifdef _MOD_FOR_NMT
#include "CommonParameters.h"
#endif

// Class of the list of the observing site
class ObservingSiteList{

public:
	// Return the the instance of the class
    static ObservingSiteList* getInstance();

	// Read data of observing site list from input file
	void readObservingSiteData();

	//// Calculate maximum length of specified coordinate
	//double calcMaximumLengthOfPoint( const CommonParameters::XYZ& coord ) const;

	// Calculate maximum length of specified coordinate on X-Y plane ( z = 0 ) 
	double calcMaximumLengthOfPoint( const CommonParameters::XY& coord ) const;

#ifdef _MOD_FOR_NMT
	// Get number of observation line
	int getNumObsLine() const;

	// Get coordinate of start point of an observation line
	CommonParameters::XY getCoordOfStartPointObsLine( const int iobs ) const;

	// Get coordinate of end point of an observation line
	CommonParameters::XY getCoordOfEndPointObsLine( const int iobs ) const;

	// Get maximum edge length of an observation line
	double getMaximumEdgeLengthObsLine( const int iobs ) const;
#endif

private:
	// Default constructer
	ObservingSiteList();

	// Destructer
	~ObservingSiteList();

	// Copy constructer
	ObservingSiteList(const ObservingSiteList& rhs);

	// Assignment operator
	ObservingSiteList& operator=(const ObservingSiteList& rhs);

	// Total number of observation point
	int m_numObservationPoint;

	// Total number of observation line
	int m_numObservationLine;

	// Array of observation point
	ObservationPoint* m_obsPoint;

	// Array of observation line
	ObservationLine* m_obsLine;

};

#endif
