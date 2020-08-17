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
#ifndef DBLDEF_TOPOGRAPHY_DATA_LIST
#define DBLDEF_TOPOGRAPHY_DATA_LIST

#include "TopographyData.h"
#include <vector>
#include <map>

// Class of the list of topography data
class TopographyDataList{

public:
	// Return the the instance of the class
    static TopographyDataList* getInstance();

	// Read altitude data
	void readAltitudeData( const std::string& fileName );

	// Read sea depth data
	void readSeaDepthData( const std::string& fileName );
		
	// Read lake bottom altitude data
	void readLakeBottomAltitudeData( const std::string& fileName );

	// Interpolate Z coordinate
	bool interpolateZCoord( const int type, const CommonParameters::XY& coord, const double distanceUsedToAvoidTooSmallDenominator, double& zCoord ) const;

private:
	// Altitude data
	TopographyData m_landAltitudeData;

	// Sea depth data
	TopographyData m_seaDepthData;

	// Lake bottom altitude data
	TopographyData m_lakeBottomAltitudeData;


	// Constructer
	TopographyDataList();

	// Destructer
	~TopographyDataList();

	// Copy constructer
	TopographyDataList(const TopographyDataList& rhs);

	// Assignment operator
	TopographyDataList& operator=(TopographyDataList& rhs);

};

#endif