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
#ifndef DBLDEF_TOPOGRAPHY_DATA
#define DBLDEF_TOPOGRAPHY_DATA

#include "TopographyData.h"
#include "CommonParameters.h"
#include <vector>
#include <string>

// Class of topography data
class TopographyData{

public:
	// Constructer
	TopographyData();

	// Destructer
	~TopographyData();

	// Read topography data
	void readTopographyData( const std::string& inFileName );

	// Interpolate Z coordinate
	bool interpolateZCoord( const CommonParameters::XY& coord, const double distanceUsedToAvoidTooSmallDenominator, double& zCoord ) const;

private:
	// Copy constructer
	TopographyData(const TopographyData& rhs);

	// Assignment operator
	TopographyData& operator=(TopographyData& rhs);

	// Coordinates
	std::vector<CommonParameters::XYZ> m_coords; 

	void outputVTK( const std::string& fname ) const;

#ifdef _TOPO_FUNC
	// Calculate z coordinate by a function
	double calcZByFunction( const CommonParameters::XY& coord ) const;
#endif

};

#endif
