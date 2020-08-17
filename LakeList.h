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
#ifndef DBLDEF_LAKE_LIST
#define DBLDEF_LAKE_LIST

#include "CommonParameters.h"
#include <vector>
#include <map>

// Class of the list of lake
class LakeList{

public:
	// Return the the instance of the class
    static LakeList* getInstance();

	// Read data of lakes from input file
	void readLakeData();

	// Get total number of lake surface points
	int getNumLakeSurfacePoints() const;

	// Get coordinate of lake surface points
	CommonParameters::XY getPointCoord( const int index ) const;

	// Get lake height
	double getLakeHeight( const int index ) const;

private:
	// Constructer
	LakeList();

	// Destructer
	~LakeList();

	// List of the coordintes of lake surface points
	std::vector<CommonParameters::XYZ> m_coordLakeSurfacePoints;

};

#endif