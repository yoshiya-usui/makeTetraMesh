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
#include "LakeList.h"
#include "OutputFiles.h"
#include <assert.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>

// Return the instance of the class
LakeList* LakeList::getInstance(){
   	static LakeList instance;// The only instance
  	return &instance;
}

// Constructer
LakeList::LakeList()
{

}

// Destructer
LakeList::~LakeList()
{
}

// Read data of lakes from input file
void LakeList::readLakeData(){

	std::ifstream ifs( "lake.dat", std::ios::in );

	if( ifs.fail() ){
		return;
	}

	OutputFiles::m_logFile << "# Read coordinates of surface points of each lake from input file." << std::endl;

	int  numLakes(0);
	ifs >> numLakes;

	OutputFiles::m_logFile << "# Total number of lake data : " << numLakes << std::endl;

	if( numLakes < 0 ){
		OutputFiles::m_logFile << "Error : Total number of lake data is negative !!" << std::endl;
		exit(1);
	}

	m_coordLakeSurfacePoints.reserve(numLakes);

	CommonParameters::XYZ coord = { 0.0, 0.0, 0.0 };
	for( int i = 0; i < numLakes; ++i ){
		ifs >> coord.X;
		ifs >> coord.Y;
		ifs >> coord.Z;
		m_coordLakeSurfacePoints.push_back(coord);
    }

	ifs.close();		

}

// Get total number of lake surface points
int LakeList::getNumLakeSurfacePoints() const{

	return static_cast<int>( m_coordLakeSurfacePoints.size() );

}

// Get coordinate of lake surface points
CommonParameters::XY LakeList::getPointCoord( const int index ) const{

	assert( index >= 0 && index < getNumLakeSurfacePoints() );

	const CommonParameters::XY coord = { m_coordLakeSurfacePoints[index].X, m_coordLakeSurfacePoints[index].Y };
	return coord;

}

// Get lake height
double LakeList::getLakeHeight( const int index ) const{

	assert( index >= 0 && index < getNumLakeSurfacePoints() );

	return m_coordLakeSurfacePoints[index].Z;

}