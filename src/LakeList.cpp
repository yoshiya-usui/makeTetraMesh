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