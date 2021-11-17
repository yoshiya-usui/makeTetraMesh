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
#include "TopographyDataList.h"
#include "Node.h"
#include "OutputFiles.h"
#include <stdlib.h> 

// Return the instance of the class
TopographyDataList* TopographyDataList::getInstance(){
   	static TopographyDataList instance;// The only instance
  	return &instance;
}

// Default constructer
TopographyDataList::TopographyDataList()
{
}

// Destructer
TopographyDataList::~TopographyDataList()
{
}

// Read altitude data
void TopographyDataList::readAltitudeData( const std::string& fileName ){
	m_landAltitudeData.readTopographyData( fileName );
}

// Read sea depth data
void TopographyDataList::readSeaDepthData( const std::string& fileName ){
	m_seaDepthData.readTopographyData( fileName );
}

// Read lake bottom altitude data
void TopographyDataList::readLakeBottomAltitudeData( const std::string& fileName ){
	m_lakeBottomAltitudeData.readTopographyData( fileName );
}

// Interpolate Z coordinate
bool TopographyDataList::interpolateZCoord( const int type, const CommonParameters::XY& coord, const double distanceUsedToAvoidTooSmallDenominator, double& zCoord ) const{

	switch( type ){
		case Node::SEA:
			return m_seaDepthData.interpolateZCoord( coord, distanceUsedToAvoidTooSmallDenominator, zCoord );
			break;
		case Node::LAND:// Go through
		case Node::LAKE:
			return m_landAltitudeData.interpolateZCoord( coord, distanceUsedToAvoidTooSmallDenominator, zCoord );
			break;
		default:
			OutputFiles::m_logFile << " Error : Location type of the node (X,Y) = ( " << coord.X << "," << coord.Y <<  " ) is  improper." << std::endl;
			exit(1);
			return false;
			break;
	}

}
