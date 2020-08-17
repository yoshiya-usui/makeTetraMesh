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
#include "TopographyDataList.h"
#include "Node.h"
#include "OutputFiles.h"

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
