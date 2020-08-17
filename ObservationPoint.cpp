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
#include "ObservationPoint.h"
#include "OutputFiles.h"
#include "Control.h"
#include "math.h"
#include <stdlib.h>
#include <fstream>
#include <iostream>

// Default constructer
ObservationPoint::ObservationPoint():
	m_numCircles(0),
	m_radius(NULL),
	m_maxEdgeLengthWithinCircle(NULL)
{

	m_pointCoord.X = 0.0;
	m_pointCoord.Y = 0.0;

}

// Destructer
ObservationPoint::~ObservationPoint(){

	if( m_radius != NULL ){
		delete [] m_radius;
		m_radius = NULL;
	}
	
	if( m_maxEdgeLengthWithinCircle != NULL ){
		delete [] m_maxEdgeLengthWithinCircle;
		m_maxEdgeLengthWithinCircle = NULL;
	}

}

// Read data of observation point from input file
void ObservationPoint::readObservationPointData( std::ifstream& ifs ){

	ifs >> m_pointCoord.X >> m_pointCoord.Y;
	if( ( Control::getInstance() )->getInvertSignYcoord() ){
		m_pointCoord.Y *= -1.0;
	}

	ifs >> m_numCircles;
	if( m_numCircles <= 0 ){
		OutputFiles::m_logFile << "Error : Total number of circles less than 1. : " << m_numCircles << std::endl;
		exit(1);
	}

	m_radius = new double[m_numCircles];
	m_maxEdgeLengthWithinCircle = new double[m_numCircles];

	for( int i = 0; i < m_numCircles; ++i ){
		ifs >> m_radius[i] >> m_maxEdgeLengthWithinCircle[i];
		if( i > 0 && ( m_radius[i] < m_radius[i-1] ) ){
			OutputFiles::m_logFile << "Error : Inner radius ( " << m_radius[i-1] << " ) is greater than outer radius ( " << m_radius[i] << " ) !! " << std::endl;
			exit(1);
		}
		if( i > 0 && ( m_maxEdgeLengthWithinCircle[i] < m_maxEdgeLengthWithinCircle[i-1] ) ){
			OutputFiles::m_logFile << "Error : Inner edge length ( " << m_maxEdgeLengthWithinCircle[i-1] << " ) is greater than outer edge length ( " << m_maxEdgeLengthWithinCircle[i] << " ) !! " << std::endl;
			exit(1);
		}
	}

	OutputFiles::m_logFile << "# Coordinate : " << m_pointCoord.X << " " << m_pointCoord.Y << std::endl;
	OutputFiles::m_logFile << "# Number of circles : " << m_numCircles << std::endl;
	for( int i = 0; i < m_numCircles; ++i ){
		OutputFiles::m_logFile << "# Radius[km] : " << m_radius[i] << ", edge length[km] : " << m_maxEdgeLengthWithinCircle[i] << std::endl;
	}
	
}

//// Calculate maximum length of specified coordinate
//double ObservationPoint::calcMaximumLengthOfPoint( const CommonParameters::XYZ& coord ) const{
//
//	const double vecX = coord.X - m_pointCoord.X;
//	const double vecY = coord.Y - m_pointCoord.Y;
//	const double vecZ = coord.Z - m_pointCoord.Z;
//
//	const double radius = sqrt( pow( vecX, 2.0 ) + pow( vecY, 2.0 ) + pow( vecZ, 2.0 ) );
//
//	const double maxEdgeLengthOut = ( Control::getInstance() )->calcMaximumLengthFromControlParamOnly( coord );
//
//	if( radius <= m_innerRadius ){
//		return m_maxEdgeLengthWithinCircle;
//	}
//	else if( radius <= m_outerRadius ){
//		return m_maxEdgeLengthWithinCircle + ( maxEdgeLengthOut - m_maxEdgeLengthWithinCircle ) * radius / m_outerRadius;
//	}
//
//	return maxEdgeLengthOut;
//
//}

// Calculate maximum length of specified coordinate on X-Y plane ( z = 0 ) 
double ObservationPoint::calcMaximumLengthOfPoint( const CommonParameters::XY& coord ) const{

	const double vecX = coord.X - m_pointCoord.X;
	const double vecY = coord.Y - m_pointCoord.Y;

	const double radius = sqrt( pow( vecX, 2.0 ) + pow( vecY, 2.0 ) );
	const double edgeLeng0 = ( Control::getInstance() )->calcMaximumLengthFromControlParamOnly( coord, CommonParameters::SURFACE );

	if( radius <= m_radius[0] ){
		return std::min( m_maxEdgeLengthWithinCircle[0], edgeLeng0 );
	}

	double edgeLeng1 = edgeLeng0;
	for( int i = m_numCircles - 1; i >= 1; --i ){
		if( radius <= m_radius[i] ){
			 edgeLeng1 = m_maxEdgeLengthWithinCircle[i-1] + ( radius - m_radius[i-1] ) / ( m_radius[i] - m_radius[i-1] ) * ( m_maxEdgeLengthWithinCircle[i] - m_maxEdgeLengthWithinCircle[i-1] );
		}else{
			break;
		}
	}

	return std::min( edgeLeng1, edgeLeng0 );

}
