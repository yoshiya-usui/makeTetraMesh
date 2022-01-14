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
#include "Ellipsoids.h"
#include "OutputFiles.h"
#include "math.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <assert.h> 
#include <stdlib.h>

// Default constructer
Ellipsoids::Ellipsoids():
	m_rotationAngle(0.0),
	m_numEllipsoids(0),
	m_radius(NULL),
	m_edgeLength(NULL),
	m_oblatenessHorizontal(NULL)
{

	m_centerCoord.X = 0.0;
	m_centerCoord.Y = 0.0;

	m_oblateness[0] = NULL;
	m_oblateness[1] = NULL;
}

// Destructer
Ellipsoids::~Ellipsoids(){

	if( m_radius != NULL ){
		delete [] m_radius;
		m_radius = NULL;
	}

	if( m_edgeLength != NULL ){
		delete [] m_edgeLength;
		m_edgeLength = NULL;
	}

	if( m_oblatenessHorizontal != NULL ){
		delete [] m_oblatenessHorizontal;
		m_oblatenessHorizontal = NULL;
	}

	if( m_oblateness[0] != NULL ){
		delete [] m_oblateness[0];
		m_oblateness[0] = NULL;
	}

	if( m_oblateness[1] != NULL ){
		delete [] m_oblateness[1];
		m_oblateness[1] = NULL;
	}

}

// Read control parameters
//void Ellipsoids::readParameters( std::ifstream& ifs ){
void Ellipsoids::readParameters( const std::string& fileName ){

	OutputFiles::m_logFile << "# Read data of ellipsoids used for specifing maximum edge length." << std::endl;

	std::ifstream ifs( fileName.c_str() );

	if( !ifs.is_open() ){
		OutputFiles::m_logFile << "Error : File open error !!" << std::endl;
		exit(1);
	}

	std::string line;
	while( getline( ifs, line ) ){

		//std::string line;
		//ifs >> line;
		//getline( ifs, line );

		if( line.find("CENTER") != std::string::npos ){
			ifs >> m_centerCoord.X >> m_centerCoord.Y >> m_centerCoord.Z;
		}
		else if( line.find("ROTATION") != std::string::npos  ){
			ifs >> m_rotationAngle;
		}
		else if( line.find("ELLIPSOIDS") != std::string::npos ){
			ifs >> m_numEllipsoids;

			if( m_numEllipsoids < 1 ){
				OutputFiles::m_logFile << "Error : Total number of ellipsoids ( " <<  m_numEllipsoids <<" > is less than 1 !!" << std::endl;
				exit(1);
			}

			m_radius = new double[m_numEllipsoids];
			m_edgeLength = new double[m_numEllipsoids];
			m_oblatenessHorizontal = new double[m_numEllipsoids];
			m_oblateness[Ellipsoids::EARTH] = new double[m_numEllipsoids];
			m_oblateness[Ellipsoids::AIR] = new double[m_numEllipsoids];
			for( int i = 0; i < m_numEllipsoids; ++i ){
				ifs >> m_radius[i] >> m_edgeLength[i] >> m_oblatenessHorizontal[i] >> m_oblateness[Ellipsoids::EARTH][i] >> m_oblateness[Ellipsoids::AIR][i];
			}
			for( int i = 1; i < m_numEllipsoids; ++i ){
				if( m_radius[i] < m_radius[i-1] ){
					std::cerr << "Radius of the region " << i << " is smaller than that of the previous region." << std::endl;
					exit(1);
				}
				if( m_edgeLength[i] < m_edgeLength[i-1] ){
					std::cerr << "Edge length of the region " << i << " is smaller than that of the previous region." << std::endl;
					exit(1);
				}
				if( m_oblatenessHorizontal[i] < 0 || m_oblatenessHorizontal[i] > 1 ){
					std::cerr << "Oblateness of horizontal ellipsoid must be smaller than 1 and larger than 0." << std::endl;
					exit(1);
				}
				if( m_oblateness[Ellipsoids::EARTH][i] < 0 || m_oblateness[Ellipsoids::EARTH][i] > 1 ){
					std::cerr << "Oblateness must be smaller than 1 and larger than 0." << std::endl;
					exit(1);
				}
				if( m_oblateness[Ellipsoids::AIR][i] < 0 || m_oblateness[Ellipsoids::AIR][i] > 1 ){
					std::cerr << "Oblateness must be smaller than 1 and larger than 0." << std::endl;
					exit(1);
				}
				//if( m_radius[i]*(1.0-m_oblatenessHorizontal[i]) < m_radius[i-1]*(1.0-m_oblatenessHorizontal[i-1] ) ){
				//	std::cerr << "Length of shorter axis of horizontal ellipsoid " << i << " is less than that of the previous ellipsoid." << std::endl;
				//	exit(1);
				//}
				if( m_radius[i]*(1.0-m_oblateness[Ellipsoids::EARTH][i]) < m_radius[i-1]*(1.0-m_oblateness[Ellipsoids::EARTH][i-1] ) ){
					std::cerr << "Depth of sphere " << i << " is shallower than that of the previous sphere in the earth." << std::endl;
					exit(1);
				}
				if( m_radius[i]*(1.0-m_oblateness[Ellipsoids::AIR][i]) < m_radius[i-1]*(1.0-m_oblateness[Ellipsoids::AIR][i-1] ) ){
					std::cerr << "Depth of sphere " << i << " is shallower than that of the previous sphere in the air." << std::endl;
					exit(1);
				}
			}
		}
		else if( line.find("END") != std::string::npos ){
			break;
		}
		
	}

	OutputFiles::m_logFile << "# Center coordinte [km] : " << m_centerCoord.X << " " << m_centerCoord.Y  << " " << m_centerCoord.Z  << std::endl;
	OutputFiles::m_logFile << "# Rotation angle [deg] : " << m_rotationAngle << std::endl;
	OutputFiles::m_logFile << "# Total number of ellipsoids : " << m_numEllipsoids << std::endl;
	OutputFiles::m_logFile << "<Radius[km]> <Edge Length[km]> <OblatenessHorizontal> <OblatenessVerticalEarth> <OblatenessVerticalAir>" << std::endl;
	OutputFiles::m_logFile.precision(6);
	for( int i = 0; i < m_numEllipsoids; ++i ){
		OutputFiles::m_logFile << std::setw(15) << std::scientific << m_radius[i];
		OutputFiles::m_logFile << std::setw(15) << std::scientific << m_edgeLength[i];
		OutputFiles::m_logFile << std::setw(15) << std::scientific << m_oblatenessHorizontal[i];
		OutputFiles::m_logFile << std::setw(15) << std::scientific << m_oblateness[Ellipsoids::EARTH][i];
		OutputFiles::m_logFile << std::setw(15) << std::scientific << m_oblateness[Ellipsoids::AIR][i] << std::endl;
	}

	// Degrees => Radians
	m_rotationAngle *= CommonParameters::DEG2RAD;

	ifs.close();

}

// Calculate maximum edge length of the specified point
double Ellipsoids::calcMaximumEdgeLength( const CommonParameters::XYZ& coord ) const{

	double minVal = 1.0e20;

	if( locateInEllipsoid( coord, 0 ) ){
		return std::min( m_edgeLength[0], minVal );
	}
	else if( !locateInEllipsoid( coord, m_numEllipsoids-1 ) ){
		return std::min( m_edgeLength[m_numEllipsoids-1], minVal );
	}

	for( int iSphere = 1; iSphere < m_numEllipsoids; ++iSphere ){
		if( locateInEllipsoid( coord, iSphere ) ){
			return std::min( calcEdgeLengthTransitionRegion( coord, iSphere ), minVal );
			break;
		}
	}

	std::cerr << "Wrong!!" << std::endl;
	exit(1);

	return -1.0;

}

// Check whether specified point is located in the Ellipsoid
bool Ellipsoids::locateInEllipsoid( const CommonParameters::XYZ& coord, const int iEllipsoid ) const{
	
	if( iEllipsoid < 0 || iEllipsoid >= m_numEllipsoids ){
		std::cerr << "Wrong shepre ID:  " << iEllipsoid << std::endl;
		exit(1);
	}

	const double vecXOrg = coord.X - m_centerCoord.X;
	const double vecYOrg = coord.Y - m_centerCoord.Y;
	// Coordinate transform
	const double vecX = vecXOrg * cos( - m_rotationAngle ) - vecYOrg * sin( - m_rotationAngle );
	const double vecY = vecXOrg * sin( - m_rotationAngle ) + vecYOrg * cos( - m_rotationAngle );
	const double vecZ = coord.Z - m_centerCoord.Z;

	const int iType = vecZ < 0 ? Ellipsoids::AIR : Ellipsoids::EARTH;

	const double longAxisLength = m_radius[iEllipsoid];
	const double shortAxisLength = longAxisLength * ( 1.0 - m_oblatenessHorizontal[iEllipsoid] );
	const double depth = longAxisLength * ( 1.0 - m_oblateness[iType][iEllipsoid] );

	double val = pow( vecX / longAxisLength, 2 )
			   + pow( vecY / shortAxisLength, 2 )
			   + pow( vecZ / depth, 2 );

	if( val <= 1.0 ){
		return true;
	}

	return false;

}

// Calculate maximum edge length between two spheres
double Ellipsoids::calcEdgeLengthTransitionRegion( const CommonParameters::XYZ& coord, const int iEllipsoid ) const{

	if( iEllipsoid < 1 || iEllipsoid >= m_numEllipsoids ){
		std::cerr << "Wrong shepre ID:  " << iEllipsoid << std::endl;
		exit(1);
	}

	const double vecXOrg = coord.X - m_centerCoord.X;
	const double vecYOrg = coord.Y - m_centerCoord.Y;
	// Coordinate transform
	const double vecX = vecXOrg * cos( - m_rotationAngle ) - vecYOrg * sin( - m_rotationAngle );
	const double vecY = vecXOrg * sin( - m_rotationAngle ) + vecYOrg * cos( - m_rotationAngle );
	const double vecZ = coord.Z - m_centerCoord.Z;

	const int iType = vecZ < 0 ? Ellipsoids::AIR : Ellipsoids::EARTH;

	const double angleHorizontal = atan2( vecY, vecX );
	const double lengthHorizontal = hypot( vecY, vecX );
	const double angleVertical = atan2( vecZ, lengthHorizontal );
	const double length = hypot( lengthHorizontal, vecZ );

	const double length0 = calculateLengthOnEllipsoid( angleHorizontal, angleVertical, iEllipsoid-1, iType );
	const double length1 = calculateLengthOnEllipsoid( angleHorizontal, angleVertical, iEllipsoid, iType );

	if( length < length0 || length > length1 ){
		std::cerr << "Specified coordinate ( " << coord.X << ", " << coord.Y << ", " << coord.Z << " ) does't locate in the trasition region" << std::endl;
		std::cerr << "iEllipsoid : " << iEllipsoid << std::endl;

		std::cerr << "vecX : " << vecX << std::endl;
		std::cerr << "vecY : " << vecY << std::endl;
		std::cerr << "vecZ : " << vecZ << std::endl;

		std::cerr << "lengthHorizontal      : " << lengthHorizontal << std::endl;
		std::cerr << "angleH : " << angleHorizontal << std::endl;
		std::cerr << "angleV : " << angleVertical << std::endl;

		std::cerr << "length  : " << length << std::endl;
		std::cerr << "length0 : " << length0 << std::endl;
		std::cerr << "length1 : " << length1 << std::endl;
		exit(1);
	}

	return ( length - length0 ) / ( length1 - length0 ) * ( m_edgeLength[iEllipsoid] - m_edgeLength[iEllipsoid-1] ) + m_edgeLength[iEllipsoid-1];

}

// Calculate length on ellipsoid
double Ellipsoids::calculateLengthOnEllipsoid( const double angleH, const double angleV, const int iEllipsoid, const int iType ) const{

	if( iEllipsoid < 0 || iEllipsoid >= m_numEllipsoids ){
		std::cerr << "Wrong shepre ID:  " << iEllipsoid << std::endl;
		exit(1);
	}

	if( angleH < -CommonParameters::PI || angleH > CommonParameters::PI ){
		std::cerr << "Horizontal angle is improper : " << angleH << std::endl;
		exit(1);
	}

	if( angleV < -CommonParameters::PI || angleV > CommonParameters::PI ){
		std::cerr << "Vertical angle is improper : " << angleV << std::endl;
		exit(1);
	}

	const double longAxisLength = m_radius[iEllipsoid];
	const double shortAxisLength = longAxisLength * ( 1.0 - m_oblatenessHorizontal[iEllipsoid] );
	const double verticalLength = longAxisLength * ( 1.0 - m_oblateness[iType][iEllipsoid] );

	const double eps = 1.0e-9;
	double lengthH(-1.0);
	if( fabs(angleH - CommonParameters::PI*0.5) < eps || fabs(angleH + CommonParameters::PI*0.5) < eps ){
		lengthH = shortAxisLength;
	}
	else{
		const double constValH = longAxisLength * shortAxisLength / hypot( shortAxisLength, longAxisLength * tan(angleH) );
		lengthH = hypot( constValH, constValH * tan(angleH) );
	}
	
	if( fabs(angleV - CommonParameters::PI*0.5) < eps || fabs(angleV + CommonParameters::PI*0.5) < eps ){
		return verticalLength;
	}
	else{
		const double constValV = lengthH * verticalLength / hypot( verticalLength, lengthH * tan(angleV) );
		return hypot( constValV, constValV * tan(angleV) );
	}

}

// Get total number of ellipsoids
int Ellipsoids::getNumOfEllipsoids() const{

	return m_numEllipsoids;

}

// Get maximum edge length within ellipsoids
double Ellipsoids::getMaxEdgeLength( const int iEllipsoid ) const{

	assert( iEllipsoid >= 0 && iEllipsoid < m_numEllipsoids );

	return m_edgeLength[iEllipsoid];

}
