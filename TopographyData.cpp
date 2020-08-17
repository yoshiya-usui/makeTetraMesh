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
#include "TopographyData.h"
#include "OutputFiles.h"
#include "Control.h"
#include "math.h"
#include "AnalysisDomain.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>

// Default constructer
TopographyData::TopographyData()
{
}

// Destructer
TopographyData::~TopographyData(){
}

// Read topography data
void TopographyData::readTopographyData( const std::string& inFileName ){

	// Open input file -----
	std::ifstream ifs( inFileName.c_str() );
	if( !ifs ) {
		std::cerr << "Cannot open file " << inFileName << std::endl;
		exit(1);
	}
	OutputFiles::m_logFile << "# Read data from " << inFileName << std::endl;
	//----------------------

	const double maxDist = ( Control::getInstance() )->getMaxDistanceForInterpolating();
	const AnalysisDomain* const ptrAnalysisDomain = AnalysisDomain::getInstance();
	const double xMin = ptrAnalysisDomain->getMinCoordX();
	const double xMax = ptrAnalysisDomain->getMaxCoordX(); 
	const double yMin = ptrAnalysisDomain->getMinCoordY(); 
	const double yMax = ptrAnalysisDomain->getMaxCoordY(); 

	const bool invertYCoord = ( Control::getInstance() )->getInvertSignYcoord();

	std::string sbuf;
	int icount(0);
    while( std::getline(ifs, sbuf) ){

		CommonParameters::XYZ coord = { 0.0, 0.0, 0.0 };

		std::istringstream iss( sbuf );
		iss >> coord.X >> coord.Y >> coord.Z;

		if( coord.X < xMin - maxDist || coord.X > xMax + maxDist || coord.Y < yMin - maxDist || coord.Y > yMax + maxDist ){
			continue;
		}

		if( invertYCoord ){
			coord.Y *= -1.0;
		}

		//coord.Z *= CommonParameters::METER2KILOMETER;
		m_coords.push_back( coord );

		++icount;
    }

	OutputFiles::m_logFile << "# Total number of data is " << icount << std::endl;

	ifs.close();

	//const std::string outputFileVTK = inFileName + ".vtk";

	//outputVTK( outputFileVTK );

}

#ifdef _TOPO_FUNC
// Calculate z coordinate by a function
double TopographyData::calcZByFunction( const CommonParameters::XY& coord ) const{

	const double height = 2.0 - hypot(coord.X, coord.Y) / 30.0;
	if( hypot(coord.X, coord.Y) > 10.0 ){
		return height;
	}
	else{
		const double lakeHeight = 2.0 - 10.0 / 30.0;
		const double lakeBottomAltitude = lakeHeight - ( 10.0 - hypot(coord.X, coord.Y) ) / 20.0;
		return lakeBottomAltitude;
	}

	//----- For Li et al. 2008 >>>>>
	//----- Type1 >>>>>
	/*
	const double maxSeaDepth = 2.0;
	if( coord.X > 0.0 ){
		return 0.0;
	}
	else if( coord.X < -2.9 ){
		return maxSeaDepth;
	}
	else{
		const double x1 = -1.95;
		const double x2 = -0.95;
		if( coord.X <= x1 ){
			const double deg = ( coord.X + 2.90 ) * 180.0 / 3.0;
			const double rad = deg * CommonParameters::DEG2RAD;
			if( rad < 0.0 ){
				return maxSeaDepth;
			}else{
				return maxSeaDepth * 0.5 * ( sin( rad + CommonParameters::PI * 0.5 ) - 1.0 ) + maxSeaDepth;
			}
		}
		else if( coord.X >= -0.05 ){
			const double deg = ( -0.05 + 2.75 ) * 180.0 / 2.8;
			const double rad = deg * CommonParameters::DEG2RAD;
			const double z = maxSeaDepth * 0.5 * ( sin( rad + CommonParameters::PI * 0.5 ) - 1.0 ) + maxSeaDepth;
			return - z / 0.05 * ( coord.X + 0.05 ) + z;
		}
		else if( coord.X >= x2 ){
			const double deg = ( coord.X + 2.75 ) * 180.0 / 2.8;
			const double rad = deg * CommonParameters::DEG2RAD;
			return maxSeaDepth * 0.5 * ( sin( rad + CommonParameters::PI * 0.5 ) - 1.0 ) + maxSeaDepth;
		}
		else{
			const double deg1 = ( x1 + 2.90 ) * 180.0 / 3.0;
			const double rad1 = deg1 * CommonParameters::DEG2RAD;
			const double z1 = maxSeaDepth * 0.5 * ( sin( rad1 + CommonParameters::PI * 0.5 ) - 1.0 ) + maxSeaDepth;
			const double deg2 = ( x2 + 2.75 ) * 180.0 / 2.8;
			const double rad2 = deg2 * CommonParameters::DEG2RAD;
			const double z2 = maxSeaDepth * 0.5 * ( sin( rad2 + CommonParameters::PI * 0.5 ) - 1.0 ) + maxSeaDepth;
			return ( z2 - z1 ) / ( x2 - x1 ) * ( coord.X - x1 ) + z1;
		}
	}
	return -9999.0;
	*/
	//----- Type1 <<<<<

	//----- Type2 >>>>>
	/*
	const double maxSeaDepth = 2.0;
	if( coord.X > 0.0 ){
		return 0.0;
	}
	else if( coord.X < -2.85 ){
		return maxSeaDepth;
	}
	else{
		const double deg = ( coord.X + 2.85 ) * 180.0 / 2.85;
		const double rad = deg * CommonParameters::DEG2RAD;
		if( rad < 0.0 ){
			return maxSeaDepth;
		}else{
			return maxSeaDepth * 0.5 * ( sin( rad + CommonParameters::PI * 0.5 ) - 1.0 ) + maxSeaDepth;
		}
	}
	return -9999.0;
	*/
	//----- Type2 <<<<<

	//----- Type3 >>>>>
	/*
	const double maxSeaDepth = 2.0;
	if( coord.X > 0.0 ){
		return 0.0;
	}
	else if( coord.X < -2.8 ){
		return maxSeaDepth;
	}
	else{
		const double x1 = -1.9;
		const double x2 = -0.9;
		if( coord.X <= x1 ){
			const double deg = ( coord.X + 2.80 ) * 180.0 / 2.75;
			const double rad = deg * CommonParameters::DEG2RAD;
			if( rad < 0.0 ){
				return maxSeaDepth;
			}else{
				return maxSeaDepth * 0.5 * ( sin( rad + CommonParameters::PI * 0.5 ) - 1.0 ) + maxSeaDepth;
			}
		}
		else if( coord.X >= -0.05 ){
			const double deg = ( -0.05 + 2.75 ) * 180.0 / 2.8;
			const double rad = deg * CommonParameters::DEG2RAD;
			const double z = maxSeaDepth * 0.5 * ( sin( rad + CommonParameters::PI * 0.5 ) - 1.0 ) + maxSeaDepth;
			return - z / 0.05 * ( coord.X + 0.05 ) + z;
		}
		else if( coord.X >= x2 ){
			const double deg = ( coord.X + 2.75 ) * 180.0 / 2.8;
			const double rad = deg * CommonParameters::DEG2RAD;
			return maxSeaDepth * 0.5 * ( sin( rad + CommonParameters::PI * 0.5 ) - 1.0 ) + maxSeaDepth;
		}
		else{
			const double deg2 = ( x2 + 2.75 ) * 180.0 / 2.8;
			const double rad2 = deg2 * CommonParameters::DEG2RAD;
			const double z2 = maxSeaDepth * 0.5 * ( sin( rad2 + CommonParameters::PI * 0.5 ) - 1.0 ) + maxSeaDepth;
			const double z1 = z2 + 1;
			return ( z2 - z1 ) / ( x2 - x1 ) * ( coord.X - x1 ) + z1;
		}
	}
	return -9999.0;
	*/
	//----- Type3 <<<<<

	//----- Type4 >>>>>
	//const double maxSeaDepth = 2.0;
	//if( coord.X > 0.0 ){
	//	return 0.0;
	//}
	//else if( coord.X < -2.9 ){
	//	return maxSeaDepth;
	//}
	//else{
	//	const double x1 = -1.9;
	//	const double x2 = -0.9;
	//	if( coord.X <= x1 ){
	//		const double deg = ( coord.X + 2.90 ) * 180.0 / 3.0555;
	//		const double rad = deg * CommonParameters::DEG2RAD;
	//		if( rad < 0.0 ){
	//			return maxSeaDepth;
	//		}else{
	//			return maxSeaDepth * 0.5 * ( sin( rad + CommonParameters::PI * 0.5 ) - 1.0 ) + maxSeaDepth;
	//		}
	//	}
	//	else if( coord.X >= -0.05 ){
	//		const double deg = ( -0.05 + 2.75 ) * 180.0 / 2.8;
	//		const double rad = deg * CommonParameters::DEG2RAD;
	//		const double z = maxSeaDepth * 0.5 * ( sin( rad + CommonParameters::PI * 0.5 ) - 1.0 ) + maxSeaDepth;
	//		return - z / 0.05 * ( coord.X + 0.05 ) + z;
	//	}
	//	else if( coord.X >= x2 ){
	//		const double deg = ( coord.X + 2.75 ) * 180.0 / 2.8;
	//		const double rad = deg * CommonParameters::DEG2RAD;
	//		return maxSeaDepth * 0.5 * ( sin( rad + CommonParameters::PI * 0.5 ) - 1.0 ) + maxSeaDepth;
	//	}
	//	else{
	//		const double deg2 = ( x2 + 2.75 ) * 180.0 / 2.8;
	//		const double rad2 = deg2 * CommonParameters::DEG2RAD;
	//		const double z2 = maxSeaDepth * 0.5 * ( sin( rad2 + CommonParameters::PI * 0.5 ) - 1.0 ) + maxSeaDepth;
	//		const double z1 = z2 + 1;
	//		return ( z2 - z1 ) / ( x2 - x1 ) * ( coord.X - x1 ) + z1;
	//	}
	//}
	//return -9999.0;
	//----- Type4 <<<<<
	//----- For Li et al. 2008 <<<<<

	//----- Not delete for future use >>>>>
	//double seaDepth = 4.0;

	//const int numMount(2);
	//const double sigma[2] = { 50.0, 50.0 };
	//const double maxHeight[2] = { 2.0, 1.6 };
	//const double centerX[2] = { -50.0, 50.0 };
	//const double centerY[2] = { -50.0, 50.0 };
	////const int numMount(1);
	////const double sigma[1] = { 50.0 };
	////const double maxHeight[1] = { 2.0 };
	////const double centerX[1] = { 0.0 };
	////const double centerY[1] = { 0.0 };
	//for( int iMount = 0; iMount < numMount; ++iMount ){
	//	const double factor = maxHeight[iMount] * sqrt( 2.0 * CommonParameters::PI ) * sigma[iMount];
	//	const double dist = hypot( coord.X - centerX[iMount], coord.Y - centerY[iMount] );
	//	const double height = factor / ( sqrt( 2.0 * CommonParameters::PI ) * sigma[iMount] ) * exp( - 0.5 * pow(dist/sigma[iMount], 2) );
	//	seaDepth -= height;
	//}

	////{
	////	const double waveLengthSmall = 5.0;
	////	const double amplitudeSmall = 0.05;
	////	const double radianX = 2.0 * CommonParameters::PI * coord.X / waveLengthSmall;
	////	const double radianY = 2.0 * CommonParameters::PI * coord.Y / waveLengthSmall;
	////	seaDepth -= amplitudeSmall * sin(radianX) * cos(radianY);
	////}	

	////{
	////	const double waveLengthSmall = 20.0;
	////	const double amplitudeSmall = 0.15;
	////	const double radianX = 2.0 * CommonParameters::PI * coord.X / waveLengthSmall;
	////	const double radianY = 2.0 * CommonParameters::PI * coord.Y / waveLengthSmall;
	////	seaDepth -= amplitudeSmall * cos(radianX) * sin(radianY);
	////}
	//
	//return seaDepth;

	////const double sigma = 12.5 * 0.5 * sqrt(2.0);
	////const double maxHeight = 2.5;
	////const double seaDepth  = 4.0;
	////const double factor = maxHeight * sqrt(2.0*CommonParameters::PI) * sigma;

	////const double dist1 = hypot(coord.X - 12.5, coord.Y - 12.5);
	////const double dist2 = hypot(coord.X + 12.5, coord.Y + 12.5);

	////const double height1 = factor / ( sqrt( 2.0 * CommonParameters::PI ) * sigma ) * exp( - 0.5 * pow(dist1/sigma,2) );
	////const double height2 = factor / ( sqrt( 2.0 * CommonParameters::PI ) * sigma ) * exp( - 0.5 * pow(dist2/sigma,2) );

	////const double waveLengthSmall = 3.0;
	////const double amplitudeSmall = 0.3;
	////const double radianX = 2.0 * CommonParameters::PI * coord.X / waveLengthSmall;
	////const double radianY = 2.0 * CommonParameters::PI * coord.Y / waveLengthSmall;
	////const double height3 = amplitudeSmall * sin(radianX) * sin(radianY);

	////return seaDepth - height1 - height2 - height3;
	//----- Not delete for future use <<<<<
}
#endif

// Interpolate Z coordinate by inverse distance weighting
bool TopographyData::interpolateZCoord( const CommonParameters::XY& coord, const double distanceUsedToAvoidTooSmallDenominator, double& zCoord ) const{

#ifdef _TOPO_FUNC
	zCoord = calcZByFunction(coord);
	return true;
#endif

	const Control* const ptrControl = Control::getInstance() ;
	const double maxDist = ptrControl->getMaxDistanceForInterpolating();
	const int maxNumPoint = ptrControl->getMaxNumOfPointsForInterpolating();
	const double factorIDW = ptrControl->getFactorInverseDistanceWeighting();

#ifdef _DEBUG_WRITE
	std::cout << "maxDist : " << maxDist << std::endl;
	std::cout << "maxNumPoint : " << maxNumPoint << std::endl;
#endif

	//double* minDistances = new double[maxNumPoint];
	//int* minIDs = new int[maxNumPoint];
	//for( int i = 0; i < maxNumPoint; ++i ){
	//	minDistances[i] = 1.0e12;
	//	minIDs[i] = -1;
	//}
	std::vector< std::pair<double,int> > stackDistances;
	for( int i = 0; i < maxNumPoint; ++i ){
		stackDistances.push_back( std::make_pair( 1.0e12, -1 ) );// Initialize 
	}

	int iDataMinDist(-1);
	double minDist(1.0e12);

	int iData(0);
	const std::vector<CommonParameters::XYZ>::const_iterator itrEnd = m_coords.end();
	for( std::vector<CommonParameters::XYZ>::const_iterator itr = m_coords.begin(); itr != itrEnd; ++itr, ++iData ){

#ifdef _DEBUG_WRITE
		std::cout << "iData : " << iData << std::endl;
#endif
		const double distance = hypot( itr->X - coord.X, itr->Y - coord.Y );

#ifdef _DEBUG_WRITE
		std::cout << "distance : " << distance << std::endl;
#endif

		if( distance < minDist ){
			iDataMinDist = iData;
			minDist = distance;
		}

		if( distance > maxDist ){
			continue;
		}

		if( distance < stackDistances.back().first ){

			stackDistances.back().first = distance;
			stackDistances.back().second = iData;

			sort( stackDistances.begin(), stackDistances.end() );
		}

	}

	if( stackDistances.front().second < 0 ){
		OutputFiles::m_logFile << " Warning : No topography data were found near the point (X,Y) = ( " << coord.X << ", " << coord.Y << " )" << std::endl;
		if( iDataMinDist < 0 ){
			OutputFiles::m_logFile << " Error : Point serial of the nearest point is negative !!" << std::endl;
			exit(1);
		}
		OutputFiles::m_logFile << "           Thus, the value of the nearest point (distance = " << minDist << ", value = " << m_coords[iDataMinDist].Z << ") is given." << std::endl;
		zCoord = m_coords[iDataMinDist].Z;
		return true;
	}

	zCoord = 0.0;
	double weightSum(0.0);
	for( std::vector< std::pair<double,int> >::iterator itr = stackDistances.begin();
		itr != stackDistances.end(); ++itr ){

		if( itr->second < 0 ){
			continue;
		}

		const double weight = 1.0 / pow( distanceUsedToAvoidTooSmallDenominator + itr->first, factorIDW );
		zCoord += weight * m_coords[itr->second].Z;
		weightSum += weight;

	}
	zCoord /= weightSum;

	return true;

}

void TopographyData::outputVTK( const std::string& fname ) const{

	std::ofstream ofsVTK( fname.c_str() );
	if( !ofsVTK ) {
		std::cerr << "Cannot open file " << fname << std::endl;
		exit(1);
	}

	ofsVTK << "# vtk DataFile Version 2.0" << std::endl;
	ofsVTK << "Altitude" << std::endl;
	ofsVTK << "ASCII" << std::endl;
	ofsVTK << "DATASET UNSTRUCTURED_GRID" << std::endl;

	// output data to vtk file -----
	const int numPoints = static_cast<int>( m_coords.size() );

	ofsVTK.precision(9);
	ofsVTK << std::fixed;

	ofsVTK << "POINTS " << numPoints << " double" << std::endl;
	ofsVTK.precision(6);
	for( std::vector<CommonParameters::XYZ>::const_iterator itr = m_coords.begin(); itr != m_coords.end(); ++itr ){
		ofsVTK << std::setw(15) << std::scientific << itr->X;
		ofsVTK << std::setw(15) << std::scientific << itr->Y;
		ofsVTK << std::setw(15) << std::scientific << itr->Z << std::endl;
	}

	ofsVTK << "CELLS " << numPoints << " " << numPoints*2 << std::endl;
	for( int i = 0; i < numPoints; ++i ){
		ofsVTK << std::setw(10) << 1;
		ofsVTK << std::setw(10) << i << std::endl;
	}

	ofsVTK << "CELL_TYPES " << numPoints << std::endl;
	for( int i = 0; i < numPoints; ++i ){
		ofsVTK << std::setw(10) << 1;
	}

	ofsVTK << "POINT_DATA " << numPoints << std::endl;
	ofsVTK << "SCALARS Height(km) float" << std::endl;
	ofsVTK << "LOOKUP_TABLE default" << std::endl;
	ofsVTK.precision(6);
	for( std::vector<CommonParameters::XYZ>::const_iterator itr = m_coords.begin(); itr != m_coords.end(); ++itr ){
		ofsVTK << std::setw(15) << std::scientific << itr->Z << std::endl;
	}
	//------------------------------

	ofsVTK.close();
}
