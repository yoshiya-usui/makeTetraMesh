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
#include "Control.h"
#include "CoastLineList.h"
#include "ObservingSiteList.h"
#include "AnalysisDomain.h"
#include "OutputFiles.h"
#include "NodeList.h"
#include "BoundaryCurveList.h"
#include "TriangleList.h"
#include "Util.h"
#include "TopographyDataList.h"
#include "LakeList.h"
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iostream>

const double Control::m_eps = 1.0e-12;

// Return the instance of the class
Control* Control::getInstance(){
   	static Control instance;// The only instance
  	return &instance;
}

// Default constructer
Control::Control():
	m_minEdgeLength(0.0),
	m_maxEdgeLength(1.0e9),
	m_thresholdEdgeLengthRatio(3.0),
	m_factorSuperTriangleSize(1.0),
	m_maxIterNumRefiningLargeAngle(0),
	m_maxIterNumRefiningEdgeLength(30),
	m_thresholdAngle(150),
	m_alpha(1.5),
	m_beta(0.2),
	m_maxIterNumCoastLineSmoothing(0),
	m_thresholdAngleCoastLineSmoothing(10.0),
	m_thresholdDistanceCoastLineSmoothing(3.0),
	m_iterNumLaplacianMethod(5),
	m_maxNumOfPointsForInterpolating(3),
	m_maxDistanceForInterpolating(100.0),
	m_distanceUsedToAvoidTooSmallDenominator(1.0e-6),
	m_thresholdDeviationOfHegith(-1.0),
	m_minAltitude(0.0),
	m_maxAltitude(1.0e10),
	m_minSeaDepth(0.01),
	m_maxSeaDepth(1.0e10),
	m_addExtendedRegion(false),
	m_numThreads(1),
	m_lengthOfExtendedRegionPlusX(0.0),
	m_lengthOfExtendedRegionMinusX(0.0),
	m_lengthOfExtendedRegionPlusY(0.0),
	m_lengthOfExtendedRegionMinusY(0.0),
	m_makeSurfMesOnAllBoundary(false),
	m_factorInverseDistanceWeighting(1.0),
	m_exhangeLandAndSeaAtStep1(false),
	m_invertSignYcoord(false),
	m_includeSedimentLayer(false),
	m_thicknessSediment(10.0),
	m_includeLayers(false),
	m_numLayers(-1),
	m_depthLayerInterfaces(NULL)
{
}

// Destructer
Control::~Control(){
}

// Read control parameters from input file
void Control::readControlData(){

	const std::string fileName = "control.dat";

	OutputFiles::m_logFile << "# Read parameters" << std::endl;

	m_ellipsoids.readParameters( fileName );

	std::ifstream ifs( fileName.c_str(), std::ios::in );
	if( ifs.fail() ){
		OutputFiles::m_logFile << "Error : File open error : " << fileName.c_str() << std::endl;
		exit(1);
	}

	bool setOuterEdgeLength = false;
	double dbuf(0.0);

	//while( !ifs.eof() ){
	std::string line;
	while( getline( ifs, line ) ){

		//getline( ifs, line );
		//ifs >> line;

#ifdef _DEBUG_WRITE
		std::cout << "line : " << line << std::endl;
#endif

		if( line.find("MAX_EDGE_RATIO") != std::string::npos ){
			ifs >> m_thresholdEdgeLengthRatio;
		}
		else if( line.find("SUPER_TRIANGLE") != std::string::npos ){
			ifs >> m_factorSuperTriangleSize;
		}
		else if( line.find("ALPHA") != std::string::npos ){
			ifs >> m_alpha;
		}
		else if( line.find("BETA") != std::string::npos ){
			ifs >> m_beta;
		}
		else if( line.find("COAST_SMOOTHING") != std::string::npos ){
			ifs >> m_maxIterNumCoastLineSmoothing;
			ifs >> m_thresholdAngleCoastLineSmoothing;
			ifs >> m_thresholdDistanceCoastLineSmoothing;
		}
		else if( line.find("MAX_ANGLE") != std::string::npos ){
			ifs >> m_thresholdAngle;
		}
		else if( line.find("ITER_MAX") != std::string::npos ){
			ifs >> m_maxIterNumRefiningEdgeLength;
			ifs >> m_maxIterNumRefiningLargeAngle;
		}
		else if( line.find("LAPLACIAN") != std::string::npos ){
			ifs >> m_iterNumLaplacianMethod;
		}
		else if( line.find("INTERPOLATE") != std::string::npos ){
			ifs >> m_maxDistanceForInterpolating;
			ifs >> m_maxNumOfPointsForInterpolating;
			ifs >> m_distanceUsedToAvoidTooSmallDenominator;
		}
		else if( line.find("ALTITUDE") != std::string::npos ){
			ifs >> m_fileNameOfAltitudeData;
			ifs >> m_minAltitude;
			ifs >> m_maxAltitude;
		}
		else if( line.find("SEA_DEPTH") != std::string::npos ){
			ifs >> m_fileNameOfSeaDepthData;
			ifs >> m_minSeaDepth;
			ifs >> m_maxSeaDepth;
		}
		else if( line.find("NUM_THREADS") != std::string::npos ){
			ifs >> m_numThreads;
			if( m_numThreads <= 0 ){
				OutputFiles::m_logFile << " Error : Total number of threads is less than 1 : " << m_numThreads << std::endl;
				exit(1);
			}
		}
		else if( line.find("EXTENDED_REGION") != std::string::npos ){
			ifs >> m_lengthOfExtendedRegionPlusX;
			ifs >> m_lengthOfExtendedRegionMinusX;
			ifs >> m_lengthOfExtendedRegionPlusY;
			ifs >> m_lengthOfExtendedRegionMinusY;
			m_addExtendedRegion = true;
		}
		else if( line.find("HEIGHT_DEVIATION") != std::string::npos ){
			ifs >> m_thresholdDeviationOfHegith;
			if( m_thresholdDeviationOfHegith <= 0.0 ){
				OutputFiles::m_logFile << " Error : Threshold value of deviation of height must be positive !! : " << m_thresholdDeviationOfHegith << std::endl;
				exit(1);
			}
		}
		else if( line.find("IDW_FACTOR") != std::string::npos ){
			ifs >> dbuf;
			if( dbuf < 0.0 ){
				OutputFiles::m_logFile << " Error : Factor of inverse distance weighting must be positive !! : " << dbuf << std::endl;
				exit(1);
			}
			m_factorInverseDistanceWeighting = dbuf;
		}
		else if( line.find("SURF_MESH") != std::string::npos ){
			m_makeSurfMesOnAllBoundary = true;
		}
		else if( line.find("OUTER_EDGE_LENGTH") != std::string::npos ){
			ifs >> dbuf;
			( CoastLineList::getInstance() )->setMaxEdgeLengthAtOuterEdges(dbuf);
			setOuterEdgeLength = true;
		}
		else if( line.find("INVERT_YCOORD") != std::string::npos ){
			m_invertSignYcoord = true;
		}
		else if( line.find("EXCHANGE") != std::string::npos ){
			m_exhangeLandAndSeaAtStep1 = true;
		}
		else if( line.find("SEDIMENT") != std::string::npos ){
			m_includeSedimentLayer = true;
			ifs >> dbuf;
			if( dbuf < 0.0 ){
				OutputFiles::m_logFile << " Error : Thickness of sediment must be positive !! : " << dbuf << std::endl;
				exit(1);
			}
			m_thicknessSediment = dbuf;
		}
		else if( line.find("LAYERS") != std::string::npos ){
			m_includeLayers = true;
			int numLayers(-1);
			ifs >> numLayers;
			if( numLayers <= 0 ){
				OutputFiles::m_logFile << " Error : Number of layers must be greater than 0 !! : " << numLayers << std::endl;
				exit(1);
			}
			m_numLayers = numLayers;
			m_depthLayerInterfaces = new double[numLayers];
			for( int i = 0; i < numLayers; ++i ){
				ifs >> dbuf;
				if( dbuf <= 0.0 ){
					OutputFiles::m_logFile << " Error : Depth of layer interface must be positive !! : " << dbuf << std::endl;
					exit(1);
				}
				m_depthLayerInterfaces[i] = dbuf;
			}
		}
		else if( line.find("END") != std::string::npos ){
			break;
		}
		
	}

	if( m_addExtendedRegion && m_makeSurfMesOnAllBoundary ){
		OutputFiles::m_logFile << " Error : Extended region cannot be added if surface meshes are made on all boundary faces." << std::endl;
		exit(1);
	}

	m_minEdgeLength = m_ellipsoids.getMaxEdgeLength(0);
	m_maxEdgeLength = m_ellipsoids.getMaxEdgeLength( m_ellipsoids.getNumOfEllipsoids() - 1 );
	if( !setOuterEdgeLength ){
		( CoastLineList::getInstance() )->setMaxEdgeLengthAtOuterEdges(m_maxEdgeLength);
	}
	else if( ( CoastLineList::getInstance() )->getMaxEdgeLengthAtOuterEdges() > m_maxEdgeLength ){
		OutputFiles::m_logFile << " Error : Maximum edge length at outer-most edges must be shorter than the maximum edge length of ellipsoids !! : " << ( CoastLineList::getInstance() )->getMaxEdgeLengthAtOuterEdges() << std::endl;
		exit(1);
	}

	OutputFiles::m_logFile << "# Threshold ratio of edge length of boundary curves : " << m_thresholdEdgeLengthRatio << std::endl;
	OutputFiles::m_logFile << "# Maximum iteration number of the coast line smoothing : " << m_maxIterNumCoastLineSmoothing << std::endl;
	OutputFiles::m_logFile << "# Threshold angle of the coast line smoothing : " << m_thresholdAngleCoastLineSmoothing << std::endl;
	OutputFiles::m_logFile << "# Threshold distance of the coast line smoothing : " << m_thresholdDistanceCoastLineSmoothing << std::endl;
	OutputFiles::m_logFile << "# Maximum edge length at outer-most edges : " << ( CoastLineList::getInstance() )->getMaxEdgeLengthAtOuterEdges() << std::endl;
	OutputFiles::m_logFile << "# Factor of the size of the super triangle : " << m_factorSuperTriangleSize << std::endl;
	OutputFiles::m_logFile << "# Parameter alpha : " << m_alpha << std::endl;
	OutputFiles::m_logFile << "# Parameter beta : " << m_beta << std::endl;
	OutputFiles::m_logFile << "# Threshold angle for refinement [deg] : " << m_thresholdAngle << std::endl;
	OutputFiles::m_logFile << "# Maximum iteration number for refining triangles having long edges : " << m_maxIterNumRefiningEdgeLength << std::endl;
	OutputFiles::m_logFile << "# Maximum iteration number for refining triangles having large angle : " << m_maxIterNumRefiningLargeAngle << std::endl;
	OutputFiles::m_logFile << "# Maximum iteration number of laplacian method : " << m_iterNumLaplacianMethod << std::endl;
	OutputFiles::m_logFile << "# Maximum number of points used for interpolating altitudes : " << m_maxNumOfPointsForInterpolating << std::endl;
	OutputFiles::m_logFile << "# Distance used to avoid too small denominator in inverse distance weighting [km] : " << m_distanceUsedToAvoidTooSmallDenominator << std::endl;
	if( m_thresholdDeviationOfHegith > 0 ){
		OutputFiles::m_logFile << "# Threshold value of deviation of height : " << m_thresholdDeviationOfHegith << std::endl;
	}
	OutputFiles::m_logFile << "# Factor of inverse distance weighting : " << m_factorInverseDistanceWeighting << std::endl;
	OutputFiles::m_logFile << "# Maximum distance  for interpolating [km] : " << m_maxDistanceForInterpolating << std::endl;
	if( !m_fileNameOfAltitudeData.empty() ){
		OutputFiles::m_logFile << "# File name of altitude of the earth : " << m_fileNameOfAltitudeData.c_str() << std::endl;
		OutputFiles::m_logFile << "# Minimum value of altitude of land [km] : " << m_minAltitude << std::endl;
		OutputFiles::m_logFile << "# Maximum value of altitude of land [km] : " << m_maxAltitude << std::endl;
	}
	if( !m_fileNameOfSeaDepthData.empty() ){
		OutputFiles::m_logFile << "# File name of depth of the sea : " << m_fileNameOfSeaDepthData.c_str() << std::endl;
		OutputFiles::m_logFile << "# Minimum value of depth of the sea [km] : " << m_minSeaDepth << std::endl;
		OutputFiles::m_logFile << "# Maximum value of depth of the sea [km] : " << m_maxSeaDepth << std::endl;
	}
	if(m_makeSurfMesOnAllBoundary){
		OutputFiles::m_logFile << "# Surface meshes are made on all boundary faces." << std::endl;
	}
	else{
		OutputFiles::m_logFile << "# Surface meshes are made only on the surface of the earth." << std::endl;
	}
	OutputFiles::m_logFile << "# Total number of threads : " << m_numThreads << std::endl;
	if(m_addExtendedRegion){ 
		OutputFiles::m_logFile << "# Add extended region" << std::endl;
		OutputFiles::m_logFile << "#  Length of extended region in + X direction [km] : " << m_lengthOfExtendedRegionPlusX << std::endl;
		OutputFiles::m_logFile << "#  Length of extended region in - X direction [km] : " << m_lengthOfExtendedRegionMinusX << std::endl;
		OutputFiles::m_logFile << "#  Length of extended region in + Y direction [km] : " << m_lengthOfExtendedRegionPlusX << std::endl;
		OutputFiles::m_logFile << "#  Length of extended region in - Y direction [km] : " << m_lengthOfExtendedRegionMinusX << std::endl;
	}
	if(m_invertSignYcoord){
		OutputFiles::m_logFile << "# Invert sign of Y coordinates." << std::endl;
	}
	if(m_exhangeLandAndSeaAtStep1){
		OutputFiles::m_logFile << "# The exchange of the types of boundary curves is performed at the step1." << std::endl;
	}

	if(m_includeSedimentLayer){
		OutputFiles::m_logFile << "# Thickness of sediment layer [km] : " << m_thicknessSediment << std::endl;
	}

	if(m_includeLayers){
		OutputFiles::m_logFile << "# Number of layers : " << m_numLayers << std::endl;
		OutputFiles::m_logFile << "# Depth of layer interfaces : " << std::endl;
		for( int i = 0; i < m_numLayers; ++i ){
			OutputFiles::m_logFile << "# " << m_depthLayerInterfaces[i] << " [km]" << std::endl;
		}
	}

	ifs.close();

}

// Run mesh generation
void Control::run( int argc, char *argv[] ){

	if( argc < 3 ){
		std::cerr << "You must specify  arguments as below" << std::endl;
		std::cerr << " " << CommonParameters::programName << " -stp [Step ID]" << std::endl;
		std::cerr << " <Step>" << std::endl;
		std::cerr << "  1 : Calculate boudary curve of the earth surface" << std::endl;
		exit(1);
	}

	OutputFiles::getInstance();

	readControlData();

	const int stepID = atoi(argv[2]);
	switch( stepID ){
		case 1:
			OutputFiles::m_logFile << "# *** Start Step " << stepID << std::endl;
			OutputFiles::m_logFile << "# *** Calculate boudary curve of the earth surface" << std::endl;
			calcBoundaryCurveEarthSurface();
			break;
		case 2:
			OutputFiles::m_logFile << "# *** Start Step " << stepID << std::endl;
			OutputFiles::m_logFile << "# *** Calculate surface mesh of the earth" << std::endl;
			calcSurfaceMesh();
			break;
		case 3:
			OutputFiles::m_logFile << "# *** Start Step " << stepID << std::endl;
			OutputFiles::m_logFile << "# *** Calculate surface mesh with topography" << std::endl;
			calcSurfaceMeshWithTopograpy();
			break;
		case 4:
			OutputFiles::m_logFile << "# *** Start Step " << stepID << std::endl;
			OutputFiles::m_logFile << "# *** Make PLCs with surface meshes of the eatrh" << std::endl;
			//makePLCsWithOnlySurfaceMeshesOfEarth();
			makePLCs();
			break;
		default:
			OutputFiles::m_logFile << " Step ID is wrong. ID = " << stepID << std::endl;
			break;
	}

}

// Calculate maximum length of specified coordinate  on the specified plane from control parameter only
double Control::calcMaximumLengthFromControlParamOnly( const CommonParameters::XY& coord, const CommonParameters::Boundary& planeType, const int iLayer ) const{
	if( planeType == CommonParameters::LAYER ){
		const double coordZ = getDepthLayerInterfaces(iLayer);
		return calcMaximumLengthFromControlParamOnly( Util::calcCoordOf3DModel(coord, planeType, coordZ) );
	}
	else{
		return calcMaximumLengthFromControlParamOnly( Util::calcCoordOf3DModel(coord, planeType) );
	}
}

// Calculate maximum length of specified coordinate from control parameter only
double Control::calcMaximumLengthFromControlParamOnly( const CommonParameters::XYZ& coord ) const{
	return m_ellipsoids.calcMaximumEdgeLength( coord );
}

//// Calculate maximum length of specified coordinate from control parameter only 
//double Control::calcMaximumLengthFromControlParamOnly( const CommonParameters::XYZ& coord ) const{
//
//	const double vecXOrg = coord.X - m_centerCoordFineZone.X;
//	const double vecYOrg = coord.Y - m_centerCoordFineZone.Y;
//
//	// Coordinate transform
//	const double vecX = vecXOrg * cos( - m_rotationAngleOfFineZone ) - vecYOrg * sin( - m_rotationAngleOfFineZone );
//	const double vecY = vecXOrg * sin( - m_rotationAngleOfFineZone ) + vecYOrg * cos( - m_rotationAngleOfFineZone );
//	const double vecZ = coord.Z - m_centerCoordFineZone.Z;
//
//	const double val1 = pow( vecX/m_axisLengthXOfFineZone, 2.0 ) + pow( vecY/m_axisLengthYOfFineZone, 2.0 ) + pow( vecZ/m_axisLengthZOfFineZone, 2.0 );
//	const double val2 = val1 * pow( 1.0/m_axisRatioOfEllipsoids, 2.0 );
//
//	if( val1 <= 1.0 ){
//		return m_maxEdgeLengthOfFineZone;
//	}
//	else if( val2 >= 1.0 ){
//		return m_maxEdgeLengthOfRoughZone;
//	}
//	
//	return m_maxEdgeLengthOfFineZone + ( m_maxEdgeLengthOfRoughZone - m_maxEdgeLengthOfFineZone ) * ( sqrt(val1) - 1.0 ) / ( m_axisRatioOfEllipsoids - 1.0 );
//
//}

// Return flag whether specified point on the specified plane locate in the finest zone
bool Control::doesLocateInFinestZone( const CommonParameters::XY& coord, const CommonParameters::Boundary& m_planeType ) const{
	if( m_planeType != CommonParameters::SURFACE ){
		return false;
	}
	return m_ellipsoids.locateInEllipsoid( Util::calcCoordOf3DModel(coord, m_planeType), 0 );
}

//// Return flag whether specified point
//bool Control::doesLocateInFineZone( const CommonParameters::XYZ& coord ) const{
//
//	// Calculate vector on cartesian coordinates
//	const double vecXOrg = coord.X - m_centerCoordFineZone.X;
//	const double vecYOrg = coord.Y - m_centerCoordFineZone.Y;
//
//	// Coordinate transform
//	const double vecX = vecXOrg * cos( - m_rotationAngleOfFineZone ) - vecYOrg * sin( - m_rotationAngleOfFineZone );
//	const double vecY = vecXOrg * sin( - m_rotationAngleOfFineZone ) + vecYOrg * cos( - m_rotationAngleOfFineZone );
//	const double vecZ = coord.Z - m_centerCoordFineZone.Z;
//
//	const double val = pow( vecX/m_axisLengthXOfFineZone, 2.0 ) + pow( vecY/m_axisLengthYOfFineZone, 2.0 ) + pow( vecZ/m_axisLengthZOfFineZone, 2.0 );
//
//	if( val <= 1.0 ){
//		return true;
//	}
//
//	return false;
//
//}

//// Get maximum edge length of fine zone 
//double Control::getMaxEdgeLengthOfFineZone() const{
//	return m_maxEdgeLengthOfFineZone;
//}
//

// Get minimum edge length
double Control::getMinEdgeLength() const{
	return m_minEdgeLength;
}

// Get maximum edge length
double Control::getMaxEdgeLength() const{
	return m_maxEdgeLength;
}

// Get factor of super triangle size
double Control::getFactorSuperTriangleSize() const{
	return m_factorSuperTriangleSize;
}

// Calculate maximum edge length of this point
double Control::calcMaximumEdgeLength( const CommonParameters::XY& coord, const CommonParameters::Boundary& planeType, const int iLayer ) const{

	double length = planeType == CommonParameters::LAYER ? calcMaximumLengthFromControlParamOnly( coord, planeType, iLayer ) : calcMaximumLengthFromControlParamOnly( coord, planeType );

	if( planeType == CommonParameters::SURFACE ){
		double dbuf = ( ObservingSiteList::getInstance() )->calcMaximumLengthOfPoint( coord );
		if( dbuf < length ){
			length = dbuf;
		}
	}

	return length;

}

// Calculate maximum edge length of this point
double Control::calcMaximumEdgeLength( const CommonParameters::XYZ& coord ) const{

	const double length = calcMaximumLengthFromControlParamOnly( coord );

	return length;

}

// Get maximum iteration number for refining triangles having large angle
int Control::getMaxIterNumRefiningTriangleWithLargeAngle() const{
	return m_maxIterNumRefiningLargeAngle;
}

// Get maximum iteration number for refining triangles having long edges
int Control::getMaxIterNumRefiningTriangleWithLongEdge() const{
	return m_maxIterNumRefiningEdgeLength;
}
	
// Get threshold angle for refining triangles having large angle
double Control::getThresholdAngle() const{
	return m_thresholdAngle;
}

// Get control parameter of refinment
double Control::getAlpha() const{
	return m_alpha;
}

// Get control parameter of refinment
double Control::getBeta() const{
	return m_beta;
}

// Get maximum iteration number of the coast line smoothing
int Control::getMaxIterNumCoastLineSmoothing() const{
	return m_maxIterNumCoastLineSmoothing;
}

// Get threshold angle of the coast line smoothing
double Control::getThresholdAngleCoastLineSmoothing() const{
	return m_thresholdAngleCoastLineSmoothing;
}

// Get threshold distance of the coast line smoothing
double Control::getThresholdDistanceCoastLineSmoothing() const{
	return m_thresholdDistanceCoastLineSmoothing;
}

// Get iteration number of laplacian method
int Control::getIterNumLaplacianMethod() const{
	return m_iterNumLaplacianMethod;
}

// Get threshold ratio for deleting points belonging two segments ratio of whose edge length is greater than the value
double Control::getThresholdEdgeLengthRatio() const{
	return m_thresholdEdgeLengthRatio;
}

// Get maximum number of points used for interpolating altitudes
int Control::getMaxNumOfPointsForInterpolating() const{
	return m_maxNumOfPointsForInterpolating;
}

// Get maximum distance used for interpolating altitudes
double Control::getMaxDistanceForInterpolating() const{
	return m_maxDistanceForInterpolating;
}

// Get distance used to avoid too small denominator in inverse distance weighting
double Control::getDistanceUsedToAvoidTooSmallDenominator() const{
	return m_distanceUsedToAvoidTooSmallDenominator;
}

// Get threshold value of deviation of height
double Control::getThresholdDeviationOfHegith() const{
	return m_thresholdDeviationOfHegith;
}

// Get minimum value of altitude of the earth
double Control::getMinAltitude() const{
	return m_minAltitude;
}

// Get maximum value of altitude of the earth
double Control::getMaxAltitude() const{
	return m_maxAltitude;
}

// Get minimum value of depth of the sea
double Control::getMinSeaDepth() const{
	return m_minSeaDepth;
}

// Get maximum value of depth of the sea
double Control::getMaxSeaDepth() const{
	return m_maxSeaDepth;
}

// Get total number of threads
int Control::getNumThreads() const{
	return m_numThreads;
}

// Flag specifing whether add extended region
bool Control::isExtendedRegionAdded() const{
	return m_addExtendedRegion;
}

// Length of extended region in + X direction
double Control::getLengthOfExtendedRegionPlusX() const{
	return m_lengthOfExtendedRegionPlusX;
}

// Length of extended region in - X direction
double Control::getLengthOfExtendedRegionMinusX() const{
	return m_lengthOfExtendedRegionMinusX;
}

// Length of extended region in + Y direction
double Control::getLengthOfExtendedRegionPlusY() const{
	return m_lengthOfExtendedRegionPlusY;
}

// Length of extended region in - Y direction
double Control::getLengthOfExtendedRegionMinusY() const{
	return m_lengthOfExtendedRegionMinusY;
}

// Get flag specifing whether surface meshes are made on all boundary faces
bool Control::getFlagWhetherMakeSurfMesOnAllBoundary() const{
	return m_makeSurfMesOnAllBoundary;
}

// Get factor of inverse distance weighting
double Control::getFactorInverseDistanceWeighting() const{
	return m_factorInverseDistanceWeighting;
}

// Get flag specifing whether invert sign of Y coordinate
bool Control::getInvertSignYcoord() const{
	return m_invertSignYcoord;
}

// Get flag specifing whether sediment layer is included
bool Control::getIncludeSedimentLayer() const{
	return m_includeSedimentLayer;
}

// Get thickness of sediment layer
double Control::getThicknessSediment() const{
	return m_thicknessSediment;
}

// Get flag specifing whether some layers are included
bool Control::getIncludeLayers() const{
	return m_includeLayers;
}

// Get number of layers
int Control::getNumLayers() const{
	return m_numLayers;
}

// Get depth of layer interfaces
double Control::getDepthLayerInterfaces( const int iLayer ) const{
	if( iLayer < 0 ||  iLayer >= m_numLayers ){
		OutputFiles::m_logFile << " Error : Index of layer is improper !! : " << iLayer << std::endl;
		exit(1);
	}
	return m_depthLayerInterfaces[iLayer];
}

// Calculate boudary curve of the earth surface
void Control::calcBoundaryCurveEarthSurface() const{

	( Control::getInstance() )->readControlData();
	( AnalysisDomain::getInstance() )->readAnalysisDomainData();
	( ObservingSiteList::getInstance() )->readObservingSiteData();

	// Read coast line data
	CoastLineList* ptrCoastLineList = CoastLineList::getInstance();
	ptrCoastLineList->readCoastLineData();
	ptrCoastLineList->writeCoastLineDataToVTK( "coastLine_fine.vtk" );

	// Roughen coast lines
	ptrCoastLineList->roughenCoastLine();
	ptrCoastLineList->writeCoastLineDataToVTK( "coastLine_rough.vtk" );

	// Make boundary curves
	BoundaryCurveList* ptrBoundaryCurveList = new BoundaryCurveList;
	ptrBoundaryCurveList->makeBoundaryCurvesCoastLine();
#ifdef _MOD_FOR_NMT
	// Add NMT dipoles to boudary curve list
	ptrBoundaryCurveList->addNMTDipolesToBoundaryCurveList();
#endif

	ptrBoundaryCurveList->relateInnerBoundToOuterBound();
	ptrBoundaryCurveList->relateSubInnerBoundToInnerBound();
	//ptrBoundaryCurveList->addBoudaryCurveOfAnomalies();

	//// Relate node to boundary curve
	//ptrBoundaryCurveList->relateNodeToBoundaryCurve();

	// Make all node of boundary curve to be fixed
	//NodeList* ptrNode2DList = NodeList::getInstance();
	
	ptrBoundaryCurveList->makeAllNodeFixed();

	if( m_exhangeLandAndSeaAtStep1 ){
		// Exchange the types of boundary curves (Sea => Land, Land => Sea)
		ptrBoundaryCurveList->exchangeTypesOfBoundaryCurves();
	}

	// Output files for paraview
	ptrBoundaryCurveList->writeBoundaryCurveListVTK( "boundary_curve_stp1.vtk" );
	ptrBoundaryCurveList->writeNode2DListToVTK( "boundary_curve_node_stp1.vtk" );

	// Output data of boundary curves to intermediate files
	ptrBoundaryCurveList->writeBoundaryCurveList( "boundary_curves.stp1" );
	ptrBoundaryCurveList->writeBoundaryCurveRelatios( "boundary_curves_relation.stp1" );
	ptrBoundaryCurveList->writeNode2DListToIntermediateFileStep1( "node_bound_curve.stp1" );

	OutputFiles::m_logFile << "# Finished !!" << std::endl;

	delete ptrBoundaryCurveList;
}

// Calculate surface mesh
void Control::calcSurfaceMesh() const{

	( AnalysisDomain::getInstance() )->readAnalysisDomainData();
	( ObservingSiteList::getInstance() )->readObservingSiteData();

	TriangleList m_triangleList;
	BoundaryCurveList* ptrBoundaryCurveList = m_triangleList.getPointerToBoundaryCurveList();

	// Read data from intermediate files
	ptrBoundaryCurveList->readBoundaryCurveList( "boundary_curves.stp1" );
	ptrBoundaryCurveList->readBoundaryCurveRelatios( "boundary_curves_relation.stp1" );
	//ptrNode2DList->readNode2DListToIntermediateFileStep1( "node2d_list.stp1" );
	ptrBoundaryCurveList->readNode2DListFromIntermediateFile( "node_bound_curve.stp1" );

#ifdef _DEBUG_WRITE
	// Debug write intermediate files
	ptrBoundaryCurveList->writeBoundaryCurveList( "boundary_curves.dbg" );
	ptrBoundaryCurveList->writeBoundaryCurveRelatios( "boundary_curves_relation.dbg" );
	//ptrNode2DList->writeNode2DListToIntermediateFileStep1( "node2d_list.dbg" );
	ptrBoundaryCurveList->writeNode2DListToIntermediateFileStep1( "node2d_list.dbg" );
#endif

	//// Relate node to boundary curve
	//ptrBoundaryCurveList->relateNodeToBoundaryCurve();

	// Create surface triangles
	TriangleList* ptrTriangleList = &m_triangleList;
	ptrTriangleList->createSurfaceTriangles();

	// Output data of triangle to intermediate files
	ptrTriangleList->writeTrianglesToIntermediateFileStep2( "triangle_list.stp2" );
	ptrTriangleList->writeNode2DListToIntermediateFileStep2( "node_mesh.stp2" );

	// Output files for paraview
	ptrBoundaryCurveList->writeBoundaryCurveListVTK( "boundary_curve_stp2.vtk" );
	ptrBoundaryCurveList->writeNode2DListToVTK( "boundary_curve_node_stp2.vtk" );

	// Output data of boundary curves to intermediate files
	ptrBoundaryCurveList->writeBoundaryCurveList( "boundary_curves.stp2" );
	ptrBoundaryCurveList->writeBoundaryCurveRelatios( "boundary_curves_relation.stp2" );
	ptrBoundaryCurveList->writeNode2DListToIntermediateFileStep1( "node_bound_curve.stp2" );

	// Write array convert node ID of boundary curve list to the one of triangle list 
	ptrTriangleList->writeNodeIDBoundCurve2Triangle( "node_curves2mesh.stp2" );

	OutputFiles::m_logFile << "# Finished !!" << std::endl;

}

// Calculate surface mesh with topograpy
void Control::calcSurfaceMeshWithTopograpy() const{

	// Read data of triangles from intermediate files
	TriangleList m_triangleList;
	TriangleList* ptrTriangleList = &m_triangleList;

	ptrTriangleList->readNode2DListFromIntermediateFile( "node_mesh.stp2" );
	ptrTriangleList->readTrianglesFromIntermediateFile( "triangle_list.stp2" );
	ptrTriangleList->relateNodesToTriangles();

	BoundaryCurveList* ptrBoundaryCurveList = m_triangleList.getPointerToBoundaryCurveList();
	ptrBoundaryCurveList->readBoundaryCurveList( "boundary_curves.stp2" );
	ptrBoundaryCurveList->readBoundaryCurveRelatios( "boundary_curves_relation.stp2" );
	ptrBoundaryCurveList->readNode2DListFromIntermediateFile( "node_bound_curve.stp2" );

	// Read array convert node ID of boundary curve list to the one of triangle list 
	ptrTriangleList->readNodeIDBoundCurve2Triangle( "node_curves2mesh.stp2" );

	( AnalysisDomain::getInstance() )->readAnalysisDomainData();

	// Read topography data
#ifdef _TOPO_FUNC
#else
	TopographyDataList* ptrTopographyDataList = TopographyDataList::getInstance();
	if( !m_fileNameOfAltitudeData.empty() ){
		ptrTopographyDataList->readAltitudeData( m_fileNameOfAltitudeData );
	}
	if( !m_fileNameOfSeaDepthData.empty() ){
		ptrTopographyDataList->readSeaDepthData( m_fileNameOfSeaDepthData );
	}
#endif

	// Assign location flag to nodes
	ptrTriangleList->assignLocationToNodes();

	// Read lake data
	( LakeList::getInstance() )->readLakeData();
	ptrBoundaryCurveList->setupRelationBetweenNodesOnBoundaryAndLakeData();

	// Interpolate altitude and sea depth to nodes
	ptrTriangleList->interpolateAltitudeToNodes();

	ptrTriangleList->writeTrianglesToIntermediateFileStep3( "triangle_list.stp3" );
	ptrTriangleList->writeNodeListToIntermediateFileStep3( "node_mesh.stp3" );

	// Output files for paraview
	ptrBoundaryCurveList->writeBoundaryCurveListVTK( "boundary_curve_stp3.vtk" );
	ptrBoundaryCurveList->writeNode2DListToVTK( "boundary_curve_node_stp3.vtk" );

	// Output data of boundary curves to intermediate files
	ptrBoundaryCurveList->writeBoundaryCurveList( "boundary_curves.stp3" );
	ptrBoundaryCurveList->writeBoundaryCurveRelatios( "boundary_curves_relation.stp3" );
	ptrBoundaryCurveList->writeNode2DListToIntermediateFileStep1( "node_bound_curve.stp3" );

	// Write array convert node ID of boundary curve list to the one of triangle list 
	ptrTriangleList->writeNodeIDBoundCurve2TriangleStep3( "node_curves2mesh.stp3" );

	OutputFiles::m_logFile << "# Finished !!" << std::endl;

}

//// Make PLCs with surface meshes of the eatrh
//void Control::makePLCsWithOnlySurfaceMeshesOfEarth() const{
// Make PLCs
void Control::makePLCs() const{

	// Read data of triangles from intermediate files
	TriangleList m_triangleList;
	TriangleList* ptrTriangleList = &m_triangleList;

	ptrTriangleList->readNodeListFromIntermediateFileStep3( "node_mesh.stp3" );
	ptrTriangleList->readTrianglesFromIntermediateFile( "triangle_list.stp3" );
	ptrTriangleList->relateNodesToTriangles();

	BoundaryCurveList* ptrBoundaryCurveList = m_triangleList.getPointerToBoundaryCurveList();
	ptrBoundaryCurveList->readBoundaryCurveList( "boundary_curves.stp3" );
	ptrBoundaryCurveList->readBoundaryCurveRelatios( "boundary_curves_relation.stp3" );
	ptrBoundaryCurveList->readNode2DListFromIntermediateFile( "node_bound_curve.stp3" );

	// Read array convert node ID of boundary curve list to the one of triangle list 
	ptrTriangleList->readNodeIDBoundCurve2Triangle( "node_curves2mesh.stp3" );

	( AnalysisDomain::getInstance() )->readAnalysisDomainData();

	// Read lake data
	( LakeList::getInstance() )->readLakeData();
	ptrBoundaryCurveList->setupRelationBetweenNodesOnBoundaryAndLakeData();

	// Make and write PLCs
	ptrTriangleList->writePLCs();

	OutputFiles::m_logFile << "# Finished !!" << std::endl;

}


