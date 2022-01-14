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
#ifndef DBLDEF_CONTROL
#define DBLDEF_CONTROL

#include "CommonParameters.h"
#include "Ellipsoids.h"

// Class of the control parameters
class Control{

public:
	// Return the the instance of the class
    static Control* getInstance();

	// Read control parameters from input file
	void readControlData();

	// Run mesh generation
	void run( int argc, char *argv[] );
	
	// Calculate maximum length of specified coordinate on the specified plane from control parameter only
	double calcMaximumLengthFromControlParamOnly( const CommonParameters::XY& coord, const CommonParameters::Boundary& planeType, const int iLayer = -1 ) const;

	// Calculate maximum length of specified coordinate from control parameter only
	double calcMaximumLengthFromControlParamOnly( const CommonParameters::XYZ& coord ) const;

	// Return flag whether specified point on X-Y plane ( z = 0 ) locate in the finest zone
	bool doesLocateInFinestZone( const CommonParameters::XY& coord, const CommonParameters::Boundary& m_planeType ) const;

	//// Return flag whether specified point
	//bool doesLocateInFineZone( const CommonParameters::XYZ& coord ) const;

	//// Get maximum edge length of fine zone 
	//double getMaxEdgeLengthOfFineZone() const;

	//// Get maximum edge length of rough zone 
	//double getMaxEdgeLengthOfRoughZone() const;

	// Get minimum edge length
	double getMinEdgeLength() const;

	// Get maximum edge length
	double getMaxEdgeLength() const;

	// Get factor of super triangle size
	double getFactorSuperTriangleSize() const;

	// Calculate maximum edge length of this points 
	double calcMaximumEdgeLength( const CommonParameters::XY& coord, const CommonParameters::Boundary& planeType, const int iLayer = -1 ) const;

	// Calculate maximum edge length of this points 
	double calcMaximumEdgeLength( const CommonParameters::XYZ& coord ) const;

	// Get maximum iteration number for refining triangles having long edges
	int getMaxIterNumRefiningTriangleWithLargeAngle() const;

	// Get maximum iteration number for refining triangles having large angle
	int getMaxIterNumRefiningTriangleWithLongEdge() const;
	
	// Get threshold angle for refining triangles having large angle
	double getThresholdAngle() const;

	// Get control parameter of refinment
	double getAlpha() const;
	double getBeta() const;

	// Get maximum iteration number of the coast line smoothing
	int getMaxIterNumCoastLineSmoothing() const;

	// Get threshold angle of the coast line smoothing
	double getThresholdAngleCoastLineSmoothing() const;
	
	// Get threshold distance of the coast line smoothing
	double getThresholdDistanceCoastLineSmoothing() const;

	// Get iteration number of laplacian method
	int getIterNumLaplacianMethod() const;

	// Get threshold ratio for deleting points belonging two segments ratio of whose edge length is greater than the value
	double getThresholdEdgeLengthRatio() const;

	// Get maximum number of points used for interpolating altitudes
	int getMaxNumOfPointsForInterpolating() const;

	// Get maximum distance used for interpolating altitudes
	double getMaxDistanceForInterpolating() const;

	// Get distance used to avoid too small denominator in inverse distance weighting
	double getDistanceUsedToAvoidTooSmallDenominator() const;

	// Get threshold value of deviation of height
	double getThresholdDeviationOfHegith() const;

	// Get minimum value of altitude of the earth
	double getMinAltitude() const;

	// Get maximum value of altitude of the earth
	double getMaxAltitude() const;

	// Get minimum value of depth of the sea
	double getMinSeaDepth() const;

	// Get maximum value of depth of the sea
	double getMaxSeaDepth() const;

	// Get total number of threads
	int getNumThreads() const;

	// Flag specifing whether add extended region
	bool isExtendedRegionAdded() const;

	// Length of extended region in + X direction
	double getLengthOfExtendedRegionPlusX() const;

	// Length of extended region in - X direction
	double getLengthOfExtendedRegionMinusX() const;

	// Length of extended region in + Y direction
	double getLengthOfExtendedRegionPlusY() const;

	// Length of extended region in - Y direction
	double getLengthOfExtendedRegionMinusY() const;

	// Get flag specifing whether surface meshes are made on all boundary faces
	bool getFlagWhetherMakeSurfMesOnAllBoundary() const;

	// Get factor of inverse distance weighting
	double getFactorInverseDistanceWeighting() const;

	// Get flag specifing whether invert sign of Y coordinate
	bool getInvertSignYcoord() const;

	// Get flag specifing whether sediment layer is included
	bool getIncludeSedimentLayer() const;

	// Get thickness of sediment layer
	double getThicknessSediment() const;

	// Get flag specifing whether some layers are included
	bool getIncludeLayers() const;

	// Get number of layers
	int getNumLayers() const;

	// Get depth of layer interfaces
	double getDepthLayerInterfaces( const int iLayer ) const;

private:
	// Constructer
	Control();

	// Destructer
	~Control();

	// Copy constructer
	Control(const Control& rhs);

	// Assignment operator
	Control& operator=(const Control& rhs);

	const static double m_eps;

	// Ellipsoids used for specifing edge lentgh
	Ellipsoids m_ellipsoids;

	// Minimum edge length
	double m_minEdgeLength;

	// Maximum edge length
	double m_maxEdgeLength;

	// Threshold ratio for deleting points belonging two segments ratio of whose edge length is greater than the value
	double m_thresholdEdgeLengthRatio;

	// Factor of super triangle size
	double m_factorSuperTriangleSize;

	// Maximum iteration number for refining triangles having large angle
	int m_maxIterNumRefiningLargeAngle;

	// Maximum iteration number for refining triangles having long edges
	int m_maxIterNumRefiningEdgeLength;

	// Threshold angle for refining triangles having large angle
	double m_thresholdAngle;

	// Control parameter of refinment
	double m_alpha;
	double m_beta;

	// Maximum iteration number of the coast line smoothing
	int m_maxIterNumCoastLineSmoothing;

	// Threshold angle of the coast line smoothing
	double m_thresholdAngleCoastLineSmoothing;

	// Threshold distance of the coast line smoothing
	double m_thresholdDistanceCoastLineSmoothing;

	// Iteration number of laplacian method
	int m_iterNumLaplacianMethod;

	// File name of altitude data
	std::string m_fileNameOfAltitudeData;

	// File name of sea depth data
	std::string m_fileNameOfSeaDepthData;

	// Maximum number of points used for interpolating altitudes
	int m_maxNumOfPointsForInterpolating;

	// Maximum distance used for interpolating altitudes
	double m_maxDistanceForInterpolating;

	// Distance used to avoid too small denominator in inverse distance weighting
	double m_distanceUsedToAvoidTooSmallDenominator;

	// Threshold value of deviation of height
	double m_thresholdDeviationOfHegith;

	// Minimum value of altitude of the earth
	double m_minAltitude;

	// Maximum value of altitude of the earth
	double m_maxAltitude;

	// Minimum value of depth of the sea
	double m_minSeaDepth;

	// Maximum value of depth of the sea
	double m_maxSeaDepth;

	// Total number of threads
	int m_numThreads;

	// Flag specifing whether add extended region
	bool m_addExtendedRegion;

	// Length of extended region in + X direction
	double m_lengthOfExtendedRegionPlusX;

	// Length of extended region in - X direction
	double m_lengthOfExtendedRegionMinusX;

	// Length of extended region in + Y direction
	double m_lengthOfExtendedRegionPlusY;

	// Length of extended region in - Y direction
	double m_lengthOfExtendedRegionMinusY;

	// Flag specifing whether surface meshes are made on all boundary faces
	bool m_makeSurfMesOnAllBoundary;

	// Flag specifing whether invert sign of Y coordinate
	bool m_invertSignYcoord;

	// Factor of inverse distance weighting
	double m_factorInverseDistanceWeighting;

	// Flag specifing whether the exhange of land with sea is performed at the step 1
	bool m_exhangeLandAndSeaAtStep1;

	// Flag specifing whether sediment layer is included
	bool m_includeSedimentLayer;

	// Thickness of sediment layer
	double m_thicknessSediment;

	// Flag specifing whether some layers are included
	bool m_includeLayers;

	// Number of layers
	int m_numLayers;

	// Depth of layer interfaces
	double* m_depthLayerInterfaces;

	// Calculate boudary curve of the earth surface
	void calcBoundaryCurveEarthSurface() const;

	// Calculate surface mesh
	void calcSurfaceMesh() const;

	// Calculate surface mesh with topograpy
	void calcSurfaceMeshWithTopograpy() const;

	// Make PLCs
	void makePLCs() const;

};

#endif
