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
#ifndef DBLDEF_ELLIPSOIDS
#define DBLDEF_ELLIPSOIDS

#include "CommonParameters.h"
#include <fstream>
#include <iostream>
#include <string>

// Class of the ellipsoids deciding maximum edge length
class Ellipsoids{

public:
	// Constructer
	Ellipsoids();

	// Destructer
	~Ellipsoids();

	// Read parameters
	//void Ellipsoids::readParameters( std::ifstream& ifs  );
	void readParameters( const std::string& fileName  );

	// Calculate maximum edge length of the specified point
	double calcMaximumEdgeLength( const CommonParameters::XYZ& coord ) const;

	// Get total number of ellipsoids
	int getNumOfEllipsoids() const;

	// Get maximum edge length within ellipsoids
	double getMaxEdgeLength( const int iEllipsoid ) const;

	// Check whether specified point is located in the Ellipsoid
	bool locateInEllipsoid( const CommonParameters::XYZ& coord, const int iEllipsoid ) const;

private:
	enum TypeID{
		EARTH = 0,
		AIR = 1
	};

	// Copy constructer
	Ellipsoids(const Ellipsoids& rhs);

	// Assignment operator
	Ellipsoids& operator=(const Ellipsoids& rhs);

	// Center coord
	CommonParameters::XYZ m_centerCoord;

	// Rotation angle of elliptical fine zone in X-Y plane
	double m_rotationAngle;

	// Total number of ellipsoids
	int m_numEllipsoids;

	// Length of long axis
	double* m_radius;

	// Maximum edge length within ellipsoids
	double* m_edgeLength;

	// Oblateness of horizontal direction
	double* m_oblatenessHorizontal;

	// Oblateness of vertical direction
	double* m_oblateness[2];

	// Calculate maximum edge length between two spheres
	double calcEdgeLengthTransitionRegion( const CommonParameters::XYZ& coord, const int iEllipsoid ) const;

	// Calculate length on ellipsoid
	double calculateLengthOnEllipsoid( const double angleH, const double angleV, const int iEllipsoid, const int iType ) const;

};

#endif
