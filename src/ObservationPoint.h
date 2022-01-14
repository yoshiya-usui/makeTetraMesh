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
#ifndef DBLDEF_OBSERVATION_POINT
#define DBLDEF_OBSERVATION_POINT

#include <fstream>
#include <iostream>
#include "CommonParameters.h"

// Class of observation point
class ObservationPoint{

public:
	// Default Constructer
	ObservationPoint();

	// Destructer
	~ObservationPoint();

	// Read data of observation point from input file
	void readObservationPointData( std::ifstream& ifs );

	//// Calculate maximum length of specified coordinate
	//double calcMaximumLengthOfPoint( const CommonParameters::XYZ& coord ) const;

	// Calculate maximum length of specified coordinate on X-Y plane ( z = 0 ) 
	double calcMaximumLengthOfPoint( const CommonParameters::XY& coord ) const;

private:
	// Copy constructer
	ObservationPoint(const ObservationPoint& rhs);

	// Assignment operator
	ObservationPoint& operator=(ObservationPoint& rhs);

	// Coordinate of the point
	//CommonParameters::XYZ m_pointCoord;
	CommonParameters::XY m_pointCoord;

	// Total number of circles
	int m_numCircles;

	// Radius for local mesh control
	double* m_radius;

	// Maximum edge length within circle around the point
	double* m_maxEdgeLengthWithinCircle;

};

#endif
