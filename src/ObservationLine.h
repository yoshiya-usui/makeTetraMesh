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
#ifndef DBLDEF_OBSERVATION_LINE
#define DBLDEF_OBSERVATION_LINE

#include <fstream>
#include <iostream>
#include <map>
#include "CommonParameters.h"
#include "ObservationPoint.h"

// Class of observation line
class ObservationLine{

public:
	// Default Constructer
	ObservationLine();

	// Destructer
	~ObservationLine();

	// Calculate maximum length of specified coordinate on X-Y plane ( z = 0 ) 
	double calcMaximumLengthOfPoint( const CommonParameters::XY& coord ) const;

	// Read data of observation points consisting line from input file
	void readObservationLineData( std::ifstream& ifs );

#ifdef _MOD_FOR_NMT
	// Get coordinate of start point
	CommonParameters::XY getCoordOfStartPoint() const;

	// Get coordinate of end point
	CommonParameters::XY getCoordOfEndPoint() const;

	// Get maximum edge length
	double getMaximumEdgeLength() const;
#endif

private:
	// Copy constructer
	ObservationLine(const ObservationLine& rhs);

	// Assignment operator
	ObservationLine& operator=(ObservationLine& rhs);

	// Pair of the points
	CommonParameters::XY m_endPoints[2];

	// Total number of layers
	int m_numLayers;

	// Radius for local mesh control
	double* m_radius;

	// Maximum edge length within layer around the line
	double* m_maxEdgeLengthWithinLayer;

};

#endif
