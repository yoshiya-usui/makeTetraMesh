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
#ifndef DBLDEF_OBSERVING_SITE_LIST
#define DBLDEF_OBSERVING_SITE_LIST

#include "ObservationPoint.h"
#include "ObservationLine.h"
#ifdef _MOD_FOR_NMT
#include "CommonParameters.h"
#endif

// Class of the list of the observing site
class ObservingSiteList{

public:
	// Return the the instance of the class
    static ObservingSiteList* getInstance();

	// Read data of observing site list from input file
	void readObservingSiteData();

	//// Calculate maximum length of specified coordinate
	//double calcMaximumLengthOfPoint( const CommonParameters::XYZ& coord ) const;

	// Calculate maximum length of specified coordinate on X-Y plane ( z = 0 ) 
	double calcMaximumLengthOfPoint( const CommonParameters::XY& coord ) const;

#ifdef _MOD_FOR_NMT
	// Get number of observation line
	int getNumObsLine() const;

	// Get coordinate of start point of an observation line
	CommonParameters::XY getCoordOfStartPointObsLine( const int iobs ) const;

	// Get coordinate of end point of an observation line
	CommonParameters::XY getCoordOfEndPointObsLine( const int iobs ) const;

	// Get maximum edge length of an observation line
	double getMaximumEdgeLengthObsLine( const int iobs ) const;
#endif

private:
	// Default constructer
	ObservingSiteList();

	// Destructer
	~ObservingSiteList();

	// Copy constructer
	ObservingSiteList(const ObservingSiteList& rhs);

	// Assignment operator
	ObservingSiteList& operator=(const ObservingSiteList& rhs);

	// Total number of observation point
	int m_numObservationPoint;

	// Total number of observation line
	int m_numObservationLine;

	// Array of observation point
	ObservationPoint* m_obsPoint;

	// Array of observation line
	ObservationLine* m_obsLine;

};

#endif
