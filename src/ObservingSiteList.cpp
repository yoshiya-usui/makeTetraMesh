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
#include "ObservingSiteList.h"
#include "Control.h"
#include "OutputFiles.h"
#include <stdlib.h>
#include <fstream>
#include <iostream>
#ifdef _MOD_FOR_NMT
#include <assert.h> 
#endif

// Return the instance of the class
ObservingSiteList* ObservingSiteList::getInstance(){
   	static ObservingSiteList instance;// The only instance
  	return &instance;
}

// Default constructer
ObservingSiteList::ObservingSiteList():
	m_numObservationPoint(0),
	m_numObservationLine(0),
	m_obsPoint(NULL),
	m_obsLine(NULL)
{

}

// Destructer
ObservingSiteList::~ObservingSiteList(){

	if( m_obsPoint != NULL ){
		delete [] m_obsPoint;
		m_obsPoint = NULL;
	}

	if( m_obsLine != NULL ){
		delete [] m_obsLine;
		m_obsLine = NULL;
	}

}

// Read data of observing site list from input file
void ObservingSiteList::readObservingSiteData(){

	std::ifstream ifs( "observing_site.dat", std::ios::in );

	if( ifs.fail() ){
		OutputFiles::m_logFile << "Error : File open error : observing_site.dat !!" << std::endl;
		exit(1);
	}

	OutputFiles::m_logFile << "# Read data of observing site list from input file." << std::endl;

	ifs >> m_numObservationPoint;

#ifdef _DEBUG_WRITE
	std::cout << "m_numObservationPoint = " << m_numObservationPoint << std::endl;
#endif

	if( m_numObservationPoint > 0 ){
		m_obsPoint = new ObservationPoint[m_numObservationPoint];
	}
	//else{
	//	OutputFiles::m_logFile << "Error : Total number of observation point is less than or equal to zero !!" << std::endl;
	//	exit(1);
	//}

	for( int i = 0; i < m_numObservationPoint; ++i ){
		m_obsPoint[i].readObservationPointData( ifs );
    }

	ifs >> m_numObservationLine;

#ifdef _DEBUG_WRITE
	std::cout << "m_numObservationLine = " << m_numObservationLine << std::endl;
#endif

	if( m_numObservationLine > 0 ){
		m_obsLine = new ObservationLine[m_numObservationLine];
	}
	//else{
	//	OutputFiles::m_logFile << "Error : Total number of observation line is less than or equal to zero !!" << std::endl;
	//	exit(1);
	//}

	for( int i = 0; i < m_numObservationLine; ++i ){
		//(m_obsLine[i].getPtrEndPoints(0))->readObservationPointData( ifs );
		//(m_obsLine[i].getPtrEndPoints(1))->readObservationPointData( ifs );
		m_obsLine[i].readObservationLineData( ifs );
    }

	ifs.close();

}

//// Calculate maximum length of specified coordinate
//double ObservingSiteList::calcMaximumLengthOfPoint( const CommonParameters::XYZ& coord ) const{
//
//	double val = ( Control::getInstance() )->calcMaximumLengthFromControlParamOnly( coord );
//
//	for( int i = 0; i < m_numObservationPoint; ++i ){
//		const double dbuf = m_obsPoint[i].calcMaximumLengthOfPoint( coord );
//		if( dbuf < val ){
//			val = dbuf;
//		}
//    }
//
//	for( int i = 0; i < m_numObservationLine; ++i ){
//		const double dbuf = m_obsLine[i].calcMaximumLengthOfPoint( coord );
//		if( dbuf < val ){
//			val = dbuf;
//		}
//    }
//
//	return val;
//
//}

// Calculate maximum length of specified coordinate on X-Y plane ( z = 0 ) 
double ObservingSiteList::calcMaximumLengthOfPoint( const CommonParameters::XY& coord ) const{

	double val = 1.0e12;

	for( int i = 0; i < m_numObservationPoint; ++i ){
		const double dbuf = m_obsPoint[i].calcMaximumLengthOfPoint( coord );
		if( dbuf < val ){
			val = dbuf;
		}
    }

	for( int i = 0; i < m_numObservationLine; ++i ){
		const double dbuf = m_obsLine[i].calcMaximumLengthOfPoint( coord );
		if( dbuf < val ){
			val = dbuf;
		}
    }

	return val;

}
#ifdef _MOD_FOR_NMT
// Get number of observation line
int ObservingSiteList::getNumObsLine() const{
	return m_numObservationLine;
}

// Get coordinate of start point of an observation line
CommonParameters::XY ObservingSiteList::getCoordOfStartPointObsLine( const int iobs ) const{
	assert( iobs >= 0 && iobs < getNumObsLine() );
	return m_obsLine[iobs].getCoordOfStartPoint();
}

// Get coordinate of end point of an observation line
CommonParameters::XY ObservingSiteList::getCoordOfEndPointObsLine( const int iobs ) const{
	assert( iobs >= 0 && iobs < getNumObsLine() );
	return m_obsLine[iobs].getCoordOfEndPoint();
}

// Get maximum edge length of an observation line
double ObservingSiteList::getMaximumEdgeLengthObsLine( const int iobs ) const{
	assert( iobs >= 0 && iobs < getNumObsLine() );
	return m_obsLine[iobs].getMaximumEdgeLength();
}
#endif
