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
