﻿//--------------------------------------------------------------------------
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
#include "CoastLineList.h"
#include "OutputFiles.h"
#include "AnalysisDomain.h"
#include "Control.h"
#include "math.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>

// Return the instance of the class
CoastLineList* CoastLineList::getInstance(){
   	static CoastLineList instance;// The only instance
  	return &instance;
}

// Constructer
CoastLineList::CoastLineList():
	m_numCoastLines(0),
	m_coastLines(NULL),
	m_maxEdgeLengthAtOuterEdges(-1.0)
{

}

// Destructer
CoastLineList::~CoastLineList(){

	if( m_coastLines != NULL ){
		delete [] m_coastLines;
		m_coastLines = NULL;
	}

}

// Read data of coast lines from input file
void CoastLineList::readCoastLineData(){

	std::ifstream ifs( "coast_line.dat", std::ios::in );

	if( ifs.fail() ){
		OutputFiles::m_logFile << "Error : File open error : coast_line.dat !!" << std::endl;
		exit(1);
	}

	OutputFiles::m_logFile << "# Read data of coast lines from input file." << std::endl;

	ifs >> m_numCoastLines;

	OutputFiles::m_logFile << "# Total number of coast line data : " << m_numCoastLines << std::endl;

#ifdef _DEBUG_WRITE
	std::cout << "m_numCoastLines = " << m_numCoastLines << std::endl;
#endif

	if( m_numCoastLines > 0 ){
		m_coastLines = new CoastLine[m_numCoastLines];
	}

	for( int i = 0; i < m_numCoastLines; ++i ){
		m_coastLines[i].readCoastLineData( ifs );
    }

	ifs.close();

}

// Get total number of coast lines
int CoastLineList::getNumCoastLines() const{
	return m_numCoastLines;
}

// Get pointer to the instances of coast lines
CoastLine* CoastLineList::getPointerToCoastLine( const int id ) const{

	if( id > m_numCoastLines || id < 0 ){
		OutputFiles::m_logFile << "Error : Coast line ID is out of range !! ID = " << id << std::endl;
		exit(1);
	}

	return &m_coastLines[id];

}

// Roughen coast lines
void CoastLineList::roughenCoastLine(){

	OutputFiles::m_logFile << "# Roughen coast lines" << std::endl;

	for( int i = 0; i < m_numCoastLines; ++i ){
		m_coastLines[i].roughenCoastLine();
    }
	
}

// Output coast line data to vtk file
void CoastLineList::writeCoastLineDataToVTK( const std::string& fileName ) const{

	// Open output vtk file -----
	std::ofstream ofsVTK( fileName.c_str() );
	if( !ofsVTK ) {
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << std::endl;
		exit(1);
	}

	ofsVTK << "# vtk DataFile Version 2.0" << std::endl;
	ofsVTK << "CoastLine" << std::endl;
	ofsVTK << "ASCII" << std::endl;
	ofsVTK << "DATASET UNSTRUCTURED_GRID" << std::endl;
	
	ofsVTK.precision(9);
	ofsVTK << std::fixed;

	// output data to vtk file -----
	std::vector<CommonParameters::XY> stack;

	for( int iLine = 0; iLine < m_numCoastLines; ++iLine ){
		if( !m_coastLines[iLine].isClosed() ){
			continue;
		}
		const int numPointOfSegment = m_coastLines[iLine].getNumOfPoints();
		if( numPointOfSegment > 0 ){
			for( int iPoint = 0; iPoint < numPointOfSegment; ++iPoint ){
				stack.push_back( m_coastLines[iLine].getCoordXYOfPoints(iPoint) );
			}
		}
	}

	const int numPointsAll = static_cast<int>( stack.size() );

	ofsVTK << "POINTS " << numPointsAll << " float" << std::endl;
	for( std::vector<CommonParameters::XY>::iterator itr = stack.begin(); itr != stack.end(); ++itr ){
		ofsVTK << itr->X << " " << itr->Y << " " << 0.0 << std::endl;
	}

	ofsVTK << "CELLS " << numPointsAll << " " << numPointsAll * 3 << std::endl;
	int icount = 0;
	for( int iLine = 0; iLine < m_numCoastLines; ++iLine ){
		if( !m_coastLines[iLine].isClosed() ){
			continue;
		}
		const int numPointOfSegment = m_coastLines[iLine].getNumOfPoints();
		if( numPointOfSegment > 0 ){
			const int IDFirst = icount; 
			for( int iPoint = 0; iPoint < numPointOfSegment - 1; ++iPoint ){
				const int ID = icount++;
				ofsVTK << 2 << " " << ID << " " << ID + 1 << std::endl;
			}
			ofsVTK << 2 << " " << icount++ << " " << IDFirst << std::endl;
		}
	}

	ofsVTK << "CELL_TYPES " << numPointsAll << std::endl;
	for( int i = 0; i < numPointsAll; ++i ){
		ofsVTK << "3" << std::endl;
	}

	ofsVTK << "CELL_DATA " << numPointsAll << std::endl;
	ofsVTK << "SCALARS PointID int" <<  std::endl;
	ofsVTK << "LOOKUP_TABLE default" <<  std::endl;
	for( int i = 0; i < numPointsAll; ++i ){
		ofsVTK << i << std::endl;
	}

	ofsVTK << "SCALARS Length float" <<  std::endl;
	ofsVTK << "LOOKUP_TABLE default" <<  std::endl;
	for( int iLine = 0; iLine < m_numCoastLines; ++iLine ){
		if( !m_coastLines[iLine].isClosed() ){
			continue;
		}
		const int numPointOfSegment = m_coastLines[iLine].getNumOfPoints();
		if( numPointOfSegment > 0 ){
			for( int iPoint = 0; iPoint < numPointOfSegment - 1; ++iPoint ){
				const double length = hypot( m_coastLines[iLine].getCoordXOfPoints(iPoint+1) - m_coastLines[iLine].getCoordXOfPoints(iPoint) , m_coastLines[iLine].getCoordYOfPoints(iPoint+1) - m_coastLines[iLine].getCoordYOfPoints(iPoint) );
				ofsVTK << length << std::endl;
			}
			const double length = hypot( m_coastLines[iLine].getCoordXOfPoints(0) - m_coastLines[iLine].getCoordXOfPoints(numPointOfSegment-1) , m_coastLines[iLine].getCoordYOfPoints(0) - m_coastLines[iLine].getCoordYOfPoints(numPointOfSegment-1) );
			ofsVTK << length << std::endl;
		}
	}

	ofsVTK.close();

}

// Output coast line data to txt file for GMT
void CoastLineList::writeCoastLineDataToText( const std::string& fineName ) const{

	std::ofstream ofs( fineName.c_str() );
	if( !ofs ) {
		OutputFiles::m_logFile << " Error : Cannot open file " << fineName.c_str() << std::endl;
		exit(1);
	}

	for( int iLine = 0; iLine < m_numCoastLines; ++iLine ){
		const int numPointOfSegment = m_coastLines[iLine].getNumOfPoints();
		if( numPointOfSegment > 0 ){
			for( int iPoint = 0; iPoint < numPointOfSegment; ++iPoint ){
				CommonParameters::XY coord = m_coastLines[iLine].getCoordXYOfPoints(iPoint);
				ofs << coord.X << " " << coord.Y << std::endl;
			}
			ofs << ">" << std::endl;
		}
	}

	ofs.close();
}

// Set maximum edge length at outer-most edges
void CoastLineList::setMaxEdgeLengthAtOuterEdges( const double length ){
	m_maxEdgeLengthAtOuterEdges = length;
}

// Get maximum edge length at outer-most edges
double CoastLineList::getMaxEdgeLengthAtOuterEdges() const{
	return m_maxEdgeLengthAtOuterEdges;
}

// Calculate maximum edge length of this points assming this is one of the point of the coast line
double CoastLineList::calcMaximumEdgeLengthForCoastLine( const CommonParameters::XY& coord ) const{

	const double length = ( Control::getInstance() )->calcMaximumEdgeLength( coord, CommonParameters::SURFACE );
	if( ( AnalysisDomain::getInstance() )->doesIntersectWithBoundary(coord) ){// Outer most edges
		return std::min( length, ( CoastLineList::getInstance() )->getMaxEdgeLengthAtOuterEdges() );
	}
	else{
		return length;
	}

}


//	//------ For debug >>>>>
//#ifdef _DEBUG_WRITE
//
//void CoastLineList::debugWriteIntersectionPoints(){
//
//	for( int i = 0; i < m_numCoastLines; ++i ){
//		m_coastLines[i].debugWriteIntersectionPoints();
//    }
//
//}
//
//#endif
//	//------ For debug <<<<<
