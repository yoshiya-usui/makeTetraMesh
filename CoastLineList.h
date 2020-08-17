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
#ifndef DBLDEF_COAST_LINE_LIST
#define DBLDEF_COAST_LINE_LIST

#include "CoastLine.h"

// Class of the list of coast lines
class CoastLineList{

public:
	// Return the the instance of the class
    static CoastLineList* getInstance();

	// Read data of coast lines from input file
	void readCoastLineData();

	// Get total number of coast lines
	int getNumCoastLines() const;

	// Get pointer to the instances of coast lines
	CoastLine* getPointerToCoastLine( const int id ) const;
	
	// Roughen coast line
	void roughenCoastLine();

	// Output coast line data to vtk file
	void writeCoastLineDataToVTK( const std::string& fineName ) const;

	// Set maximum edge length at outer-most edges
	void setMaxEdgeLengthAtOuterEdges( const double length );

	// Get maximum edge length at outer-most edges
	double getMaxEdgeLengthAtOuterEdges() const;

	// Calculate maximum edge length of this points assming this is one of the point of the coast line
	double calcMaximumEdgeLengthForCoastLine( const CommonParameters::XY& coord ) const;

	//	//------ For debug >>>>>
//#ifdef _DEBUG_WRITE
//	void debugWriteIntersectionPoints();
//#endif
//	//------ For debug <<<<<

private:
	// Constructer
	CoastLineList();

	// Destructer
	~CoastLineList();

	// Copy constructer
	CoastLineList(const CoastLineList& rhs);

	// Assignment operator
	CoastLineList& operator=(const CoastLineList& rhs);

	// Total number of coast lines
	int m_numCoastLines;

	// List of coast lines
	CoastLine* m_coastLines;

	// Maximum edge length at outer-most edges
	double m_maxEdgeLengthAtOuterEdges;

};

#endif
