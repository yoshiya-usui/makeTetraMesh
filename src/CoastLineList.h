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

	// Output coast line data to txt file for GMT
	void writeCoastLineDataToText( const std::string& fineName ) const;

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
