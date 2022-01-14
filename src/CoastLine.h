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
#ifndef DBLDEF_COAST_LINE
#define DBLDEF_COAST_LINE

#include "CommonParameters.h"
#include "Node.h"
#include <vector>
#include <fstream>

// Class of the coast line
class CoastLine{

public:
	// Default Constructer
	CoastLine();

	// Constructer
	CoastLine( const bool isClosed );

	// Destructer
	~CoastLine();

	// Read data of coast lines from input file
	void readCoastLineData( std::ifstream& ifs );

	// Get total number of points
	int getNumOfPoints() const;

	// Get coordinate X of points
	double getCoordXOfPoints( const int id ) const;

	// Get coordinate Y of points
	double getCoordYOfPoints( const int id ) const;

	// Get coordinate of points
	CommonParameters::XY getCoordXYOfPoints( const int id ) const;

	// Get first coordinate of points
	CommonParameters::XY getFirstCoordOfPoints() const;

	// Get last coordinate of points
	CommonParameters::XY getLastCoordOfPoints() const;

	// Get node data
	Node getNode( const int id ) const;

	// Get flags specifing whether the point must be included or not
	bool isFix( const int id ) const;

	// Get flag specifing the coast line is closed or not
	bool isClosed() const;

	// Remove the nodes consisting to too sharp angle and the nodes whose distances are too short
	void coastLineSmoothing( std::vector<Node>& stack );

	// Roughen coast line
	void roughenCoastLine();

private:
	// Copy constructer
	CoastLine(const CoastLine& rhs);

	// Assignment operator
	CoastLine& operator=(CoastLine& rhs);

	// Flag specifing the coast line is closed or not
	bool m_isClosed;

	// Coordinates of nodes
	std::vector< Node > m_nodes;

	// Flags specifing whether the point must be included or not
	std::vector< bool > m_fix;

};

#endif
