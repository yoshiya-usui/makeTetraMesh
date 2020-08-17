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
