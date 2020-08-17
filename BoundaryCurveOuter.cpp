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
#include "BoundaryCurveOuter.h"
#include "NodeList.h"
#include "OutputFiles.h"
#include <fstream>
#include <iostream>
#include <algorithm>

// Constructer
BoundaryCurveOuter::BoundaryCurveOuter( const std::vector<int>& coords, const int itype ):
	BoundaryCurve( coords, itype )
{
}

// Constructer
BoundaryCurveOuter::BoundaryCurveOuter( const int itype ):
	BoundaryCurve( itype )
{
}

// Default constructer
BoundaryCurveOuter::BoundaryCurveOuter():
	BoundaryCurve()
{
}

// Destructer
BoundaryCurveOuter::~BoundaryCurveOuter(){
}

// Copy constructer
BoundaryCurveOuter::BoundaryCurveOuter(const BoundaryCurveOuter& rhs):
	BoundaryCurve( rhs )
{
}

// Assignment operator
BoundaryCurveOuter& BoundaryCurveOuter::operator=(BoundaryCurveOuter& rhs){

	OutputFiles::m_logFile << "Error : " << __FUNCTION__  << " has not been implemented."<< std::endl;
	exit(1);

}

// Get flag specifing whether object locate within the boundary from 
// number of segment on which the object locate left hand side or right hand side
bool BoundaryCurveOuter::include( const int iLeft, const int iRight ) const{

	if( iRight > iLeft ){
		return true;
	}

	return false;
	
}
