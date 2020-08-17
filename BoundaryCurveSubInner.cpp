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
#include "BoundaryCurveSubInner.h"
#include "NodeList.h"
#include "OutputFiles.h"
#include <fstream>
#include <iostream>
#include <algorithm>

// Constructer
BoundaryCurveSubInner::BoundaryCurveSubInner( const std::vector<int>& coords, const int itype ):
	BoundaryCurve( coords, itype )
{
}

// Constructer
BoundaryCurveSubInner::BoundaryCurveSubInner( const int itype ):
	BoundaryCurve( itype )
{
}

// Default constructer
BoundaryCurveSubInner::BoundaryCurveSubInner():
	BoundaryCurve()
{
}

// Destructer
BoundaryCurveSubInner::~BoundaryCurveSubInner(){
}

// Copy constructer
BoundaryCurveSubInner::BoundaryCurveSubInner(const BoundaryCurveSubInner& rhs):
	BoundaryCurve( rhs )
{
}

// Assignment operator
BoundaryCurveSubInner& BoundaryCurveSubInner::operator=(BoundaryCurveSubInner& rhs){

	OutputFiles::m_logFile << "Error : " << __FUNCTION__  << " has not been implemented."<< std::endl;
	exit(1);

}

// Get flag specifing whether object locate within the boundary from 
// number of segment on which the object locate left hand side or right hand side
bool BoundaryCurveSubInner::include( const int iLeft, const int iRight ) const{

	if( iLeft > iRight ){
		return true;
	}

	return false;
	
}
