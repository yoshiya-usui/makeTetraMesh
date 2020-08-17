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
#ifndef DBLDEF_BOUNDARY_CURVE_OUTER
#define DBLDEF_BOUNDARY_CURVE_OUTER

#include "BoundaryCurve.h"
#include "CommonParameters.h"

// Class of outer boundary curve
class BoundaryCurveOuter : public BoundaryCurve{

public:
	// Constructer
	explicit BoundaryCurveOuter( const std::vector<int>& coords, const int itype );

	// Default constructer
	BoundaryCurveOuter();

	// Constructer
	explicit BoundaryCurveOuter( const int itype );

	// Destructer
	virtual ~BoundaryCurveOuter();

	// Copy constructer
	BoundaryCurveOuter(const BoundaryCurveOuter& rhs);

	// Assignment operator
	BoundaryCurveOuter& operator=(BoundaryCurveOuter& rhs);
	
	// Get flag specifing whether object locate within the boundary from 
	// number of segment on which the object locate left hand side or right hand side
	virtual bool include( const int iLeft, const int iRight ) const;


};

#endif
