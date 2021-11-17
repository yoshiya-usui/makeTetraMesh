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
#include "BoundaryCurveInner.h"
#include "NodeList.h"
#include "OutputFiles.h"
#include <fstream>
#include <iostream>
#include <algorithm>

// Constructer
BoundaryCurveInner::BoundaryCurveInner( const std::vector<int>& coords, const int itype ):
	BoundaryCurve( coords, itype )
{
}

// Constructer
BoundaryCurveInner::BoundaryCurveInner( const int itype ):
	BoundaryCurve( itype )
{
}

// Default constructer
BoundaryCurveInner::BoundaryCurveInner():
	BoundaryCurve()
{
}

// Destructer
BoundaryCurveInner::~BoundaryCurveInner(){
}

// Copy constructer
BoundaryCurveInner::BoundaryCurveInner(const BoundaryCurveInner& rhs):
	BoundaryCurve( rhs )
{
}

// Assignment operator
BoundaryCurveInner& BoundaryCurveInner::operator=(BoundaryCurveInner& rhs){

	OutputFiles::m_logFile << "Error : " << __FUNCTION__  << " has not been implemented."<< std::endl;
	exit(1);

}


// Get flag specifing whether object locate within the boundary from 
// number of segment on which the object locate left hand side or right hand side
bool BoundaryCurveInner::include( const int iLeft, const int iRight ) const{

	if( iLeft > iRight ){
		return true;
	}

	return false;
	
}
