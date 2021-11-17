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
#ifndef DBLDEF_LAKE_LIST
#define DBLDEF_LAKE_LIST

#include "CommonParameters.h"
#include <vector>
#include <map>

// Class of the list of lake
class LakeList{

public:
	// Return the the instance of the class
    static LakeList* getInstance();

	// Read data of lakes from input file
	void readLakeData();

	// Get total number of lake surface points
	int getNumLakeSurfacePoints() const;

	// Get coordinate of lake surface points
	CommonParameters::XY getPointCoord( const int index ) const;

	// Get lake height
	double getLakeHeight( const int index ) const;

private:
	// Constructer
	LakeList();

	// Destructer
	~LakeList();

	// List of the coordintes of lake surface points
	std::vector<CommonParameters::XYZ> m_coordLakeSurfacePoints;

};

#endif