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
#ifndef DBLDEF_COMMON_PARAMETERS
#define DBLDEF_COMMON_PARAMETERS

namespace CommonParameters{

struct XY{
	double X;
	double Y;
};

struct XYZ{
	double X;
	double Y;
	double Z;
};

enum DomainType{
	UNKNOWN = -1,
	SEA = 0,
	LAND,
	LAKE,
	AIR,
#ifdef _MOD_FOR_NMT
	NMT_DIPOLE,
#endif
	OUTSIDE_OF_DOMAIN
};

enum Boundary{
	UNDEFINED_BOUNDARY = -1,
	SURFACE = 0,
	TOP,
	BOT,
	YZ_MINUS,
	YZ_PLUS,
	ZX_MINUS,
	ZX_PLUS,
	LAYER,
};

// Circular constant
const static double PI = 3.14159265359;

// Factor converting values from radians to degrees
const static double RAD2DEG = 180.0 / PI;

// Factor converting values from degrees to radians
const static double DEG2RAD = PI / 180.0;

static char programName[]="makeTetraMesh";

static char version[]="v1.3";

}

#endif
