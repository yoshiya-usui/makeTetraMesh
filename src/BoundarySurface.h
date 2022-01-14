#ifndef DBLDEF_BOUNDARY_SURFACE
#define DBLDEF_BOUNDARY_SURFACE

#include <vector>
#include "CommonParameters.h"

// Class of boundary surface
class BoundarySurface{

public:

	enum GeologicalType{
		UNKNOWN = -1,
		SEA = 0,
		LAND,
		ANOMALY,
	};

	// Default constructer
	explicit BoundarySurface();

	// Destructer
	~BoundarySurface();

	// Make boundary surface by delaunay triagulation
	void makeBoundarySurface();

private:

	// Copy constructer
	explicit BoundarySurface(const BoundarySurface& rhs);

	// Assignment operator
	BoundarySurface& operator=(BoundarySurface& rhs);

	// Triangles of boundary surface
	std::vector<int> m_triangleIDs;

	// Geological type of boundary
	int m_geologicalType;


};

#endif
