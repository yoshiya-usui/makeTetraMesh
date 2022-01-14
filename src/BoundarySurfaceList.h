#ifndef DBLDEF_BOUNDARY_SURFACE_LIST
#define DBLDEF_BOUNDARY_SURFACE_LIST

#include "BoundarySurface.h"
#include "CommonParameters.h"
#include <vector>

// Class of the list of boundary surface
class BoundarySurfaceList{

public:
	// Return the the instance of the class
    static BoundarySurfaceList* getInstance();

private:
	// Constructer
	BoundarySurfaceList();

	// Destructer
	~BoundarySurfaceList();

	// Copy constructer
	BoundarySurfaceList(const BoundarySurfaceList& rhs);

	// Assignment operator
	BoundarySurfaceList& operator=(BoundarySurfaceList& rhs);

	std::vector<BoundarySurface> m_boundarySurface;

};

#endif
