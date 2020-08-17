#ifndef DBLDEF_LAKE
#define DBLDEF_LAKE

#include "CommonParameters.h"
#include <vector>;

// Class of lake
class Lake{

public:
	// Constructer
	Lake();

	// Destructer
	~Lake();

private:
	// Copy constructer
	Lake(const Lake& rhs);

	// Assignment operator
	Lake& operator=(Lake& rhs);

	// IDs of triangles share the node
	std::vector< CommonParameters::XY > m_nodes;

	// Flag specifing the coast line is closed or not
	const static bool m_isClosed = true;

};

#endif
