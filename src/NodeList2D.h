#ifndef DBLDEF_NODE_LIST_2D
#define DBLDEF_NODE_LIST_2D

#include <vector>
#include "Node2D.h"

// Class of node of two dimensinal field
class NodeList2D{

public:
	// Constructer
	NodeList2D();

	// Destructer
	~NodeList2D();

private:
	// Copy constructer
	NodeList2D(const NodeList2D& rhs);

	// Assignment operator
	NodeList2D& operator=(const NodeList2D& rhs);

	// IDs of nodes of two dimensional field
	std::vector<Node2D> m_nodes;

};

#endif
