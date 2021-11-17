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
#ifndef DBLDEF_NODE_2D_LIST
#define DBLDEF_NODE_2D_LIST

#include "CommonParameters.h"
#include "Node.h"
#include <vector>
#include <string>
#include <set>

// Class of the list of lake
class NodeList{

public:
	//// Return the the instance of the class
 //   static NodeList* getInstance();

	// Constructer
	NodeList();

	// Destructer
	~NodeList();

	// Clear list of nodes
	void clearList();

	// Add new node and return the ID of the node
	int addNewNode( const Node& node );

	// Get total number of nodes
	int getTotalNumberOfNode() const;

	// Set coordinates of node
	void setCoordOfPoints( const int id, const CommonParameters::XY& coord );

	// Get coordinates of node
	CommonParameters::XY getCoordXYOfPoints( const int id ) const;

	// Get coordinates of node
	CommonParameters::XYZ getCoordXYZOfPoints( const int id ) const;

	// Get pointer to the node of the boundary
	Node* getPointerToNode( const int id );

	// Make all node to be fixed
	void makeAllNodeFixed();

	// Write node list of 2D mesh to VTK file
	void writeNode2DListToVTK( const std::string& fileName ) const;
	
	// Write node list of 2D mesh to intermediate file after step 1
	void writeNode2DListToIntermediateFileStep1( const std::string& fileName ) const;

	// Write node list of 2D mesh without the ones of supert triangle to intermediate file after step 2
	void writeNode2DListToIntermediateFileStep2( const std::string& fileName ) const;

	// Write node list of mesh without the ones of supert triangle to intermediate file after step 3
	void writeNodeListToIntermediateFileStep3( const std::string& fileName ) const;

	// Read node list of 2D mesh from intermediate file
	void readNode2DListFromIntermediateFile( const std::string& fileName );

	// Set node list of 2D mesh
	void setNode2DList( const int numNodes, const CommonParameters::XY* coord2D );

	// Read node list of mesh from intermediate file after step 3
	void readNodeListFromIntermediateFileStep3( const std::string& fileName );

	// Add IDs of triangles sharing the specified node
	void addTriangleIDs( const int nodeID, const int triangleID );

	// Get IDs of triangles sharing the specified node
	std::vector<int> getTriangleIDs( const int nodeID ) const;

	// Initialize node list
	void initializeNodeList();

	// Remove specified node
	void removeNode( const int nodeID );

private:
	// Copy constructer
	NodeList(const NodeList& rhs);

	// Assignment operator
	NodeList& operator=(const NodeList& rhs);

	// List of coast lines
	std::vector<Node> m_nodes;

};

#endif
