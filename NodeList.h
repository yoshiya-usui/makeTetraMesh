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
