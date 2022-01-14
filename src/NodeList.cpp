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
#include "NodeList.h"
#include "OutputFiles.h"
#include <iostream>
#include <iomanip>
#include <stdlib.h>

//// Return the instance of the class
//NodeList* NodeList::getInstance(){
//   	static NodeList instance;// The only instance
//  	return &instance;
//}

// Default constructer
NodeList::NodeList()
{

}

// Destructer
NodeList::~NodeList(){

}

// Clear list of nodes
void NodeList::clearList(){

	m_nodes.clear();

}

// Add new node and return the ID of the node
int NodeList::addNewNode( const Node& node ){

	m_nodes.push_back( node );

	return static_cast<int>( m_nodes.size() ) - 1;

}

// Get total number of nodes
int NodeList::getTotalNumberOfNode() const{

	return static_cast<int>( m_nodes.size() );

}

// Set coordinates of node
void NodeList::setCoordOfPoints( const int id, const CommonParameters::XY& coord ){

	if( id < 0 || id >= static_cast<int>( m_nodes.size() ) ){
		OutputFiles::m_logFile << "Error : ID of point is out of range at the function " << __FUNCTION__  << " !! ID = " << id << std::endl;
		exit(1);
	}

	m_nodes[id].setCoordXY( coord );

}

// Get coordinates of node
CommonParameters::XY NodeList::getCoordXYOfPoints( const int id ) const{

	if( id < 0 || id >= static_cast<int>( m_nodes.size() ) ){
		//OutputFiles::m_logFile << "Error : ID of point is out of range !! ID = " << id << std::endl;
		OutputFiles::m_logFile << "Error : ID of point is out of range at the function " << __FUNCTION__  << " !! ID = " << id << std::endl;
		exit(1);
	}

	return m_nodes[id].getCoordXY();
	
}

// Get coordinates of node
CommonParameters::XYZ NodeList::getCoordXYZOfPoints( const int id ) const{

	if( id < 0 || id >= static_cast<int>( m_nodes.size() ) ){
		OutputFiles::m_logFile << "Error : ID of point is out of range at the function " << __FUNCTION__  << " !! ID = " << id << std::endl;
		exit(1);
	}

	return m_nodes[id].getCoordXYZ();
	
}

// Get pointer to the node of the boundary
Node* NodeList::getPointerToNode( const int id ){

	if( id < 0 || id >= static_cast<int>( m_nodes.size() ) ){
		OutputFiles::m_logFile << "Error : ID of point is out of range at the function " << __FUNCTION__  << " !! ID = " << id << std::endl;
		exit(1);
	}

	return &m_nodes[id];

}

// Make all node to be fixed
void NodeList::makeAllNodeFixed(){

	std::vector<Node>::iterator itrEnd = m_nodes.end();
	for( std::vector<Node>::iterator itr = m_nodes.begin(); itr != itrEnd; ++itr ){
		itr->fixThisNode();
	}

}

// Write nodes of 2D mesh to VTK file
void NodeList::writeNode2DListToVTK( const std::string& fileName ) const{

	std::ofstream ofsVTK( fileName.c_str() );
	if( !ofsVTK ) {
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << std::endl;
		exit(1);
	}

	ofsVTK << "# vtk DataFile Version 2.0" << std::endl;
	ofsVTK << "NodeList" << std::endl;
	ofsVTK << "ASCII" << std::endl;
	ofsVTK << "DATASET UNSTRUCTURED_GRID" << std::endl;
	
	ofsVTK.precision(9);
	ofsVTK << std::fixed;

	const int numPointsAll = static_cast<int>( m_nodes.size() );

	ofsVTK << "POINTS " << numPointsAll << " float" << std::endl;
	std::vector<Node>::const_iterator itrEnd = m_nodes.end();
	for( std::vector<Node>::const_iterator itr = m_nodes.begin(); itr != itrEnd; ++itr ){
		ofsVTK << itr->getCoordX() << " " << itr->getCoordY() << " " << 0.0 << std::endl;
	}
	
	ofsVTK << "CELLS " << numPointsAll << " " << numPointsAll * 2 << std::endl;
	for( int i = 0; i < numPointsAll; ++i ){
		ofsVTK << 1 << " " << i << std::endl;
	}

	ofsVTK << "CELL_TYPES " << numPointsAll << std::endl;
	for( int i = 0; i < numPointsAll; ++i ){
		ofsVTK << "1" << std::endl;
	}

	ofsVTK << "CELL_DATA " << numPointsAll << std::endl;
	ofsVTK << "SCALARS Location int" << std::endl;
	ofsVTK << "LOOKUP_TABLE default" << std::endl;
	for( std::vector<Node>::const_iterator itr = m_nodes.begin(); itr != itrEnd; ++itr ){
		ofsVTK << itr->getLocation() << std::endl;
	}

	ofsVTK.close();

}

// Write node list of 2D mesh to intermediate file after step 1
void NodeList::writeNode2DListToIntermediateFileStep1( const std::string& fileName ) const{

	std::ofstream ofs( fileName.c_str() );
	if( !ofs.is_open() ) {
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << " !!" << std::endl;
		exit(1);
	}

	const int numNodes = static_cast<int>( m_nodes.size() );

	ofs << std::setw(10) << numNodes << std::endl;

	int icount(0);

	ofs.precision(12);

	for( std::vector<Node>::const_iterator itr = m_nodes.begin(); itr != m_nodes.end(); ++itr ){

		ofs << std::setw(10) << icount++;
		ofs << std::setw(25) << std::scientific << itr->getCoordX()
			<< std::setw(25) << std::scientific << itr->getCoordY()
			<< std::setw(25) << std::scientific << 0.0 << std::endl;

	}

	ofs.close();

}

// Write node list of 2D mesh without the ones of supert triangle to intermediate file after step 2
void NodeList::writeNode2DListToIntermediateFileStep2( const std::string& fileName ) const{

	std::ofstream ofs( fileName.c_str() );
	if( !ofs.is_open() ) {
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << " !!" << std::endl;
		exit(1);
	}

	const int numNodes = static_cast<int>( m_nodes.size() ) - 3;

	ofs << std::setw(10) << numNodes << std::endl;

	int icount(0);
	ofs.precision(12);

	std::vector<Node>::const_iterator itr = m_nodes.begin();
	advance( itr, 3 );
	for( ; itr != m_nodes.end(); ++itr ){

		ofs << std::setw(10) << icount++;
		ofs << std::setw(25) << std::scientific << itr->getCoordX()
			<< std::setw(25) << std::scientific << itr->getCoordY()
			<< std::setw(25) << std::scientific << 0.0 << std::endl;
	}

	ofs.close();

}

// Write node list of mesh without the ones of supert triangle to intermediate file after step 3
void NodeList::writeNodeListToIntermediateFileStep3( const std::string& fileName ) const{

	std::ofstream ofs( fileName.c_str() );
	if( !ofs.is_open() ) {
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << " !!" << std::endl;
		exit(1);
	}

	const int numNodes = static_cast<int>( m_nodes.size() );
	ofs << std::setw(10) << numNodes << std::endl;

	int icount(0);
	ofs.precision(12);

	for( std::vector<Node>::const_iterator itr = m_nodes.begin(); itr != m_nodes.end(); ++itr ){
		ofs << std::setw(10) << icount++;
		ofs << std::setw(25) << std::scientific << itr->getCoordX()
			<< std::setw(25) << std::scientific << itr->getCoordY()
			<< std::setw(25) << std::scientific << itr->getCoordZ()
			<< std::setw(10) << std::scientific << itr->getLocation() << std::endl;
	}

	ofs.close();

}

// Read node list of 2D mesh from intermediate file
void NodeList::readNode2DListFromIntermediateFile( const std::string& fileName ){

	std::ifstream ifs( fileName.c_str() );
	if( !ifs.is_open() ) {
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << " !!" << std::endl;
		exit(1);
	}

	OutputFiles::m_logFile << "# Read node list of 2D mesh from " << fileName.c_str() << std::endl;

	int numNodes(0);

	ifs >> numNodes;

	int ibuf(0);
	double dbuf(0.0);
	for( int iNode = 0; iNode < numNodes; ++iNode ){

		ifs >> ibuf;
		if( ibuf != iNode ){
			OutputFiles::m_logFile << " Error : Node ID must be sequence number from 0 in " << fileName << " !!" << std::endl;
			exit(1);
		}
		
		CommonParameters::XY coord = { 0.0, 0.0 };
		ifs >> coord.X >> coord.Y >> dbuf;

		m_nodes.push_back( Node( coord, true ) );

	}

	ifs.close();

}

// Set node list of 2D mesh
void NodeList::setNode2DList( const int numNodes, const CommonParameters::XY* coord2D ){

	for( int iNode = 0; iNode < numNodes; ++iNode ){
		m_nodes.push_back( Node( coord2D[iNode], true ) );
	}

}

// Read node list of mesh from intermediate file after step 3
void NodeList::readNodeListFromIntermediateFileStep3( const std::string& fileName ){

	std::ifstream ifs( fileName.c_str() );
	if( !ifs.is_open() ) {
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << " !!" << std::endl;
		exit(1);
	}

	OutputFiles::m_logFile << "# Read node list of mesh from " << fileName.c_str() << std::endl;

	int numNodes(0);

	ifs >> numNodes;

	int ibuf(0);
	for( int iNode = 0; iNode < numNodes; ++iNode ){

		ifs >> ibuf;
		if( ibuf != iNode ){
			OutputFiles::m_logFile << " Error : Node ID must be sequence number from 0 in " << fileName << " !!" << std::endl;
			exit(1);
		}
		
		CommonParameters::XYZ coord = { 0.0, 0.0, 0.0 };
		ifs >> coord.X >> coord.Y >> coord.Z;

		int location(CommonParameters::UNKNOWN);
		ifs >> location;

		m_nodes.push_back( Node( coord, true, location ) );

	}

	ifs.close();

}

// Add IDs of triangles sharing the specified node
void NodeList::addTriangleIDs( const int nodeID, const int triangleID ){

	if( nodeID < 0 || nodeID >= static_cast<int>( m_nodes.size() ) ){
		OutputFiles::m_logFile << " Error : Node ID specified to be added ( " << nodeID << " ) is out of range !!" << std::endl;
		exit(1);
	}

	m_nodes[nodeID].addTriangleIDs( triangleID );

}

// Get IDs of triangles sharing the specified node
std::vector<int> NodeList::getTriangleIDs( const int nodeID ) const{

	if( nodeID < 0 || nodeID >= static_cast<int>( m_nodes.size() ) ){
		OutputFiles::m_logFile << " Error : Node ID ( " << nodeID << " ) is out of range !!" << std::endl;
		exit(1);
	}

	return m_nodes[nodeID].getTriangleIDs();

}

// Initialize node list
void NodeList::initializeNodeList(){

	m_nodes.clear();

}

// Remove specified node
void NodeList::removeNode( const int nodeID ){

	if( nodeID < 0 || nodeID >= static_cast<int>( m_nodes.size() ) ){
		OutputFiles::m_logFile << " Error : Node ID specified to be removed ( " << nodeID << " ) is out of range !!" << std::endl;
		exit(1);
	}

	std::vector<Node>::iterator itr = m_nodes.begin();
	advance( itr, nodeID ); 
	m_nodes.erase( itr );

}
