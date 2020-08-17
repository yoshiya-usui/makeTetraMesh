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
#include "TriangleList.h"
#include "NodeList.h"
#include "OutputFiles.h"
#include "AnalysisDomain.h"
#include "Control.h"
#include "BoundaryCurveList.h"
#include "BoundaryCurveOuter.h"
#include "BoundaryCurveInner.h"
#include "BoundaryCurveSubInner.h"
#include "TopographyDataList.h"
#include "LakeList.h"
#include "Util.h"
#include "math.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <set>
#include <algorithm>
#include <assert.h>

#ifdef _USE_OMP
#include <omp.h>
#endif

// Default constructer
TriangleList::TriangleList():
	m_lastSearchedTriangle(-1),
	m_locationOfPlane(CommonParameters::SURFACE),
	m_indexLayerInterface(-1)
{

}

// Destructer
TriangleList::~TriangleList(){

}


// Create surface triangles
void TriangleList::createSurfaceTriangles(){

	addSuperTriangle();

	makeRoughTriangles();
	std::ostringstream oss0;
	oss0 << "surface_triangle.rough.vtk";
	writeTrianglesToVTK( oss0.str().c_str() );

	assignDomainTypeToTriangles();

	writeTrianglesToVTK( "surface_triangle.rough.domain.vtk" );

	refineTrianglesHavingLongEdges();

	std::ostringstream oss1;
	oss1 << "surface_triangle.refine1_last.vtk";
	writeTrianglesToVTK( oss1.str().c_str() );

	refineTrianglesHavingLargeAngle();

	std::ostringstream oss2;
	oss2 << "surface_triangle.refine2_last.vtk";
	writeTrianglesToVTK( oss2.str().c_str() );

	assignLocationToNodes();

	// Insert node to the triangles more than two nodes of which locate on the coast line
	insertNodeToTriangleWithAllNodesOnCoast();

	std::ostringstream oss3;
	oss3 << "surface_triangle.refine3_last.vtk";
	writeTrianglesToVTK( oss3.str().c_str() );

	laplacian();

	std::ostringstream oss4;
	oss4 << "surface_triangle.laplacian_last.vtk";
	writeTrianglesToVTK( oss4.str().c_str() );
}

// Create surface triangles on additional plane
void TriangleList::createSurfaceTrianglesOnAdditionalPlane(){

	addSuperTriangle();

	makeRoughTriangles();
	std::ostringstream oss0;
	oss0 << "surface_triangle_additional.rough.vtk";
	writeTrianglesToVTK( oss0.str().c_str() );

	assignDomainTypeToTriangles();

	refineTrianglesHavingLongEdges();

	std::ostringstream oss1;
	oss1 << "surface_triangle_additional.refine1_last.vtk";
	writeTrianglesToVTK( oss1.str().c_str() );

	refineTrianglesHavingLargeAngle(true);

	std::ostringstream oss2;
	oss2 << "surface_triangle_additional.refine2_last.vtk";
	writeTrianglesToVTK( oss2.str().c_str() );

	laplacian();

	std::ostringstream oss4;
	oss4 << "surface_triangle_additional.laplacian_last.vtk";
	writeTrianglesToVTK( oss4.str().c_str() );

}

// Add super triangle
void TriangleList::addSuperTriangle(){

	m_triangles.clear();
	m_nodeList.initializeNodeList();

	const AnalysisDomain* const ptrAnalysisDomain = AnalysisDomain::getInstance();

	double maxLength = ( Control::getInstance() )->getFactorSuperTriangleSize();
	double centerX = -9999.0;
	double centerY = -9999.0;

	switch(m_locationOfPlane){
		case CommonParameters::SURFACE:// Go through
		case CommonParameters::TOP:// Go through
		case CommonParameters::BOT:// Go through
		case CommonParameters::LAYER:
			maxLength *= std::max( ptrAnalysisDomain->calcXLength(), ptrAnalysisDomain->calcYLength() );
			centerX = ptrAnalysisDomain->calcCenterCoordX();
			centerY = ptrAnalysisDomain->calcCenterCoordY();
			break;
		case CommonParameters::YZ_MINUS:
			// Go through
		case CommonParameters::YZ_PLUS:
			maxLength *= std::max( ptrAnalysisDomain->calcYLength(), ptrAnalysisDomain->calcZLength() );
			centerX = ptrAnalysisDomain->calcCenterCoordY();
			centerY = ptrAnalysisDomain->calcCenterCoordZ();
			break;
		case CommonParameters::ZX_MINUS:
			// Go through
		case CommonParameters::ZX_PLUS:
			maxLength *= std::max( ptrAnalysisDomain->calcZLength(), ptrAnalysisDomain->calcXLength() );
			centerX = ptrAnalysisDomain->calcCenterCoordZ();
			centerY = ptrAnalysisDomain->calcCenterCoordX();
			break;
		default:
			OutputFiles::m_logFile << " Error : Unknown plane type : " << m_locationOfPlane << std::endl;
			exit(1);
			break;
	}

	const CommonParameters::XY coord0 = { centerX - 1.0 * maxLength , centerY - 1.73 * maxLength }; 
	const CommonParameters::XY coord1 = { centerX + 2.0 * maxLength , centerY }; 
	const CommonParameters::XY coord2 = { centerX - 1.0 * maxLength , centerY + 1.73 * maxLength }; 

	const int nodeIDs[3] ={
		m_nodeList.addNewNode( Node( coord0, false ) ),
		m_nodeList.addNewNode( Node( coord1, false ) ),
		m_nodeList.addNewNode( Node( coord2, false ) )
	};

	const int dummy[3] = { -1, -1, -1 };

	m_triangles.push_back( Triangle( nodeIDs, dummy, dummy, false, CommonParameters::UNKNOWN ) );

	for( int iNode = 0; iNode < 3; ++iNode ){
		m_nodeList.getPointerToNode( nodeIDs[iNode] )->addTriangleIDs(0);
	}

	//----- NOT DELETE FOR FURTURE USE >>>>>
	//const double lengXHalf = ptrAnalysisDomain->calcXLength() * 0.5;
	//const double lengYHalf = ptrAnalysisDomain->calcYLength() * 0.5;
	//const double centerX = ptrAnalysisDomain->calcCenterCoordX();
	//const double centerY = ptrAnalysisDomain->calcCenterCoordY();
	//const double centerZ = ptrAnalysisDomain->calcCenterCoordZ();
	//const double minX = centerX - lengXHalf;
	//const double maxX = centerX + lengXHalf;
	//const double minY = centerY - lengYHalf;
	//const double maxY = centerY + lengYHalf;

	////NodeList* ptrNode2DList = NodeList::getInstance();

	//const CommonParameters::XY coord0 = { minX , minY };
	//const CommonParameters::XY coord1 = { maxX , minY }; 
	//const CommonParameters::XY coord2 = { maxX , maxY }; 
	//const CommonParameters::XY coord3 = { minX , maxY }; 
	//m_nodeList.addNewNode( Node( coord0, false ) );
	//m_nodeList.addNewNode( Node( coord1, false ) );
	//m_nodeList.addNewNode( Node( coord2, false ) );
	//m_nodeList.addNewNode( Node( coord3, false ) );

	//{// Triangle 0
	//	const int nodeIDs[3] ={ 2, 0, 1	};
	//	const int adjTriangles[3] = { 1, -1, -1 };
	//	const int edgeOfAdjTriangles[3] = { 0, -1, -1 };
	//	m_triangles.push_back( Triangle( nodeIDs, adjTriangles, edgeOfAdjTriangles, Triangle::UNKNOWN ) );
	//	for( int iNode = 0; iNode < 3; ++iNode ){
	//		m_nodeList.getPointerToNode( nodeIDs[iNode] )->addTriangleIDs( 0 );
	//	}
	//}

	//{// Triangle 1
	//	const int nodeIDs[3] ={ 0, 2, 3	};
	//	const int adjTriangles[3] = { 0, -1, -1 };
	//	const int edgeOfAdjTriangles[3] = { 0, -1, -1 };
	//	m_triangles.push_back( Triangle( nodeIDs, adjTriangles, edgeOfAdjTriangles, Triangle::UNKNOWN ) );
	//	for( int iNode = 0; iNode < 3; ++iNode ){
	//		m_nodeList.getPointerToNode( nodeIDs[iNode] )->addTriangleIDs( 1 );
	//	}
	//}
	//----- NOT DELETE FOR FURTURE USE <<<<<

}

// Edge flipping
void TriangleList::edgeFlipping( const CommonParameters::XY& coord, std::vector<int>& stack ){

	//const  NodeList* const ptrNode2DList = NodeList::getInstance();

#ifdef _DEBUG_WRITE
	std::cout << "In edgeFlipping" << std::endl;
#endif

	const int countMax = static_cast<int>( m_triangles.size() ) * 10;

#ifdef _DEBUG_WRITE
	std::cout << "countMax : " << countMax << std::endl;
	std::cout << "stack : " << stack.size() << std::endl;
#endif

	int icount(0);
	while(1){

		if( stack.empty() ){
			break;
		}

		if( ++icount > countMax ){
			OutputFiles::m_logFile << " Error : Loop count reach maximum value ( " << countMax << " ) . " << std::endl;
			exit(1);
		}

		const int triangleID = stack.back();

#ifdef _DEBUG_WRITE
	std::cout << "triangleID : " << triangleID << std::endl;
#endif

		const int triangleIDAdj = m_triangles[ triangleID ].getAdjacentTriangleID(1);
		const int edgeIDAdj = m_triangles[ triangleID ].getEdgeIDOfAdjacentTriangle(1);

#ifdef _DEBUG_WRITE
		std::cout << "triangleIDAdj : " << triangleIDAdj << std::endl;
		std::cout << "edgeIDAdj : " << edgeIDAdj << std::endl;
#endif

		if( triangleIDAdj < 0 ){
			stack.pop_back();
			continue;
		}

		const int nodeIDsOrg[3] = {
			m_triangles[triangleID].getNodeID( 0 ),
			m_triangles[triangleID].getNodeID( 1 ),
			m_triangles[triangleID].getNodeID( 2 )
		};

		// For constrained delauny
		std::map<int,int>::iterator itrTemp1 = m_nodeIDTriangle2BoundCurve.find( nodeIDsOrg[1] );
		std::map<int,int>::iterator itrTemp2 = m_nodeIDTriangle2BoundCurve.find( nodeIDsOrg[2] );
		if( itrTemp1 != m_nodeIDTriangle2BoundCurve.end() &&
			itrTemp2 != m_nodeIDTriangle2BoundCurve.end() &&
			m_boundaryCurveList.nodePairOfBoundaryCurve( itrTemp1->second, itrTemp2->second ) ){

			// This edge belongs to the boundary curve
			stack.pop_back();
			continue;
		}

		const int adjTriangleIDsOrg[3] = {
			m_triangles[triangleID].getAdjacentTriangleID( 0 ),
			m_triangles[triangleID].getAdjacentTriangleID( 1 ),
			m_triangles[triangleID].getAdjacentTriangleID( 2 )
		};

		//const int edgeIDsOfAdjTriangleOrg[3] = {
		//	m_triangles[triangleID].getEdgeIDOfAdjacentTriangle( 0 ),
		//	m_triangles[triangleID].getEdgeIDOfAdjacentTriangle( 1 ),
		//	m_triangles[triangleID].getEdgeIDOfAdjacentTriangle( 2 )
		//};

		const int nodeIDsAdjOrg[3] = {
			m_triangles[triangleIDAdj].getNodeID( 0 ),
			m_triangles[triangleIDAdj].getNodeID( 1 ),
			m_triangles[triangleIDAdj].getNodeID( 2 )
		};

		const int adjTriangleIDsAdjOrg[3] = {
			m_triangles[triangleIDAdj].getAdjacentTriangleID( 0 ),
			m_triangles[triangleIDAdj].getAdjacentTriangleID( 1 ),
			m_triangles[triangleIDAdj].getAdjacentTriangleID( 2 )
		};

		const int edgeIDsOfAdjTriangleAdjOrg[3] = {
			m_triangles[triangleIDAdj].getEdgeIDOfAdjacentTriangle( 0 ),
			m_triangles[triangleIDAdj].getEdgeIDOfAdjacentTriangle( 1 ),
			m_triangles[triangleIDAdj].getEdgeIDOfAdjacentTriangle( 2 )
		};

#ifdef _DEBUG_WRITE
		std::cout << "Original triangle first" << std::endl;
		std::cout << "triangleID : " << triangleID << std::endl;
		std::cout << "NodeID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getNodeID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "AdjacentTriangleID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getAdjacentTriangleID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "EdgeIDOfAdjacentTriangle : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getEdgeIDOfAdjacentTriangle(iNode) << " ";
		}
		std::cout << std::endl;

		std::cout << "Original triangle second" << std::endl;
		std::cout << "triangleID : " << triangleIDAdj << std::endl;
		std::cout << "NodeID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleIDAdj].getNodeID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "AdjacentTriangleID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleIDAdj].getAdjacentTriangleID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "EdgeIDOfAdjacentTriangle : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleIDAdj].getEdgeIDOfAdjacentTriangle(iNode) << " ";
		}
		std::cout << std::endl;
#endif

		const CommonParameters::XY coord0 = m_nodeList.getCoordXYOfPoints( nodeIDsAdjOrg[0] );
		const CommonParameters::XY coord1 = m_nodeList.getCoordXYOfPoints( nodeIDsAdjOrg[1] );
		const CommonParameters::XY coord2 = m_nodeList.getCoordXYOfPoints( nodeIDsAdjOrg[2] );

		if( Util::inCircle( coord0, coord1, coord2, coord ) ){
			// The inserted node belongs to the circumcircule of adjacent triangle

			//---------------------------------------------------------------------
			m_triangles[triangleID].setNodeID( 2, nodeIDsAdjOrg[( edgeIDAdj + 2 ) % 3]  );

			m_triangles[triangleID].setAdjacentTriangleID( 1, adjTriangleIDsAdjOrg[( edgeIDAdj + 1 ) % 3] );
			m_triangles[triangleID].setAdjacentTriangleID( 2, triangleIDAdj );

			m_triangles[triangleID].setEdgeIDOfAdjacentTriangle( 1, edgeIDsOfAdjTriangleAdjOrg[( edgeIDAdj + 1 ) % 3] );
			m_triangles[triangleID].setEdgeIDOfAdjacentTriangle( 2, 0 );
			//---------------------------------------------------------------------

			//---------------------------------------------------------------------
			m_triangles[triangleIDAdj].setNodeID( 0, nodeIDsOrg[0] );
			m_triangles[triangleIDAdj].setNodeID( 1, nodeIDsAdjOrg[ ( edgeIDAdj + 2 ) % 3 ]  );
			m_triangles[triangleIDAdj].setNodeID( 2, nodeIDsOrg[2] );

			m_triangles[triangleIDAdj].setAdjacentTriangleID( 0, triangleID );
			m_triangles[triangleIDAdj].setAdjacentTriangleID( 1, adjTriangleIDsAdjOrg[( edgeIDAdj + 2 ) % 3] );
			m_triangles[triangleIDAdj].setAdjacentTriangleID( 2, adjTriangleIDsOrg[2] );

			m_triangles[triangleIDAdj].setEdgeIDOfAdjacentTriangle( 0, 2 );
			m_triangles[triangleIDAdj].setEdgeIDOfAdjacentTriangle( 1, edgeIDsOfAdjTriangleAdjOrg[( edgeIDAdj + 2 ) % 3] );
			m_triangles[triangleIDAdj].setEdgeIDOfAdjacentTriangle( 2, 0 );
			//---------------------------------------------------------------------

			//---------------------------------------------------------------------
			m_triangles[adjTriangleIDsOrg[0]].setAdjacentTriangleID( 2, triangleID );
			m_triangles[adjTriangleIDsOrg[0]].setEdgeIDOfAdjacentTriangle( 2, 0 );
			//---------------------------------------------------------------------

			//---------------------------------------------------------------------
			m_triangles[adjTriangleIDsOrg[2]].setAdjacentTriangleID( 0, triangleIDAdj );
			m_triangles[adjTriangleIDsOrg[2]].setEdgeIDOfAdjacentTriangle( 0, 2 );
			//---------------------------------------------------------------------

			//---------------------------------------------------------------------
			if( adjTriangleIDsAdjOrg[(edgeIDAdj+1) % 3] >= 0 ){
				m_triangles[adjTriangleIDsAdjOrg[(edgeIDAdj+1) % 3]].setAdjacentTriangleID( edgeIDsOfAdjTriangleAdjOrg[(edgeIDAdj+1) % 3], triangleID );
				m_triangles[adjTriangleIDsAdjOrg[(edgeIDAdj+1) % 3]].setEdgeIDOfAdjacentTriangle( edgeIDsOfAdjTriangleAdjOrg[(edgeIDAdj+1) % 3], 1 );
			}
			//---------------------------------------------------------------------

			//---------------------------------------------------------------------
			if( adjTriangleIDsAdjOrg[(edgeIDAdj+2) % 3] >= 0 ){
				m_triangles[adjTriangleIDsAdjOrg[(edgeIDAdj+2) % 3]].setAdjacentTriangleID( edgeIDsOfAdjTriangleAdjOrg[(edgeIDAdj+2) % 3], triangleIDAdj );
				m_triangles[adjTriangleIDsAdjOrg[(edgeIDAdj+2) % 3]].setEdgeIDOfAdjacentTriangle( edgeIDsOfAdjTriangleAdjOrg[(edgeIDAdj+2) % 3], 1 );
			}
			//---------------------------------------------------------------------

#ifdef _DEBUG_WRITE
			std::cout << "Modfied triangle first" << std::endl;
			std::cout << "triangleID : " << triangleID << std::endl;
			std::cout << "NodeID : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[triangleID].getNodeID(iNode) << " ";
			}
			std::cout << std::endl;
			std::cout << "AdjacentTriangleID : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[triangleID].getAdjacentTriangleID(iNode) << " ";
			}
			std::cout << std::endl;
			std::cout << "EdgeIDOfAdjacentTriangle : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[triangleID].getEdgeIDOfAdjacentTriangle(iNode) << " ";
			}
			std::cout << std::endl;

			std::cout << "Modfied triangle second" << std::endl;
			std::cout << "triangleID : " << triangleIDAdj << std::endl;
			std::cout << "NodeID : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[triangleIDAdj].getNodeID(iNode) << " ";
			}
			std::cout << std::endl;
			std::cout << "AdjacentTriangleID : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[triangleIDAdj].getAdjacentTriangleID(iNode) << " ";
			}
			std::cout << std::endl;
			std::cout << "EdgeIDOfAdjacentTriangle : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[triangleIDAdj].getEdgeIDOfAdjacentTriangle(iNode) << " ";
			}
			std::cout << std::endl;
#endif

			//---------------------------------------------------------------------
			m_nodeList.getPointerToNode( nodeIDsOrg[0] )->addTriangleIDs( triangleIDAdj );
			m_nodeList.getPointerToNode( nodeIDsOrg[1] )->eraseTriangleIDs( triangleIDAdj );
			m_nodeList.getPointerToNode( nodeIDsOrg[2] )->eraseTriangleIDs( triangleID );
			m_nodeList.getPointerToNode( nodeIDsAdjOrg[( edgeIDAdj + 2 ) % 3] )->addTriangleIDs( triangleID );
			//---------------------------------------------------------------------
			
#ifdef _DEBUG_WRITE
			for( int i = 0; i < 3; ++i ){
				std::cout << "nodeIDsOrg[" << i << "] : " << nodeIDsOrg[i] << " ";
				std::vector<int> tempVec = m_nodeList.getPointerToNode( nodeIDsOrg[i] )->getTriangleIDs();
				for( std::vector<int>::iterator itr = tempVec.begin(); itr != tempVec.end(); ++itr ){
					std::cout << *itr << " ";
				}
				std::cout << std::endl;
			}
			std::cout << "nodeIDsAdjOrg[( edgeIDAdj + 2 ) % 3] : " << nodeIDsAdjOrg[( edgeIDAdj + 2 ) % 3] << " ";
			std::vector<int> tempVec = m_nodeList.getPointerToNode( nodeIDsAdjOrg[( edgeIDAdj + 2 ) % 3] )->getTriangleIDs();
			for( std::vector<int>::iterator itr = tempVec.begin(); itr != tempVec.end(); ++itr ){
				std::cout << *itr << " ";
			}
			std::cout << std::endl;
#endif

			if( std::find( stack.begin(), stack.end(), triangleIDAdj ) == stack.end() ){
				stack.push_back( triangleIDAdj );
			}

		}
		else{
			// The inserted node locate outside of circumcircule of adjacent triangle
			stack.pop_back();
		}


	}

}

// Search triangle enclosing the specified point
int TriangleList::searchTriangleEnclosingPoint( const CommonParameters::XY& coord, int& onSide, int& locType ){

#ifdef _DEBUG_WRITE
	std::cout << "In searchTriangleEnclosingPoint" << std::endl;
#endif

	const int numTriangles = static_cast<int>( m_triangles.size() );
	const int maxCount = numTriangles;

#ifdef _DEBUG_WRITE
	std::cout << "numTriangles = " << numTriangles << std::endl;
	std::cout << "maxCount = " << maxCount << std::endl;
#endif

	if( m_lastSearchedTriangle < 0 ){
		m_lastSearchedTriangle = 0;
	}
	
	int iTriangle = m_lastSearchedTriangle;

#ifdef _DEBUG_WRITE
	std::cout << "m_lastSearchedTriangle = " << m_lastSearchedTriangle << std::endl;
#endif

	int icount(0);
	while(1){

		++icount;

		//if( icount  > maxCount ){
		//	OutputFiles::m_logFile << " Error : Loop count reach maximum value " << maxCount << " in searching triangle enclosing point ( " << coord.X << ", " << coord.Y << " )" << std::endl;
		//	exit(1);
		//}
		if( icount  > maxCount ){
#ifdef _DEBUG_WRITE
			std::cout << "Loop count reach maximum value " << maxCount << " in searching triangle enclosing point ( " << coord.X << ", " << coord.Y << " )" << std::endl;
#endif
			break;
		}

		//if( icount % numTriangles ){
		//	//srand((unsigned int)time(NULL));
		//	iTriangle = rand() % numTriangles;
		//}

		onSide = -1;
		locType = m_triangles[iTriangle].doesLocateInTriangle( coord, onSide, &m_nodeList );
		switch( locType ){
			case Triangle::OUTSIDE_TRIANGLE:
				iTriangle = onSide;// Update next triangle
				if( iTriangle < 0 ){
					iTriangle = rand() % numTriangles;
				}
				break;
			case Triangle::ON_NODE:
				OutputFiles::m_logFile << " Warning : Specified point (X,Y) = ( " << coord.X << ", " << coord.Y << " ) is a duplicate of node " << onSide << " of triangle " << iTriangle << "." << std::endl;
				return m_triangles[iTriangle].getNodeID(onSide);
				break;
			case Triangle::INSIDE_TRIANGLE:
				onSide = -1;
				return iTriangle; 
				break;
			case Triangle::ON_EDGE:
				return iTriangle; 
				break;
			default:
				OutputFiles::m_logFile << " Error : Unsupported flag of result of searching procedure." << std::endl;
				exit(1);
				break;
		}


	}

	for( int iTriangle = 0; iTriangle < numTriangles; ++iTriangle ){

#ifdef _DEBUG_WRITE
		std::cout << "iTriangle = " << iTriangle << std::endl;
#endif

		onSide = -1;
		locType = m_triangles[iTriangle].doesLocateInTriangle( coord, onSide, &m_nodeList );
		switch( locType ){
			case Triangle::OUTSIDE_TRIANGLE:
				break;
			case Triangle::ON_NODE:
				OutputFiles::m_logFile << " Warning : Specified point (X,Y) = ( " << coord.X << ", " << coord.Y << " ) is a duplicate of node " << onSide << " of triangle " << iTriangle << "." << std::endl;
				return m_triangles[iTriangle].getNodeID(onSide);
				break;
			case Triangle::INSIDE_TRIANGLE:
				onSide = -1;
				return iTriangle; 
				break;
			case Triangle::ON_EDGE:
				return iTriangle; 
				break;
			default:
				OutputFiles::m_logFile << " Error : Unsupported flag of result of searching procedure." << std::endl;
				exit(1);
				break;
		}

	}

	OutputFiles::m_logFile << " Error : Triangle enclosing the secified point (X,Y) = ( " << coord.X << ", " << coord.Y << " ) is not found." << std::endl;
	exit(1);

	return -1;

}

// Insert new node and swap
int TriangleList::insertNewNodeAndFlip( const CommonParameters::XY& coord, int& locType, const int nodeIDBoundCurveList ){

	//***********************
	//***** Insert node *****
	//***********************

	int onSide(0);
	const int triangleID = searchTriangleEnclosingPoint( coord, onSide, locType );

	if( locType == Triangle::ON_NODE ){
		return triangleID;
	}
	else if( triangleID < 0 ){
		return -1;
	}

#ifdef _DEBUG_WRITE
	std::cout << "triangleID : " << triangleID << std::endl;
	std::cout << "locType : " << locType << std::endl;
#endif
	
	//const int iDomainType = m_triangles[triangleID].getDomainType();

	// Insert node to list
	const int nodeIDInserted = m_nodeList.addNewNode( Node( coord, false ) );
	if( nodeIDBoundCurveList >= 0 ){
		m_nodeIDTriangle2BoundCurve.insert( std::make_pair( nodeIDInserted, nodeIDBoundCurveList ) );
		m_nodeIDBoundCurve2Triangle.insert( std::make_pair( nodeIDBoundCurveList, nodeIDInserted ) );
	}

#ifdef _DEBUG_WRITE
	std::cout << "      nodeIDInserted : " << nodeIDInserted << std::endl;
	std::cout << "nodeIDBoundCurveList : " << nodeIDBoundCurveList << std::endl;
#endif

	const int triangleIDNew1 = static_cast<int>( m_triangles.size() );
	const int triangleIDNew2 = triangleIDNew1 + 1;

#ifdef _DEBUG_WRITE
	std::cout << "triangleIDNew1 : " << triangleIDNew1 << std::endl;
	std::cout << "triangleIDNew2 : " << triangleIDNew2 << std::endl;
#endif

	const int nodeIDsOrg[3] = {
		m_triangles[triangleID].getNodeID( 0 ),
		m_triangles[triangleID].getNodeID( 1 ),
		m_triangles[triangleID].getNodeID( 2 )
	};

#ifdef _DEBUG_WRITE
	std::cout << "nodeIDsOrg : ";
	for( int iNode = 0; iNode < 3; ++iNode ){
		std::cout << nodeIDsOrg[iNode] << " ";
	}
	std::cout << std::endl;
#endif

	const int adjTriangleIDsOrg[3] = {
		m_triangles[triangleID].getAdjacentTriangleID( 0 ),
		m_triangles[triangleID].getAdjacentTriangleID( 1 ),
		m_triangles[triangleID].getAdjacentTriangleID( 2 )
	};

#ifdef _DEBUG_WRITE
	std::cout << "adjTriangleIDsOrg : ";
	for( int iNode = 0; iNode < 3; ++iNode ){
		std::cout << adjTriangleIDsOrg[iNode] << " ";
	}
	std::cout << std::endl;
#endif

	const int edgeIDsOfAdjTriangleOrg[3] = {
		m_triangles[triangleID].getEdgeIDOfAdjacentTriangle( 0 ),
		m_triangles[triangleID].getEdgeIDOfAdjacentTriangle( 1 ),
		m_triangles[triangleID].getEdgeIDOfAdjacentTriangle( 2 )
	};

#ifdef _DEBUG_WRITE
	std::cout << "edgeIDsOfAdjTriangleOrg : ";
	for( int iNode = 0; iNode < 3; ++iNode ){
		std::cout << edgeIDsOfAdjTriangleOrg[iNode] << " ";
	}
	std::cout << std::endl;
	std::cout << "onSide = " << onSide << std::endl;
#endif

	const int domainTypeID = m_triangles[triangleID].getDomainType();

	std::vector<int> stack;// Stack of IDs of triangle to be checked

	if( onSide < 0 ){// Locate inside of triangle

#ifdef _DEBUG_WRITE
		std::cout << "Locate inside of triangle" << std::endl;
#endif

		//---------------------------------------------------------------------
		// For triangle containing the node
		//---------------------------------------------------------------------
		m_triangles[triangleID].setNodeID( 0, nodeIDInserted );
		m_triangles[triangleID].setNodeID( 1, nodeIDsOrg[0]  );
		m_triangles[triangleID].setNodeID( 2, nodeIDsOrg[1]  );

		m_triangles[triangleID].setAdjacentTriangleID( 0, triangleIDNew2 );
		m_triangles[triangleID].setAdjacentTriangleID( 1, adjTriangleIDsOrg[0] );
		m_triangles[triangleID].setAdjacentTriangleID( 2, triangleIDNew1 );

		m_triangles[triangleID].setEdgeIDOfAdjacentTriangle( 0, 2 );
		m_triangles[triangleID].setEdgeIDOfAdjacentTriangle( 1, edgeIDsOfAdjTriangleOrg[0] );
		m_triangles[triangleID].setEdgeIDOfAdjacentTriangle( 2, 0 );

		m_triangles[triangleID].setSatisfyCriteria( false );

#ifdef _DEBUG_WRITE
		std::cout << "Original triangle" << std::endl;
		std::cout << "triangleID : " << triangleID << std::endl;
		std::cout << "NodeID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getNodeID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "AdjacentTriangleID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getAdjacentTriangleID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "EdgeIDOfAdjacentTriangle : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getEdgeIDOfAdjacentTriangle(iNode) << " ";
		}
		std::cout << std::endl;
#endif
		//---------------------------------------------------------------------

		//---------------------------------------------------------------------
		// For new triangle ( first )
		//---------------------------------------------------------------------
		int nodeIDsNew[3] = {
			nodeIDInserted,
			nodeIDsOrg[1],
			nodeIDsOrg[2]
		};
		int adjTriangleIDsNew[3] = { 
			triangleID,
			adjTriangleIDsOrg[1],
			triangleIDNew2
		};
		int edgeIDsOfAdjTriangleNew[3] = { 
			2,
			edgeIDsOfAdjTriangleOrg[1],
			0
		};
		m_triangles.push_back( Triangle( nodeIDsNew, adjTriangleIDsNew, edgeIDsOfAdjTriangleNew, false, domainTypeID ) );

#ifdef _DEBUG_WRITE
		std::cout << "new triangle ( first )" << std::endl;
		std::cout << "triangleID : " << m_triangles.size() - 1 << std::endl;
		std::cout << "NodeID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles.back().getNodeID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "AdjacentTriangleID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles.back().getAdjacentTriangleID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "EdgeIDOfAdjacentTriangle : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles.back().getEdgeIDOfAdjacentTriangle(iNode) << " ";
		}
		std::cout << std::endl;
#endif

		//---------------------------------------------------------------------

		//---------------------------------------------------------------------
		// For new triangle ( second )
		//---------------------------------------------------------------------
		nodeIDsNew[0] = nodeIDInserted;
		nodeIDsNew[1] = nodeIDsOrg[2];
		nodeIDsNew[2] = nodeIDsOrg[0];
		adjTriangleIDsNew[0] = triangleIDNew1;
		adjTriangleIDsNew[1] = adjTriangleIDsOrg[2];
		adjTriangleIDsNew[2] = triangleID;
		edgeIDsOfAdjTriangleNew[0] = 2;
		edgeIDsOfAdjTriangleNew[1] = edgeIDsOfAdjTriangleOrg[2];
		edgeIDsOfAdjTriangleNew[2] = 0;
		m_triangles.push_back( Triangle( nodeIDsNew, adjTriangleIDsNew, edgeIDsOfAdjTriangleNew, false, domainTypeID ) );

#ifdef _DEBUG_WRITE
		std::cout << "new triangle ( second )" << std::endl;
		std::cout << "triangleID : " << m_triangles.size() - 1 << std::endl;
		std::cout << "NodeID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles.back().getNodeID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "AdjacentTriangleID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles.back().getAdjacentTriangleID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "EdgeIDOfAdjacentTriangle : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles.back().getEdgeIDOfAdjacentTriangle(iNode) << " ";
		}
		std::cout << std::endl;
#endif
		//---------------------------------------------------------------------

		//---------------------------------------------------------------------
		// For adjacent triangle ( first )
		//---------------------------------------------------------------------
		if( adjTriangleIDsOrg[0] >= 0 ){
			m_triangles[ adjTriangleIDsOrg[0] ].setAdjacentTriangleID( edgeIDsOfAdjTriangleOrg[0], triangleID );
			m_triangles[ adjTriangleIDsOrg[0] ].setEdgeIDOfAdjacentTriangle( edgeIDsOfAdjTriangleOrg[0], 1 );

#ifdef _DEBUG_WRITE
			std::cout << "adjacent triangle ( first )" << std::endl;
			std::cout << "triangleID : " << adjTriangleIDsOrg[0] << std::endl;
			std::cout << "NodeID : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsOrg[0]].getNodeID(iNode) << " ";
			}
			std::cout << std::endl;
			std::cout << "AdjacentTriangleID : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsOrg[0]].getAdjacentTriangleID(iNode) << " ";
			}
			std::cout << std::endl;
			std::cout << "EdgeIDOfAdjacentTriangle : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsOrg[0]].getEdgeIDOfAdjacentTriangle(iNode) << " ";
			}
			std::cout << std::endl;
#endif

		}
		//---------------------------------------------------------------------

		//---------------------------------------------------------------------
		// For adjacent triangle ( second )
		//---------------------------------------------------------------------
		if( adjTriangleIDsOrg[1] >= 0 ){
			m_triangles[ adjTriangleIDsOrg[1] ].setAdjacentTriangleID( edgeIDsOfAdjTriangleOrg[1], triangleIDNew1 );
			m_triangles[ adjTriangleIDsOrg[1] ].setEdgeIDOfAdjacentTriangle( edgeIDsOfAdjTriangleOrg[1], 1 );

#ifdef _DEBUG_WRITE
			std::cout << "adjacent triangle ( second )" << std::endl;
			std::cout << "triangleID : " << adjTriangleIDsOrg[1] << std::endl;
			std::cout << "NodeID : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsOrg[1]].getNodeID(iNode) << " ";
			}
			std::cout << std::endl;
			std::cout << "AdjacentTriangleID : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsOrg[1]].getAdjacentTriangleID(iNode) << " ";
			}
			std::cout << std::endl;
			std::cout << "EdgeIDOfAdjacentTriangle : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsOrg[1]].getEdgeIDOfAdjacentTriangle(iNode) << " ";
			}
			std::cout << std::endl;
#endif

		}
		//---------------------------------------------------------------------

		//---------------------------------------------------------------------
		// For adjacent triangle ( third )
		//---------------------------------------------------------------------
		if( adjTriangleIDsOrg[2] >= 0 ){
			m_triangles[ adjTriangleIDsOrg[2] ].setAdjacentTriangleID( edgeIDsOfAdjTriangleOrg[2], triangleIDNew2 );
			m_triangles[ adjTriangleIDsOrg[2] ].setEdgeIDOfAdjacentTriangle( edgeIDsOfAdjTriangleOrg[2], 1 );

#ifdef _DEBUG_WRITE
			std::cout << "adjacent triangle ( thrid )" << std::endl;
			std::cout << "triangleID : " << adjTriangleIDsOrg[2] << std::endl;
			std::cout << "NodeID : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsOrg[2]].getNodeID(iNode) << " ";
			}
			std::cout << std::endl;
			std::cout << "AdjacentTriangleID : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsOrg[2]].getAdjacentTriangleID(iNode) << " ";
			}
			std::cout << std::endl;
			std::cout << "EdgeIDOfAdjacentTriangle : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsOrg[2]].getEdgeIDOfAdjacentTriangle(iNode) << " ";
			}
			std::cout << std::endl;
#endif

		}
		//---------------------------------------------------------------------
		
		stack.push_back( triangleID );
		stack.push_back( triangleIDNew1 );
		stack.push_back( triangleIDNew2 );

		//---------------------------------------------------------------------
		// Array relating node to elements
		//---------------------------------------------------------------------
		m_nodeList.getPointerToNode( nodeIDInserted )->addTriangleIDs( triangleID );
		m_nodeList.getPointerToNode( nodeIDInserted )->addTriangleIDs( triangleIDNew1 );
		m_nodeList.getPointerToNode( nodeIDInserted )->addTriangleIDs( triangleIDNew2 );

		m_nodeList.getPointerToNode( nodeIDsOrg[0] )->addTriangleIDs( triangleIDNew2 );

		m_nodeList.getPointerToNode( nodeIDsOrg[1] )->addTriangleIDs( triangleIDNew1 );

		m_nodeList.getPointerToNode( nodeIDsOrg[2] )->eraseTriangleIDs( triangleID );
		m_nodeList.getPointerToNode( nodeIDsOrg[2] )->addTriangleIDs( triangleIDNew1 );
		m_nodeList.getPointerToNode( nodeIDsOrg[2] )->addTriangleIDs( triangleIDNew2 );

#ifdef _DEBUG_WRITE
		std::cout << "nodeIDInserted : " << nodeIDInserted << " ";
		std::vector<int> tempVec = m_nodeList.getPointerToNode( nodeIDInserted )->getTriangleIDs();
		for( std::vector<int>::iterator itr = tempVec.begin(); itr != tempVec.end(); ++itr ){
			std::cout << *itr << " ";
		}
		std::cout << std::endl;

		for( int i = 0; i < 3; ++i ){
			std::cout << "nodeIDsOrg[" << i << "] : " << nodeIDsOrg[i] << " ";
			std::vector<int> tempVec = m_nodeList.getPointerToNode( nodeIDsOrg[i] )->getTriangleIDs();
			for( std::vector<int>::iterator itr = tempVec.begin(); itr != tempVec.end(); ++itr ){
				std::cout << *itr << " ";
			}
			std::cout << std::endl;
		}
#endif

	}
	else{// Locate on an edge of triangle

#ifdef _DEBUG_WRITE
		std::cout << "Locate on an edge of triangle" << std::endl;
#endif

		const int triangleIDAdj = adjTriangleIDsOrg[onSide];
		const int onSideAdj = edgeIDsOfAdjTriangleOrg[onSide];
		
		const int nodeIDsAdjOrg[3] = {
			m_triangles[triangleIDAdj].getNodeID( 0 ),
			m_triangles[triangleIDAdj].getNodeID( 1 ),
			m_triangles[triangleIDAdj].getNodeID( 2 )
		};

		const int adjTriangleIDsAdjOrg[3] = {
			m_triangles[triangleIDAdj].getAdjacentTriangleID( 0 ),
			m_triangles[triangleIDAdj].getAdjacentTriangleID( 1 ),
			m_triangles[triangleIDAdj].getAdjacentTriangleID( 2 )
		};

		const int edgeIDsOfAdjTriangleAdjOrg[3] = {
			m_triangles[triangleIDAdj].getEdgeIDOfAdjacentTriangle( 0 ),
			m_triangles[triangleIDAdj].getEdgeIDOfAdjacentTriangle( 1 ),
			m_triangles[triangleIDAdj].getEdgeIDOfAdjacentTriangle( 2 )
		};

		const int domainTypeIDAdj = m_triangles[triangleIDAdj].getDomainType();

#ifdef _DEBUG_WRITE
		std::cout << "Original adjacent triangle" << std::endl;
		std::cout << "triangleIDAdj : " << triangleIDAdj << std::endl;
		std::cout << "NodeID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleIDAdj].getNodeID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "AdjacentTriangleID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleIDAdj].getAdjacentTriangleID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "EdgeIDOfAdjacentTriangle : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleIDAdj].getEdgeIDOfAdjacentTriangle(iNode) << " ";
		}
		std::cout << std::endl;
#endif

		//---------------------------------------------------------------------
		// For triangle containing the node ( first )
		//---------------------------------------------------------------------
		m_triangles[triangleID].setNodeID( 0, nodeIDInserted );
		m_triangles[triangleID].setNodeID( 1, nodeIDsOrg[(onSide+2)%3]  );
		m_triangles[triangleID].setNodeID( 2, nodeIDsOrg[onSide%3]  );

		m_triangles[triangleID].setAdjacentTriangleID( 0, triangleIDNew1 );
		m_triangles[triangleID].setAdjacentTriangleID( 1, adjTriangleIDsOrg[(onSide+2)%3] );
		m_triangles[triangleID].setAdjacentTriangleID( 2, triangleIDAdj );

		m_triangles[triangleID].setEdgeIDOfAdjacentTriangle( 0, 2 );
		m_triangles[triangleID].setEdgeIDOfAdjacentTriangle( 1, edgeIDsOfAdjTriangleOrg[(onSide+2)%3] );
		m_triangles[triangleID].setEdgeIDOfAdjacentTriangle( 2, 0 );

		m_triangles[triangleID].setSatisfyCriteria( false );

#ifdef _DEBUG_WRITE
		std::cout << "Modified triangle" << std::endl;
		std::cout << "triangleID : " << triangleID << std::endl;
		std::cout << "NodeID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getNodeID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "AdjacentTriangleID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getAdjacentTriangleID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "EdgeIDOfAdjacentTriangle : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getEdgeIDOfAdjacentTriangle(iNode) << " ";
		}
		std::cout << std::endl;
#endif
		//---------------------------------------------------------------------

		//---------------------------------------------------------------------
		// For triangle containing the node ( second )
		//---------------------------------------------------------------------
		m_triangles[triangleIDAdj].setNodeID( 0, nodeIDInserted );
		m_triangles[triangleIDAdj].setNodeID( 1, nodeIDsAdjOrg[(onSideAdj+1)%3] );
		m_triangles[triangleIDAdj].setNodeID( 2, nodeIDsAdjOrg[(onSideAdj+2)%3] );

		m_triangles[triangleIDAdj].setAdjacentTriangleID( 0, triangleID );
		m_triangles[triangleIDAdj].setAdjacentTriangleID( 1, adjTriangleIDsAdjOrg[(onSideAdj+1)%3] );
		m_triangles[triangleIDAdj].setAdjacentTriangleID( 2, triangleIDNew2 );

		m_triangles[triangleIDAdj].setEdgeIDOfAdjacentTriangle( 0, 2 );
		m_triangles[triangleIDAdj].setEdgeIDOfAdjacentTriangle( 1, edgeIDsOfAdjTriangleAdjOrg[(onSideAdj+1)%3] );
		m_triangles[triangleIDAdj].setEdgeIDOfAdjacentTriangle( 2, 0 );

		m_triangles[triangleIDAdj].setSatisfyCriteria( false );

#ifdef _DEBUG_WRITE
		std::cout << "Modified adjacent triangle" << std::endl;
		std::cout << "triangleIDAdj : " << triangleIDAdj << std::endl;
		std::cout << "NodeID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleIDAdj].getNodeID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "AdjacentTriangleID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleIDAdj].getAdjacentTriangleID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "EdgeIDOfAdjacentTriangle : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleIDAdj].getEdgeIDOfAdjacentTriangle(iNode) << " ";
		}
		std::cout << std::endl;
#endif
		//---------------------------------------------------------------------

		//---------------------------------------------------------------------
		// For new triangle ( first )
		//---------------------------------------------------------------------
		int nodeIDsNew[3] = {
			nodeIDInserted,
			nodeIDsOrg[(onSide+1)%3],
			nodeIDsOrg[(onSide+2)%3]
		};
		int adjTriangleIDsNew[3] = { 
			triangleIDNew2,
			adjTriangleIDsOrg[(onSide+1)%3],
			triangleID
		};
		int edgeIDsOfAdjTriangleNew[3] = { 
			2,
			edgeIDsOfAdjTriangleOrg[(onSide+1)%3],
			0
		};
		m_triangles.push_back( Triangle( nodeIDsNew, adjTriangleIDsNew, edgeIDsOfAdjTriangleNew, false, domainTypeID ) );

#ifdef _DEBUG_WRITE
		std::cout << "new triangle ( first )" << std::endl;
		std::cout << "triangleID : " << m_triangles.size() - 1 << std::endl;
		std::cout << "NodeID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles.back().getNodeID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "AdjacentTriangleID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles.back().getAdjacentTriangleID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "EdgeIDOfAdjacentTriangle : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles.back().getEdgeIDOfAdjacentTriangle(iNode) << " ";
		}
		std::cout << std::endl;
#endif
		//---------------------------------------------------------------------

		//---------------------------------------------------------------------
		// For new triangle ( second )
		//---------------------------------------------------------------------
		nodeIDsNew[0] = nodeIDInserted;
		nodeIDsNew[1] = nodeIDsAdjOrg[(onSideAdj+2)%3];
		nodeIDsNew[2] = nodeIDsAdjOrg[onSideAdj%3];
		adjTriangleIDsNew[0] = triangleIDAdj;
		adjTriangleIDsNew[1] = adjTriangleIDsAdjOrg[(onSideAdj+2)%3];
		adjTriangleIDsNew[2] = triangleIDNew1;
		edgeIDsOfAdjTriangleNew[0] = 2;
		edgeIDsOfAdjTriangleNew[1] = edgeIDsOfAdjTriangleAdjOrg[(onSideAdj+2)%3];
		edgeIDsOfAdjTriangleNew[2] = 0;
		m_triangles.push_back( Triangle( nodeIDsNew, adjTriangleIDsNew, edgeIDsOfAdjTriangleNew, false, domainTypeIDAdj ) );

#ifdef _DEBUG_WRITE
		std::cout << "new triangle ( second )" << std::endl;
		std::cout << "triangleID : " << m_triangles.size() - 1 << std::endl;
		std::cout << "NodeID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles.back().getNodeID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "AdjacentTriangleID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles.back().getAdjacentTriangleID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "EdgeIDOfAdjacentTriangle : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles.back().getEdgeIDOfAdjacentTriangle(iNode) << " ";
		}
		std::cout << std::endl;
#endif
		//---------------------------------------------------------------------

		//---------------------------------------------------------------------
		// For adjacent triangle ( first )
		//---------------------------------------------------------------------
		if( adjTriangleIDsOrg[(onSide+2)%3] >= 0 ){
			m_triangles[ adjTriangleIDsOrg[(onSide+2)%3] ].setAdjacentTriangleID( edgeIDsOfAdjTriangleOrg[(onSide+2)%3], triangleID );
			m_triangles[ adjTriangleIDsOrg[(onSide+2)%3] ].setEdgeIDOfAdjacentTriangle( edgeIDsOfAdjTriangleOrg[(onSide+2)%3], 1 );

#ifdef _DEBUG_WRITE
			std::cout << "adjacent triangle ( first )" << std::endl;
			std::cout << "triangleID : " << adjTriangleIDsOrg[(onSide+2)%3] << std::endl;
			std::cout << "NodeID : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsOrg[(onSide+2)%3]].getNodeID(iNode) << " ";
			}
			std::cout << std::endl;
			std::cout << "AdjacentTriangleID : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsOrg[(onSide+2)%3]].getAdjacentTriangleID(iNode) << " ";
			}
			std::cout << std::endl;
			std::cout << "EdgeIDOfAdjacentTriangle : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsOrg[(onSide+2)%3]].getEdgeIDOfAdjacentTriangle(iNode) << " ";
			}
			std::cout << std::endl;
#endif
		}
		//---------------------------------------------------------------------

		//---------------------------------------------------------------------
		// For adjacent triangle ( second )
		//---------------------------------------------------------------------
		if( adjTriangleIDsAdjOrg[(onSideAdj+1)%3] >= 0 ){
			m_triangles[ adjTriangleIDsAdjOrg[(onSideAdj+1)%3] ].setAdjacentTriangleID( edgeIDsOfAdjTriangleAdjOrg[(onSideAdj+1)%3], triangleIDAdj );
			m_triangles[ adjTriangleIDsAdjOrg[(onSideAdj+1)%3] ].setEdgeIDOfAdjacentTriangle( edgeIDsOfAdjTriangleAdjOrg[(onSideAdj+1)%3], 1 );


#ifdef _DEBUG_WRITE
			std::cout << "adjacent triangle ( second )" << std::endl;
			std::cout << "triangleID : " << adjTriangleIDsAdjOrg[(onSideAdj+1)%3] << std::endl;
			std::cout << "NodeID : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsAdjOrg[(onSideAdj+1)%3]].getNodeID(iNode) << " ";
			}
			std::cout << std::endl;
			std::cout << "AdjacentTriangleID : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsAdjOrg[(onSideAdj+1)%3]].getAdjacentTriangleID(iNode) << " ";
			}
			std::cout << std::endl;
			std::cout << "EdgeIDOfAdjacentTriangle : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsAdjOrg[(onSideAdj+1)%3]].getEdgeIDOfAdjacentTriangle(iNode) << " ";
			}
			std::cout << std::endl;
#endif
		}
		//---------------------------------------------------------------------

		//---------------------------------------------------------------------
		// For adjacent triangle ( third )
		//---------------------------------------------------------------------
		if( adjTriangleIDsOrg[(onSide+1)%3] >= 0 ){
			m_triangles[ adjTriangleIDsOrg[(onSide+1)%3] ].setAdjacentTriangleID( edgeIDsOfAdjTriangleOrg[(onSide+1)%3], triangleIDNew1 );
			m_triangles[ adjTriangleIDsOrg[(onSide+1)%3] ].setEdgeIDOfAdjacentTriangle( edgeIDsOfAdjTriangleOrg[(onSide+1)%3], 1 );

#ifdef _DEBUG_WRITE
			std::cout << "adjacent triangle ( third )" << std::endl;
			std::cout << "triangleID : " << adjTriangleIDsOrg[(onSide+1)%3] << std::endl;
			std::cout << "NodeID : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsOrg[(onSide+1)%3]].getNodeID(iNode) << " ";
			}
			std::cout << std::endl;
			std::cout << "AdjacentTriangleID : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsOrg[(onSide+1)%3]].getAdjacentTriangleID(iNode) << " ";
			}
			std::cout << std::endl;
			std::cout << "EdgeIDOfAdjacentTriangle : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsOrg[(onSide+1)%3]].getEdgeIDOfAdjacentTriangle(iNode) << " ";
			}
			std::cout << std::endl;
#endif
		}
		//---------------------------------------------------------------------

		//---------------------------------------------------------------------
		// For adjacent triangle ( fourth )
		//---------------------------------------------------------------------
		if( adjTriangleIDsAdjOrg[(onSideAdj+2)%3] >= 0 ){
			m_triangles[ adjTriangleIDsAdjOrg[(onSideAdj+2)%3] ].setAdjacentTriangleID( edgeIDsOfAdjTriangleAdjOrg[(onSideAdj+2)%3], triangleIDNew2 );
			m_triangles[ adjTriangleIDsAdjOrg[(onSideAdj+2)%3] ].setEdgeIDOfAdjacentTriangle( edgeIDsOfAdjTriangleAdjOrg[(onSideAdj+2)%3], 1 );

#ifdef _DEBUG_WRITE
			std::cout << "adjacent triangle ( fourth )" << std::endl;
			std::cout << "triangleID : " << adjTriangleIDsAdjOrg[(onSideAdj+2)%3]  << std::endl;
			std::cout << "NodeID : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsAdjOrg[(onSideAdj+2)%3] ].getNodeID(iNode) << " ";
			}
			std::cout << std::endl;
			std::cout << "AdjacentTriangleID : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsAdjOrg[(onSideAdj+2)%3] ].getAdjacentTriangleID(iNode) << " ";
			}
			std::cout << std::endl;
			std::cout << "EdgeIDOfAdjacentTriangle : ";
			for( int iNode = 0; iNode < 3; ++iNode ){
				std::cout << m_triangles[adjTriangleIDsAdjOrg[(onSideAdj+2)%3] ].getEdgeIDOfAdjacentTriangle(iNode) << " ";
			}
			std::cout << std::endl;
#endif
		}
		//---------------------------------------------------------------------

		stack.push_back( triangleID );
		stack.push_back( triangleIDAdj );
		stack.push_back( triangleIDNew1 );
		stack.push_back( triangleIDNew2 );

		//---------------------------------------------------------------------
		// Array relating node to elements
		//---------------------------------------------------------------------
		m_nodeList.getPointerToNode( nodeIDInserted )->addTriangleIDs( triangleID );
		m_nodeList.getPointerToNode( nodeIDInserted )->addTriangleIDs( triangleIDAdj );
		m_nodeList.getPointerToNode( nodeIDInserted )->addTriangleIDs( triangleIDNew1 );
		m_nodeList.getPointerToNode( nodeIDInserted )->addTriangleIDs( triangleIDNew2 );

		m_nodeList.getPointerToNode( nodeIDsOrg[(onSide+1)%3] )->eraseTriangleIDs( triangleID );
		m_nodeList.getPointerToNode( nodeIDsOrg[(onSide+1)%3] )->eraseTriangleIDs( triangleIDAdj );
		m_nodeList.getPointerToNode( nodeIDsOrg[(onSide+1)%3] )->addTriangleIDs( triangleIDNew1 );
		m_nodeList.getPointerToNode( nodeIDsOrg[(onSide+1)%3] )->addTriangleIDs( triangleIDNew2 );

		m_nodeList.getPointerToNode( nodeIDsOrg[(onSide+2)%3] )->addTriangleIDs( triangleIDNew1 );

		m_nodeList.getPointerToNode( nodeIDsAdjOrg[(onSideAdj+2)%3] )->addTriangleIDs( triangleIDNew2 );

#ifdef _DEBUG_WRITE
		std::cout << "nodeIDInserted : " << nodeIDInserted << " ";
		std::vector<int> tempVec = m_nodeList.getPointerToNode( nodeIDInserted )->getTriangleIDs();
		for( std::vector<int>::iterator itr = tempVec.begin(); itr != tempVec.end(); ++itr ){
			std::cout << *itr << " ";
		}
		std::cout << std::endl;

		for( int i = 0; i < 3; ++i ){
			std::cout << "nodeIDsOrg[" << (onSide+i)%3 << "] : " << nodeIDsOrg[(onSide+i)%3] << " ";
			std::vector<int> tempVec = m_nodeList.getPointerToNode( nodeIDsOrg[(onSide+i)%3] )->getTriangleIDs();
			for( std::vector<int>::iterator itr = tempVec.begin(); itr != tempVec.end(); ++itr ){
				std::cout << *itr << " ";
			}
			std::cout << std::endl;
		}
		std::cout << "nodeIDsAdjOrg[(onSideAdj+2)%3] : " << nodeIDsAdjOrg[(onSideAdj+2)%3] << " ";
		tempVec = m_nodeList.getPointerToNode( nodeIDsAdjOrg[(onSideAdj+2)%3] )->getTriangleIDs();
		for( std::vector<int>::iterator itr = tempVec.begin(); itr != tempVec.end(); ++itr ){
			std::cout << *itr << " ";
		}
		std::cout << std::endl;
#endif

	}

	//*************************
	//***** Edge flipping *****
	//*************************
	//if( static_cast<int>( stack.size() ) > 2 ){
	//	edgeFlipping( coord, stack );
	//	if( !stack.empty() ){
	//		OutputFiles::m_logFile << " Error : Stack is not empty !!" << std::endl;
	//		exit(1);
	//	}
	//}
	edgeFlipping( coord, stack );
	if( !stack.empty() ){
		OutputFiles::m_logFile << " Error : Stack is not empty !!" << std::endl;
		exit(1);
	}

	return nodeIDInserted;

}

// Search for triangles between specified two nodes
void TriangleList::searchElemBetweenNodes( const int node0, const int node1, std::vector<int>& triangleIDs ) const{

	//const  NodeList* const ptrNode2DList = NodeList::getInstance();

#ifdef _DEBUG_WRITE
	std::cout << "node0 node1 : " << node0 << " " << node1 << std::endl;
#endif

	const CommonParameters::XY coord0 = m_nodeList.getCoordXYOfPoints( node0 );
	const CommonParameters::XY coord1 = m_nodeList.getCoordXYOfPoints( node1 );

	// Element including the first node
	std::vector<int> trianglesContainingNode0 = m_nodeList.getTriangleIDs( node0 );
	
	int triangleID(-1);
	int edgeID(-1);

#ifdef _DEBUG_WRITE
	std::cout << "trianglesContainingNode0.size() : " << trianglesContainingNode0.size() << std::endl;
#endif

	// Element including the first node which have edge intersect with the segment connecting specified two nodes
	for( std::vector<int>::iterator itr = trianglesContainingNode0.begin(); itr != trianglesContainingNode0.end(); ++itr ){

		triangleID = *itr;

#ifdef _DEBUG_WRITE
		std::cout << "triangleID : " << triangleID << std::endl;
#endif
		//int onSide(-1);
		//if( m_triangles[triangleID].doesLocateInTriangle( coord1, onSide, &m_nodeList ) ){
		//	return;
		//}

		int icount(0);
		int otherNodes[2] = { -1, -1 };
		for( int iNode = 0; iNode < 3; ++iNode  ){
			const int nodeID = m_triangles[triangleID].getNodeID( iNode );

#ifdef _DEBUG_WRITE
			std::cout << "nodeID : " << nodeID << std::endl;
#endif

			if( nodeID == node0 ){
				continue;
			}
			else if( nodeID == node1 ){
				return;
			}
			otherNodes[icount++] = iNode;
		}

		if( icount != 2 ){
			OutputFiles::m_logFile << " Error : icount ( " << icount << " ) must be 2 !!" << std::endl;
			exit(1);
		}

		if( Util::intersectTwoSegments( coord0, coord1,
			m_nodeList.getCoordXYOfPoints( m_triangles[triangleID].getNodeID( otherNodes[0] ) ),
			m_nodeList.getCoordXYOfPoints( m_triangles[triangleID].getNodeID( otherNodes[1] ) ) ) ){
			triangleIDs.push_back( triangleID );
			switch( otherNodes[0] + otherNodes[1] ){
				case 1:
					edgeID = 0;
					break;
				case 2:
					edgeID = 2;
					break;
				case 3:
					edgeID = 1;
					break;
				default:
					OutputFiles::m_logFile << " Error : Wrong node ID pair ( " << node0 << ", " << node1 << " )." << std::endl;
					exit(1);
					break;
			}
			break;
		}
		
	}

	if( triangleIDs.empty() ){
		OutputFiles::m_logFile << " Error : No element includes the first node which have edge intersect with the segment connecting node " << node0  << " and " << node1 << "." << std::endl;
		exit(1);
	}

	if( triangleID < -1 ){
		OutputFiles::m_logFile << " Error : Triangle ID is not assigned." << std::endl;
		exit(1);
	}

	if( edgeID < -1 ){
		OutputFiles::m_logFile << " Error : Edge ID is not assigned." << std::endl;
		exit(1);
	}

	const int countMax = static_cast<int>( m_triangles.size() );

	int icount(0);
	while(1){

		if( ++icount > countMax ){
			OutputFiles::m_logFile << " Error : Loop count reach maximum value ( " << countMax << " ) . " << std::endl;
			exit(1);
		}

		const int edgeIDPre = m_triangles[triangleID].getEdgeIDOfAdjacentTriangle( edgeID );
		triangleID = m_triangles[triangleID].getAdjacentTriangleID( edgeID );

		if( m_triangles[triangleID].getNodeID( (edgeIDPre+2) % 3 ) == node1 ){
			triangleIDs.push_back( triangleID );
			return;
		}
		
		for( int iEdge = 0; iEdge < 3; ++iEdge ){
			if( iEdge == edgeIDPre ){
				continue;
			}

			const int nodeIDs[2] = { iEdge, (iEdge+1) % 3 };

			const int nodeIDGlobal[2] = { m_triangles[triangleID].getNodeID( nodeIDs[0] ), m_triangles[triangleID].getNodeID( nodeIDs[1] ) };

			if( nodeIDGlobal[0] == node1 || nodeIDGlobal[1] == node1 ){
				return;
			}

			if( Util::intersectTwoSegments( coord0, coord1,	m_nodeList.getCoordXYOfPoints( nodeIDGlobal[0] ), m_nodeList.getCoordXYOfPoints( nodeIDGlobal[1] ) ) ){
				triangleIDs.push_back( triangleID );
				switch( nodeIDs[0] + nodeIDs[1] ){
					case 1:
						edgeID = 0;
						break;
					case 2:
						edgeID = 2;
						break;
					case 3:
						edgeID = 1;
						break;
					default:
						OutputFiles::m_logFile << " Error : Wrong node ID pair ( " << node0 << ", " << node1 << " )." << std::endl;
						exit(1);
						break;
				}
			}
			
		}

	}

}

// Subdivide polygon constructed from specified triangles 
void TriangleList::subDividePolygon( const int node0, const int node1, const std::vector<int>& triangleIDs ){

	//const CommonParameters::XY coord0 = m_nodeList.getCoordXYOfPoints( node0 );
	//const CommonParameters::XY coord1 = m_nodeList.getCoordXYOfPoints( node1 );

	std::vector< std::pair< std::pair< int, int >, std::pair< int, int > > > segmentVec;

	int triangleID(0);
	int iNode(0);
	for( std::vector<int>::const_iterator itr = triangleIDs.begin(); itr != triangleIDs.end(); ++itr ){

		triangleID = *itr;

		for( iNode = 0; iNode < 3; ++iNode  ){
			const int nodeID = m_triangles[triangleID].getNodeID( iNode );
			if( node0 == nodeID ){
				break;
			}
		}
		
	}

	const int edgeID = iNode;
	const int triangleIDAdj =  m_triangles[triangleID].getAdjacentTriangleID( edgeID );

	if( find( triangleIDs.begin(), triangleIDs.end(), triangleIDAdj ) != triangleIDs.end() ){
		// Adjacent element belongs to polygon
		triangleID = triangleIDAdj;
		const int edgeIDAdj = m_triangles[triangleID].getEdgeIDOfAdjacentTriangle( edgeID );
		iNode = ( edgeIDAdj + 1 ) % 3;
	}
	else{
		// Adjacent element does not belong to polygon
		const std::pair<int,int> nodePair( m_triangles[triangleID].getNodeID( iNode ), m_triangles[triangleID].getNodeID( (iNode+1)%3 ) );
		const std::pair<int,int> elemAndEdge( triangleID, edgeID );
		segmentVec.push_back( std::make_pair( nodePair, elemAndEdge ) );
	}


	//const int countMax = static_cast<int>( triangleIDs.size() );

	//int icount(0);

	//while(1){

	//	if( ++icount > countMax ){
	//		OutputFiles::m_logFile << " Error : Loop count reach maximum value ( " << countMax << " ) . " << std::endl;
	//		exit(1);
	//	}




	//}

}

// Write triangles to VTK
void TriangleList::writeTrianglesToVTK( const std::string& fileName ){

	// Open output vtk file -----
	std::ofstream ofsVTK( fileName.c_str() );
	if( !ofsVTK ) {
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << std::endl;
		exit(1);
	}

	ofsVTK << "# vtk DataFile Version 2.0" << std::endl;
	ofsVTK << "SurfaceTriangles" << std::endl;
	ofsVTK << "ASCII" << std::endl;
	ofsVTK << "DATASET UNSTRUCTURED_GRID" << std::endl;
	
	ofsVTK.precision(6);
	ofsVTK << std::fixed;

	//std::set<int> nodesOfTriangles;
	//for( std::vector< Triangle >::iterator itr = m_triangles.begin(); itr != m_triangles.end(); ++itr ){
	//	for( int iNode = 0; iNode < 3; ++iNode ){
	//		nodesOfTriangles.insert( itr->getNodeID(iNode) );
	//	}
	//}

	//NodeList* ptrNode2DList = NodeList::getInstance();

	const int numNodes = m_nodeList.getTotalNumberOfNode();
	ofsVTK << "POINTS " << numNodes<< " float" << std::endl;
	for( int iNode = 0; iNode < numNodes; ++iNode ){
		const CommonParameters::XYZ coord = m_nodeList.getCoordXYZOfPoints( iNode );
		ofsVTK << std::setw(15) << std::scientific << coord.X
			   << std::setw(15) << std::scientific << coord.Y
			   << std::setw(15) << std::scientific << coord.Z << std::endl;
	}

	const int numTriangles = static_cast<int>( m_triangles.size() );
	ofsVTK << "CELLS " << numTriangles << " " << numTriangles * 4 << std::endl;

	for( std::vector< Triangle >::iterator itr = m_triangles.begin(); itr != m_triangles.end(); ++itr ){
		ofsVTK << std::setw(10) << 3;
		for( int iNode = 0; iNode < 3; ++iNode ){
			ofsVTK << std::setw(10) << itr->getNodeID(iNode);
		}
		ofsVTK << std::endl;
	}

	ofsVTK << "CELL_TYPES " << numTriangles << std::endl;
	for( int i = 0; i < numTriangles; ++i ){
		ofsVTK << std::setw(10) << 5 << std::endl;
	}

	ofsVTK << "CELL_DATA " << numTriangles << std::endl;
	ofsVTK << "SCALARS DomainType int" << std::endl;
	ofsVTK << "LOOKUP_TABLE default" << std::endl;
	for( int i = 0; i < numTriangles; ++i ){
		ofsVTK << std::setw(10) << m_triangles[i].getDomainType() << std::endl;
	}

	ofsVTK << "SCALARS TriangleSerial int" << std::endl;
	ofsVTK << "LOOKUP_TABLE default" << std::endl;
	for( int i = 0; i < numTriangles; ++i ){
		ofsVTK << std::setw(10) << i << std::endl;
	}

	ofsVTK << "POINT_DATA " << numNodes << std::endl;
	ofsVTK << "SCALARS Height float" << std::endl;
	ofsVTK << "LOOKUP_TABLE default" << std::endl;
	for( int iNode = 0; iNode < numNodes; ++iNode ){
		ofsVTK << std::setw(15) << std::scientific << static_cast<float>( m_nodeList.getCoordXYZOfPoints(iNode).Z ) << std::endl;
	}

	double* deviation = new double[numNodes];
	calcDeviationOfHeight(deviation);
	ofsVTK << "SCALARS HeightDeviation float" << std::endl;
	ofsVTK << "LOOKUP_TABLE default" << std::endl;
	for( int iNode = 0; iNode < numNodes; ++iNode ){
		ofsVTK << std::setw(15) << std::scientific << static_cast<float>( deviation[iNode] ) << std::endl;
	}
	delete [] deviation;

	ofsVTK << "SCALARS NodeSerial int" << std::endl;
	ofsVTK << "LOOKUP_TABLE default" << std::endl;
	for( int iNode = 0; iNode < numNodes; ++iNode ){
		ofsVTK << std::setw(15) << iNode << std::endl;
	}

	ofsVTK.close();

}

// Write PLCs to VTK
void TriangleList::writePLCsToVTK( const std::string& fileName, const std::vector<CommonParameters::XYZ>& nodeList, const std::vector<TriangleList::Facet>& facetList ){

	// Open output vtk file -----
	std::ofstream ofsVTK( fileName.c_str() );
	if( !ofsVTK ) {
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << std::endl;
		exit(1);
	}

	ofsVTK << "# vtk DataFile Version 2.0" << std::endl;
	ofsVTK << "PLC" << std::endl;
	ofsVTK << "ASCII" << std::endl;
	ofsVTK << "DATASET UNSTRUCTURED_GRID" << std::endl;
	
	ofsVTK.precision(6);
	ofsVTK << std::fixed;

	//std::set<int> nodesOfTriangles;
	//for( std::vector< Triangle >::iterator itr = m_triangles.begin(); itr != m_triangles.end(); ++itr ){
	//	for( int iNode = 0; iNode < 3; ++iNode ){
	//		nodesOfTriangles.insert( itr->getNodeID(iNode) );
	//	}
	//}

	//NodeList* ptrNode2DList = NodeList::getInstance();

	const int numNodesExceptHoles = static_cast<int>( nodeList.size() );
	int numHoles = 0;
	for( std::vector<TriangleList::Facet>::const_iterator itr = facetList.begin(); itr != facetList.end(); ++itr ){
		for( std::vector<CommonParameters::XYZ>::const_iterator itrHoles = itr->holes.begin(); itrHoles != itr->holes.end(); ++itrHoles ){
			++numHoles;
		}
	}

	const int numNodes = numHoles + numNodesExceptHoles;

	ofsVTK << "POINTS " << numNodes<< " float" << std::endl;
	for( std::vector<CommonParameters::XYZ>::const_iterator itr = nodeList.begin(); itr != nodeList.end(); ++itr ){
		ofsVTK << std::setw(15) << std::scientific << itr->X
			   << std::setw(15) << std::scientific << itr->Y
			   << std::setw(15) << std::scientific << itr->Z << std::endl;
	}
	for( std::vector<TriangleList::Facet>::const_iterator itr = facetList.begin(); itr != facetList.end(); ++itr ){
		for( std::vector<CommonParameters::XYZ>::const_iterator itrHoles = itr->holes.begin(); itrHoles != itr->holes.end(); ++itrHoles ){
			ofsVTK << std::setw(15) << std::scientific << itrHoles->X
				   << std::setw(15) << std::scientific << itrHoles->Y
				   << std::setw(15) << std::scientific << itrHoles->Z << std::endl;
		}
	}

	int numPolygons = 0;
	int numData = 0;
	for( std::vector<TriangleList::Facet>::const_iterator itr = facetList.begin(); itr != facetList.end(); ++itr ){
		numPolygons += static_cast<int>( itr->polygons.size() );
		for( std::vector<PolygonNodes>::const_iterator itrPoly = itr->polygons.begin(); itrPoly != itr->polygons.end(); ++itrPoly ){
			numData += static_cast<int>( itrPoly->size() ) + 1;
		}
	}

	ofsVTK << "CELLS " << numPolygons + numHoles << " " << numData + numHoles * 2 << std::endl;
	for( std::vector<TriangleList::Facet>::const_iterator itr = facetList.begin(); itr != facetList.end(); ++itr ){
		for( std::vector<PolygonNodes>::const_iterator itrPoly = itr->polygons.begin(); itrPoly != itr->polygons.end(); ++itrPoly ){
			ofsVTK << std::setw(10) << static_cast<int>( itrPoly->size() );
			for( std::vector<int>::const_iterator itrNode = itrPoly->begin(); itrNode != itrPoly->end(); ++itrNode ){
				ofsVTK << std::setw(10) << *itrNode;
			}
			ofsVTK << std::endl;
		}
	}
	for( int i = numNodesExceptHoles; i < numNodes; ++i ){
		ofsVTK << std::setw(10) << 1 << std::setw(10) << i << std::endl;
	}

	ofsVTK << "CELL_TYPES " << numPolygons + numHoles << std::endl;
	for( int i = 0; i < numPolygons; ++i ){
		ofsVTK << std::setw(10) << 7 << std::endl;
	}
	for( int i = 0; i < numHoles; ++i ){
		ofsVTK << std::setw(10) << 1 << std::endl;
	}

	ofsVTK << "CELL_DATA " << numPolygons + numHoles << std::endl;
	ofsVTK << "SCALARS ID int" << std::endl;
	ofsVTK << "LOOKUP_TABLE default" << std::endl;
	int icount(1);
	for( std::vector<TriangleList::Facet>::const_iterator itr = facetList.begin(); itr != facetList.end(); ++itr ){
		for( std::vector<PolygonNodes>::const_iterator itrPoly = itr->polygons.begin(); itrPoly != itr->polygons.end(); ++itrPoly ){
			ofsVTK << std::setw(10) << icount++ << std::endl;
		}
	}
	for( int i = 0; i < numHoles; ++i ){
		ofsVTK << std::setw(10) << icount++ << std::endl;
	}

	ofsVTK << "POINT_DATA " << numNodes << std::endl;
	ofsVTK << "SCALARS NodeID int" << std::endl;
	ofsVTK << "LOOKUP_TABLE default" << std::endl;
	icount = 1;
	for( int iNode = 0; iNode < numNodes; ++iNode ){
		ofsVTK << std::setw(10) << icount++ << std::endl;
	}


	ofsVTK.close();

}

// Write triangles without the ones belonging to supert triangle to intermediate file after step 2
void TriangleList::writeTrianglesToIntermediateFileStep2( const std::string& fileName ) const{

	OutputFiles::m_logFile << "# Write triangles without super triangle" << std::endl;

	std::ofstream ofs( fileName.c_str() );
	if( !ofs.is_open() ) {
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << " !!" << std::endl;
		exit(1);
	}

	const int numTriangleWithSupperTriangle = static_cast<int>( m_triangles.size() );
	int* convTriangleIDs = new int[ numTriangleWithSupperTriangle ];

	int numTriangles(0);
	for( int iTriangle = 0; iTriangle < numTriangleWithSupperTriangle; ++iTriangle ){
		if( m_triangles[iTriangle].getDomainType() == CommonParameters::OUTSIDE_OF_DOMAIN ||
			m_triangles[iTriangle].getDomainType() == CommonParameters::UNKNOWN ){
			continue;
		}
		convTriangleIDs[iTriangle] = numTriangles++;
	}

	ofs << std::setw(10) << numTriangles << std::endl;

	ofs.precision(9);
	for( int iTriangle = 0; iTriangle < numTriangleWithSupperTriangle; ++iTriangle ){
		if( m_triangles[iTriangle].getDomainType() == CommonParameters::OUTSIDE_OF_DOMAIN ||
			m_triangles[iTriangle].getDomainType() == CommonParameters::UNKNOWN ){
			continue;
		}

		ofs << std::setw(10) << convTriangleIDs[iTriangle];

		for( int i = 0; i < 3; ++i ){
			ofs << std::setw(10) << m_triangles[iTriangle].getNodeID(i) - 3;
		}

		for( int i = 0; i < 3; ++i ){
			const int triangleIDAdj = m_triangles[iTriangle].getAdjacentTriangleID(i);
			if( m_triangles[triangleIDAdj].getDomainType() == CommonParameters::OUTSIDE_OF_DOMAIN ||
				m_triangles[triangleIDAdj].getDomainType() == CommonParameters::UNKNOWN ){
				ofs << std::setw(10) << -1;
			}else{
				ofs << std::setw(10) << convTriangleIDs[ triangleIDAdj ];
			}
		}

		for( int i = 0; i < 3; ++i ){
			const int triangleIDAdj = m_triangles[iTriangle].getAdjacentTriangleID(i);
			if( m_triangles[triangleIDAdj].getDomainType() == CommonParameters::OUTSIDE_OF_DOMAIN ||
				m_triangles[triangleIDAdj].getDomainType() == CommonParameters::UNKNOWN ){
				ofs << std::setw(10) << -1;
			}else{
				ofs << std::setw(10) << m_triangles[iTriangle].getEdgeIDOfAdjacentTriangle(i);
			}
		}

		ofs << std::setw(10) << m_triangles[iTriangle].getDomainType() << std::endl;

	}

	delete [] convTriangleIDs;

	ofs.close();

}

// Write triangles without the ones belonging to supert triangle to intermediate file after step 3
void TriangleList::writeTrianglesToIntermediateFileStep3( const std::string& fileName ) const{

	OutputFiles::m_logFile << "# Write triangles" << std::endl;

	std::ofstream ofs( fileName.c_str() );
	if( !ofs.is_open() ) {
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << " !!" << std::endl;
		exit(1);
	}

	const int numTriangles = static_cast<int>( m_triangles.size() );

	ofs << std::setw(10) << numTriangles << std::endl;

	ofs.precision(9);
	for( int iTriangle = 0; iTriangle < numTriangles; ++iTriangle ){

		ofs << std::setw(10) << iTriangle;

		for( int i = 0; i < 3; ++i ){
			ofs << std::setw(10) << m_triangles[iTriangle].getNodeID(i);
		}

		for( int i = 0; i < 3; ++i ){
			const int triangleIDAdj = m_triangles[iTriangle].getAdjacentTriangleID(i);
			if( triangleIDAdj < 0 || 
				m_triangles[triangleIDAdj].getDomainType() == CommonParameters::OUTSIDE_OF_DOMAIN ||
				m_triangles[triangleIDAdj].getDomainType() == CommonParameters::UNKNOWN ){
				ofs << std::setw(10) << -1;
			}else{
				ofs << std::setw(10) << triangleIDAdj;
			}
		}

		for( int i = 0; i < 3; ++i ){
			const int triangleIDAdj = m_triangles[iTriangle].getAdjacentTriangleID(i);
			if( triangleIDAdj < 0 || 
				m_triangles[triangleIDAdj].getDomainType() == CommonParameters::OUTSIDE_OF_DOMAIN ||
				m_triangles[triangleIDAdj].getDomainType() == CommonParameters::UNKNOWN ){
				ofs << std::setw(10) << -1;
			}else{
				ofs << std::setw(10) << m_triangles[iTriangle].getEdgeIDOfAdjacentTriangle(i);
			}
		}

		ofs << std::setw(10) << m_triangles[iTriangle].getDomainType() << std::endl;

	}

	ofs.close();

}

// Write nodes of triangles without the ones of supert triangle to intermediate file after step 2
void TriangleList::writeNode2DListToIntermediateFileStep2( const std::string& fileName ) const{

	OutputFiles::m_logFile << "# Write nodes" << std::endl;
	m_nodeList.writeNode2DListToIntermediateFileStep2( fileName );

}

// Write nodes of triangles without the ones of supert triangle to intermediate file after step 3
void TriangleList::writeNodeListToIntermediateFileStep3( const std::string& fileName ) const{

	OutputFiles::m_logFile << "# Write nodes" << std::endl;
	m_nodeList.writeNodeListToIntermediateFileStep3( fileName );

}

// Read triangles without the ones belonging to supert triangle from intermediate file
void TriangleList::readTrianglesFromIntermediateFile( const std::string& fileName ){

	OutputFiles::m_logFile << "# Read triangles" << std::endl;

	std::ifstream ifs( fileName.c_str() );
	if( !ifs.is_open() ) {
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << " !!" << std::endl;
		exit(1);
	}

	int numTriangles(0);
	ifs >> numTriangles;

	for( int iTriangle = 0; iTriangle < numTriangles; ++iTriangle ){

		int nodeIDs[3] = { 0, 0, 0 };
		int triangleIDAdj[3] = { 0, 0, 0 };
		int edgeIDAdj[3] = { 0, 0, 0 };
		int domainType = CommonParameters::UNKNOWN;

		int ibuf(0);
		ifs >> ibuf;

		if( ibuf != iTriangle ){
			OutputFiles::m_logFile << " Error : Triangle IDs must be sequence number from zero" << std::endl;
			exit(1);
		}

		for( int i = 0; i < 3; ++i ){
			ifs >> nodeIDs[i];
		}

		for( int i = 0; i < 3; ++i ){
			ifs >> triangleIDAdj[i];
		}

		for( int i = 0; i < 3; ++i ){
			ifs >> edgeIDAdj[i];
		}

		ifs >> domainType;

		m_triangles.push_back( Triangle( nodeIDs, triangleIDAdj, edgeIDAdj, false, domainType ) );

	}

	ifs.close();

}

// Read nodes of triangles without the ones of supert triangle from intermediate file
void TriangleList::readNode2DListFromIntermediateFile( const std::string& fileName ){

	OutputFiles::m_logFile << "# Read nodes" << std::endl;

	m_nodeList.readNode2DListFromIntermediateFile( fileName );

}

// Read nodes of triangles without the ones of supert triangle from intermediate file after step 3
void TriangleList::readNodeListFromIntermediateFileStep3( const std::string& fileName ){

	OutputFiles::m_logFile << "# Read nodes" << std::endl;

	m_nodeList.readNodeListFromIntermediateFileStep3( fileName );

}

// Relate nodes to triangles
void TriangleList::relateNodesToTriangles(){

	OutputFiles::m_logFile << "# Relate nodes to triangles" << std::endl;

	int icount(0);
	for( std::vector<Triangle>::iterator itr = m_triangles.begin(); itr != m_triangles.end(); ++itr, ++icount ){

		for( int iNode = 0; iNode < 3; ++iNode ){
			m_nodeList.getPointerToNode( itr->getNodeID(iNode) )->addTriangleIDs(icount);
		}
		
	}

#ifdef _DEBUG_WRITE
	writeTrianglesToVTK("triangles_debug.vtk");
#endif

}

// Assign location flag to nodes
void TriangleList::assignLocationToNodes(){

	OutputFiles::m_logFile << "# Assign location flag to nodes" << std::endl;

	int icount(0);
	for( std::vector<Triangle>::iterator itr = m_triangles.begin(); itr != m_triangles.end(); ++itr, ++icount ){

		Node::Location locOfNode = Node::UNKNOWN;

		const int domainType = itr->getDomainType();

		if( domainType == CommonParameters::UNKNOWN ){
			OutputFiles::m_logFile << " Error : Domain type of the element " << icount << " is unknowon." << std::endl;
			exit(1);
		}
		else if( domainType == CommonParameters::SEA ){
			locOfNode = Node::SEA;
		}
		else if( domainType == CommonParameters::LAND ){
			locOfNode = Node::LAND;
		}
		else if( domainType == CommonParameters::LAKE ){
			locOfNode = Node::LAKE;
		}
		else if( domainType == CommonParameters::OUTSIDE_OF_DOMAIN ){
			locOfNode = Node::OUT_OF_DOMAIN;
		}

		for( int iNode = 0; iNode < 3; ++iNode ){
			if( m_nodeList.getPointerToNode( itr->getNodeID(iNode) )->getLocation() != Node::UNKNOWN ){
				continue;// Do not overwrite
			}
			m_nodeList.getPointerToNode( itr->getNodeID(iNode) )->setLocation( locOfNode );
		}

		for( int iAdj = 0; iAdj < 3; ++iAdj ){
			if( itr->getAdjacentTriangleID(iAdj) >= 0 &&  m_triangles[ itr->getAdjacentTriangleID(iAdj) ].getDomainType() != domainType ){
				if( domainType == CommonParameters::LAKE ||  m_triangles[ itr->getAdjacentTriangleID(iAdj) ].getDomainType() == CommonParameters::LAKE ){
					m_nodeList.getPointerToNode( itr->getNodeID( iAdj       ) )->setLocation( Node::LAKE_LINE );
					m_nodeList.getPointerToNode( itr->getNodeID( (iAdj+1)%3 ) )->setLocation( Node::LAKE_LINE );
				}else{
					m_nodeList.getPointerToNode( itr->getNodeID( iAdj       ) )->setLocation( Node::COAST_LINE );
					m_nodeList.getPointerToNode( itr->getNodeID( (iAdj+1)%3 ) )->setLocation( Node::COAST_LINE );
				}
			}
		}
		
	}

#ifdef _DEBUG_WRITE
	m_nodeList.writeNode2DListToVTK( "node_location_assigned.dbg.vtk" );
#endif

}

// Interpolate altitude and sea depth to nodes
void TriangleList::interpolateAltitudeToNodes(){
	
	OutputFiles::m_logFile << "# Interpolate altitude and sea depth to nodes" << std::endl;

	const double minSeaDepth = ( Control::getInstance() )->getMinSeaDepth();
	const double minAltitude = ( Control::getInstance() )->getMinAltitude();
	const double maxSeaDepth = ( Control::getInstance() )->getMaxSeaDepth();
	const double maxAltitude = ( Control::getInstance() )->getMaxAltitude();
	const double distanceUsedToAvoidTooSmallDenominator = ( Control::getInstance() )->getDistanceUsedToAvoidTooSmallDenominator();
	const LakeList* const ptrLakeList = LakeList::getInstance(); 
	const TopographyDataList* const ptrTopographyDataList = TopographyDataList::getInstance();

#ifdef _USE_OMP
	const int numThread = Control::getInstance()->getNumThreads();	
	omp_set_num_threads(numThread);
#endif

	const int numNodes = m_nodeList.getTotalNumberOfNode();

	int iNode = 0;
	double value= 0.0;
	std::map<int,int>::iterator itr;
#ifdef _USE_OMP
	#pragma omp parallel for default(shared) \
		private( iNode, value, itr )
#endif
	for( iNode = 0; iNode < numNodes; ++iNode ){
		value= 0.0;
		switch( m_nodeList.getPointerToNode( iNode )->getLocation() ){
			case Node::UNKNOWN:
				OutputFiles::m_logFile << " Error : Location of the node " << iNode << " is unknowon." << std::endl;
				exit(1);
				break;
			case Node::COAST_LINE:
				m_nodeList.getPointerToNode( iNode )->setCoordZ( 0.0 );
				break;
			case Node::LAKE_LINE:
				itr = m_nodeIDTriangle2BoundCurve.find( iNode );
				if( itr != m_nodeIDTriangle2BoundCurve.end() ){
					value = ptrLakeList->getLakeHeight( m_boundaryCurveList.getLakeIndexFromNodeID(itr->second) );
					m_nodeList.getPointerToNode( iNode )->setCoordZ( - value );
				}else{
					OutputFiles::m_logFile << " Error : Node ID " << iNode << " is not located on boundary curves." << std::endl;
					exit(1);
				}
				break;
			case Node::SEA:
				if( ptrTopographyDataList->interpolateZCoord( Node::SEA, m_nodeList.getPointerToNode( iNode )->getCoordXY(), distanceUsedToAvoidTooSmallDenominator, value ) ){
					// Topology data found near the point
					if( value < minSeaDepth ){
						value = minSeaDepth;
					}
					else if( value > maxSeaDepth ){
						value = maxSeaDepth;
					}
					m_nodeList.getPointerToNode( iNode )->setCoordZ( value );
				}
				else{
					// Topology data are NOT found near the point
					OutputFiles::m_logFile << " Warning : No topography data were found around node " << iNode << " ." << std::endl;
					OutputFiles::m_logFile << "           Thus, minimum depth " << minSeaDepth << " [km] give to the node ." << std::endl;
				}
				break;
			case Node::LAND:// Go through
			case Node::LAKE:
				if( ptrTopographyDataList->interpolateZCoord( Node::LAND, m_nodeList.getPointerToNode( iNode )->getCoordXY(), distanceUsedToAvoidTooSmallDenominator, value ) ){
					// Topology data found near the point
					if( value < minAltitude ){
						value = minAltitude;
					}
					else if( value > maxAltitude ){
						value = maxAltitude;
					}
					m_nodeList.getPointerToNode( iNode )->setCoordZ( - value );
				}
				else{
					// Topology data are NOT found near the point
					OutputFiles::m_logFile << " Warning : No topography data were found around node " << iNode << " ." << std::endl;
					OutputFiles::m_logFile << "           Thus, minimum altitude " << minAltitude << " [km] give to the node ." << std::endl;
				}
				break;
			default:
				OutputFiles::m_logFile << " Error : Location of the node " << iNode << " is improper." << std::endl;
				exit(1);
				break;
		}

	}

	editDeviatingHeight();

//#ifdef _DEBUG_WRITE
	writeTrianglesToVTK("triangles_with_height.vtk");
//#endif

}

// Write array convert node ID of boundary curve list to the one of triangle list 
void TriangleList::writeNodeIDBoundCurve2Triangle( const std::string& fileName ) const{

	std::ofstream ofs( fileName.c_str() ); 

	if( !ofs.is_open() ) {
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << " !!" << std::endl;
		exit(1);
	}

	ofs << std::setw(10) << static_cast<int>( m_nodeIDBoundCurve2Triangle.size() ) << std::endl;

	for( std::map<int,int>::const_iterator itr = m_nodeIDBoundCurve2Triangle.begin(); itr != m_nodeIDBoundCurve2Triangle.end(); ++itr ){
		ofs << std::setw(10) << itr->first << std::setw(10) << itr->second - 3 << std::endl;
	}
	
	ofs.close();

}

// Write array convert node ID of boundary curve list to the one of triangle list after step3
void TriangleList::writeNodeIDBoundCurve2TriangleStep3( const std::string& fileName ) const{

	std::ofstream ofs( fileName.c_str() ); 

	if( !ofs.is_open() ) {
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << " !!" << std::endl;
		exit(1);
	}

	ofs << std::setw(10) << static_cast<int>( m_nodeIDBoundCurve2Triangle.size() ) << std::endl;

	for( std::map<int,int>::const_iterator itr = m_nodeIDBoundCurve2Triangle.begin(); itr != m_nodeIDBoundCurve2Triangle.end(); ++itr ){
		ofs << std::setw(10) << itr->first << std::setw(10) << itr->second << std::endl;
	}
	
	ofs.close();

}

// Read array convert node ID of boundary curve list to the one of triangle list 
void TriangleList::readNodeIDBoundCurve2Triangle( const std::string& fileName ){

	std::ifstream ifs( fileName.c_str() ); 

	if( !ifs.is_open() ) {
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << " !!" << std::endl;
		exit(1);
	}

	int numPairs(-1);
	ifs >> numPairs;

	for( int i = 0; i < numPairs; ++i ){
		std::pair<int,int> nodePair( -1, -1 );
		ifs >> nodePair.first >> nodePair.second;
		m_nodeIDBoundCurve2Triangle.insert( nodePair );
		m_nodeIDTriangle2BoundCurve.insert( std::make_pair(nodePair.second, nodePair.first) );
	}
	
	ifs.close();

}

// Write PLC for TETGEN
void TriangleList::writePLCs(){

	OutputFiles::m_logFile << "# Make PLCs" << std::endl;

	std::vector<CommonParameters::XYZ> nodeList;
	nodeList.reserve(  m_nodeList.getTotalNumberOfNode() + 12 );

	std::vector<Facet> facetList;
	facetList.reserve( static_cast<int>( m_triangles.size() ) * 2 );

	makePLCs( nodeList, facetList );

	OutputFiles::m_logFile << "# Write PLCs" << std::endl;

	writePLCsToVTK( "plc.vtk", nodeList, facetList ); 

	writePLCsToPolyFile( "output.poly", nodeList, facetList ); 

}

// Get pointer to the class of boundary curve list
BoundaryCurveList* TriangleList::getPointerToBoundaryCurveList(){
	return &m_boundaryCurveList;
}

// Set data of outer boundary curve
void TriangleList::setDataOfOneOuterBoundaryCurve( const int numNodes, const int* nodeIDs, const CommonParameters::XY* coord2D, const int geologicalType ){

#ifdef _DEBUG_WRITE
	for( int iNode = 0; iNode < numNodes; ++iNode ){
		std::cout << "nodeIDs : " << nodeIDs[iNode] << std::endl; 
		std::cout << "coord2D : " << coord2D[iNode].X << " " << coord2D[iNode].Y << std::endl; 
	}
#endif

	getPointerToBoundaryCurveList()->setNode2DList( numNodes, coord2D );

	getPointerToBoundaryCurveList()->setOneOuterBoundaryCurve( numNodes, nodeIDs, geologicalType );

}

// Set type of the plane on which 2D mesh is created
void TriangleList::setPlaneLocation( const CommonParameters::Boundary& planeLocation, const int iLayer ){
	m_locationOfPlane = planeLocation;
	if( planeLocation == CommonParameters::LAYER ){
		m_indexLayerInterface = iLayer;
	}
}

//
//// Add new triangle
//void TriangleList::addNewTriangle( const int* nodes, const int* adjTriangles, const int* edgeOfAdjTriangles, const int location ){
//
//	m_triangles.push_back( Triangle( nodes, adjTriangles, edgeOfAdjTriangles, Triangle::SUPPER_TRIANGLE ) );
//
//	const int triangleID = m_triangles.size() - 1;
//	for( int iNode = 0; iNode < 3; ++iNode ){
//		m_nodeList.getPointerToNode( nodes[iNode] )->addTriangleIDs( triangleID );
//	}
//
//}

// Make rough triangles from nodes of boundaries
void TriangleList::makeRoughTriangles(){

	OutputFiles::m_logFile << "# Make rough triangles from the node of boundaries" << std::endl;

	const BoundaryCurveList* const ptrBoundaryCurveList = &m_boundaryCurveList;
	const int numOuterBoundary = ptrBoundaryCurveList->getTotalNumberOuterBoundaries();
	//const NodeList* const ptrNode2DList = ptrBoundaryCurveList->getPointerToNodeList();

	for( int iBounOuter = 0; iBounOuter < numOuterBoundary; ++iBounOuter ){

		OutputFiles::m_logFile << "# Insert nodes of the outer boundary " << iBounOuter << std::endl;

		//const int offset = m_nodeList.getTotalNumberOfNode();

		const int numNodes = ptrBoundaryCurveList->getNumNodeOuterBoundary( iBounOuter );
		
#ifdef _DEBUG_WRITE
		std::cout << "numNodes = " << numNodes << std::endl;
#endif

		int locType = Triangle::INSIDE_TRIANGLE;
		//int iNodeFirst = 0;
		//for( int iNode = 0; iNode < numNodes; ++iNode ){
		//	const int nodeIDBoundCurveList = ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getNodeID(iNode);
		//	insertNewNodeAndFlip( ptrBoundaryCurveList->getPointCoordOuterBoundary( iBounOuter, 0 ), locType, nodeIDBoundCurveList );
		//	if( locType != Triangle::ON_NODE ){
		//		iNodeFirst = iNode;
		//		break;
		//	}
		//}
		const int nodeIDBoundCurveList = ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getNodeID(0);
		const int nodeIDTriangleList = insertNewNodeAndFlip( ptrBoundaryCurveList->getPointCoordOuterBoundary( iBounOuter, 0 ), locType, nodeIDBoundCurveList );
		if( nodeIDTriangleList < 0 ){
			OutputFiles::m_logFile << " Error : Node " << 0 << " of the boundary curve ( " << nodeIDBoundCurveList << " ) can't be inserted." << std::endl;
			exit(1);
		}

		const int iNodeFirst(nodeIDTriangleList);
		int iNodePre(nodeIDTriangleList);

#ifdef _DEBUG_WRITE
		std::cout << "locType = " << locType << std::endl;
		std::cout << "iNodeFirst : " << iNodeFirst << std::endl;
#endif

//#ifdef _DEBUG_WRITE
		//std::ostringstream oss;
		//oss << "temp_triangle.outer." << iBounOuter << "_" << 0 << ".vtk";
		//writeTrianglesToVTK( oss.str().c_str() );
//#endif

		//for( int iNode = iNodeFirst + 1; iNode < numNodes; ++iNode ){
		for( int iNode = 1; iNode < numNodes; ++iNode ){
#ifdef _DEBUG_WRITE
			std::cout << "iNode = " << iNode << std::endl;
			std::cout << "coord = " << ptrBoundaryCurveList->getPointCoordOuterBoundary( iBounOuter, iNode ).X << " " << ptrBoundaryCurveList->getPointCoordOuterBoundary( iBounOuter, iNode ).Y << std::endl;
#endif
			const int nodeIDBoundCurveList = ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getNodeID(iNode);
			const int nodeIDTriangleList = insertNewNodeAndFlip( ptrBoundaryCurveList->getPointCoordOuterBoundary( iBounOuter, iNode ), locType, nodeIDBoundCurveList );

#ifdef _DEBUG_WRITE
			std::cout << "nodeIDBoundCurveList = " << nodeIDBoundCurveList << std::endl;
			std::cout << "nodeIDTriangleList = " << nodeIDTriangleList << std::endl;
#endif

			if( nodeIDTriangleList < 0 ){
				OutputFiles::m_logFile << " Error : Node " << iNode << " of the boundary curve ( " << nodeIDBoundCurveList << " ) can't be inserted." << std::endl;
				exit(1);
			}

#ifdef _DEBUG_WRITE
			std::cout << "locType = " << locType << std::endl;
#endif

			if( locType == Triangle::ON_NODE ){
				iNodePre = nodeIDTriangleList;
				continue;
			}

#ifdef _DEBUG_WRITE
			std::ostringstream oss;
			oss << "temp_triangle.outer." << iBounOuter << "_" << iNode << ".vtk";
			writeTrianglesToVTK( oss.str().c_str() );
#endif

			std::vector<int> triangleIDs;
			//searchElemBetweenNodes( iNodePre + offset, iNodePre + 1 + offset, triangleIDs );
			//iNodePre = iNode;
			searchElemBetweenNodes( iNodePre, nodeIDTriangleList, triangleIDs );
			if( !triangleIDs.empty() ){
				// Subdivide polygon to keep boundary
				subdividePolygon( triangleIDs, iNodePre, nodeIDTriangleList );
			}
			iNodePre = nodeIDTriangleList;
		}

#ifdef _DEBUG_WRITE
		std::cout << "iNodePre iNodeFirst : " << iNodePre << " " << iNodeFirst << std::endl;
#endif

		std::vector<int> triangleIDs;
		searchElemBetweenNodes( iNodePre, iNodeFirst, triangleIDs );
		if( !triangleIDs.empty() ){
			// Subdivide polygon to keep boundary
			subdividePolygon( triangleIDs, iNodePre, iNodeFirst );
		}

		//=================================================================================================================
		// Inner boundaries inclued
		//=================================================================================================================
		const int numBounInnerIncluded = ptrBoundaryCurveList->getNumInnerBoundaryIncluded( iBounOuter );
		if( numBounInnerIncluded > 0 ){
			OutputFiles::m_logFile << "# Insert nodes of the inner boundaries within the outer boundary " << iBounOuter << std::endl;
			OutputFiles::m_logFile << "# Total number of the inner boundaries within the outer boundary : " << numBounInnerIncluded << std::endl;
		}

		for( int iBounInnerIncluded = 0; iBounInnerIncluded < numBounInnerIncluded; ++iBounInnerIncluded ){
			const int iBounInner = ptrBoundaryCurveList->getInnerBoundaryIDIncluded(iBounOuter,iBounInnerIncluded);

			OutputFiles::m_logFile << "# Insert nodes of the inner boundary " << iBounInner << std::endl;

			//const int offset = m_nodeList.getTotalNumberOfNode();

			const int numNodes = ptrBoundaryCurveList->getNumNodeInnerBoundary( iBounInner );
		
			int locType = Triangle::INSIDE_TRIANGLE;
			//int iNodeFirst = 0;
			//for( int iNode = 0; iNode < numNodes; ++iNode ){
			//	const int nodeIDBoundCurveList = ( ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner ) )->getNodeID(iNode);
			//	insertNewNodeAndFlip( ptrBoundaryCurveList->getPointCoordInnerBoundary( iBounInner, 0 ), locType, nodeIDBoundCurveList );
			//	if( locType != Triangle::ON_NODE ){
			//		iNodeFirst = iNode;
			//		break;
			//	}
			//}
#ifdef _MOD_FOR_NMT
			const BoundaryCurveInner* ptrInnerBoundary = ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner );
			int nodeIDBoundCurveList = ptrInnerBoundary->getNodeID(0);
			int nodeIDTriangleList = -1;
			if( ptrInnerBoundary->getGeologicalType() == CommonParameters::NMT_DIPOLE ){
				std::map<int,int>::const_iterator itr = m_nodeIDBoundCurve2Triangle.find( nodeIDBoundCurveList );
				if( itr != m_nodeIDBoundCurve2Triangle.end() ){
					// Point of the NMT dipole is already inserted
					nodeIDTriangleList = itr->second;
					if( nodeIDTriangleList < 0 ){
						OutputFiles::m_logFile << " Error : Node ID (" << nodeIDTriangleList << ") is wrong." << std::endl;
						exit(1);
					}
				}
			}
			if( nodeIDTriangleList < 0 ){
				nodeIDTriangleList = insertNewNodeAndFlip( ptrBoundaryCurveList->getPointCoordInnerBoundary( iBounInner, 0 ), locType, nodeIDBoundCurveList );
				if( nodeIDTriangleList < 0 ){
					OutputFiles::m_logFile << " Error : Node " << 0 << " of the boundary curve ( " << nodeIDBoundCurveList << " ) can't be inserted." << std::endl;
					exit(1);
				}
			}
#else
			const int nodeIDBoundCurveList = ( ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner ) )->getNodeID(0);
			const int nodeIDTriangleList = insertNewNodeAndFlip( ptrBoundaryCurveList->getPointCoordInnerBoundary( iBounInner, 0 ), locType, nodeIDBoundCurveList );
			if( nodeIDTriangleList < 0 ){
				OutputFiles::m_logFile << " Error : Node " << 0 << " of the boundary curve ( " << nodeIDBoundCurveList << " ) can't be inserted." << std::endl;
				exit(1);
			}
#endif
			const int iNodeFirst(nodeIDTriangleList);
			int iNodePre(nodeIDTriangleList);

//#ifdef _DEBUG_WRITE
//			std::ostringstream oss;
//			oss << "temp_triangle.inner." << iBounInner << "_" << 0 << ".vtk";
//			writeTrianglesToVTK( oss.str().c_str() );
//#endif

			for( int iNode = 1; iNode < numNodes; ++iNode ){

#ifdef _DEBUG_WRITE
			std::cout << "* iNode = " << iNode << std::endl;
			std::cout << "* coord = " << ptrBoundaryCurveList->getPointCoordInnerBoundary( iBounInner, iNode ).X << " " << ptrBoundaryCurveList->getPointCoordInnerBoundary( iBounInner, iNode ).Y << std::endl;
#endif

				const int nodeIDBoundCurveList = ( ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner ) )->getNodeID(iNode);
				const int nodeIDTriangleList = insertNewNodeAndFlip( ptrBoundaryCurveList->getPointCoordInnerBoundary( iBounInner, iNode ), locType, nodeIDBoundCurveList );
				if( nodeIDTriangleList < 0 ){
					OutputFiles::m_logFile << " Error : Node " << iNode << " of the boundary curve ( " << nodeIDBoundCurveList << " ) can't be inserted." << std::endl;
					exit(1);
				}

#ifdef _DEBUG_WRITE
				std::cout << "* locType = " << locType << std::endl;
#endif

				if( locType == Triangle::ON_NODE ){
					iNodePre = nodeIDTriangleList;
					continue;
				}

#ifdef _DEBUG_WRITE
				std::ostringstream oss;
				oss << "temp_triangle.inner." << iBounInner << "_" << iNode << ".vtk";
				writeTrianglesToVTK( oss.str().c_str() );
#endif

				std::vector<int> triangleIDs;
				//searchElemBetweenNodes( iNodePre + offset, iNodePre + 1 + offset, triangleIDs );
				//iNodePre = iNode;
				searchElemBetweenNodes( iNodePre, nodeIDTriangleList, triangleIDs );
				if( !triangleIDs.empty() ){
					// Subdivide polygon to keep boundary
					subdividePolygon( triangleIDs, iNodePre, nodeIDTriangleList );
				}
				iNodePre = nodeIDTriangleList;

			}

#ifdef _MOD_FOR_NMT
			if( ( ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner ) )->getGeologicalType() != CommonParameters::NMT_DIPOLE ){
				// Exclude NMT dipole because line for NMT dipole is not closed
				std::vector<int> triangleIDs;
				searchElemBetweenNodes( iNodePre, iNodeFirst, triangleIDs );
				if( !triangleIDs.empty() ){
					// Subdivide polygon to keep boundary
					subdividePolygon( triangleIDs, iNodePre, iNodeFirst );
				}
			}
#else
			std::vector<int> triangleIDs;
			//searchElemBetweenNodes( iNodePre + offset, iNodeFirst + offset, triangleIDs );
			searchElemBetweenNodes( iNodePre, iNodeFirst, triangleIDs );
			if( !triangleIDs.empty() ){
				// Subdivide polygon to keep boundary
				subdividePolygon( triangleIDs, iNodePre, iNodeFirst );
			}
#endif

			//-----------------------------------------------------------------------------------------------------------------
			// Sub-inner boundaries inclued
			//-----------------------------------------------------------------------------------------------------------------
			const int numBounSubInnerIncluded = ptrBoundaryCurveList->getNumSubInnerBoundaryIncluded( iBounInner );

			if( numBounSubInnerIncluded > 0 ){
				OutputFiles::m_logFile << "# Insert nodes of the sub-inner boundaries within the inner boundary " << iBounInner << std::endl;
				OutputFiles::m_logFile << "# Total number of the sub-inner boundaries within the inner boundary : " << numBounSubInnerIncluded << std::endl;
			}

			for( int iBounSubInnerIncluded = 0; iBounSubInnerIncluded < numBounSubInnerIncluded; ++iBounSubInnerIncluded ){
				const int iBounSubInner = ptrBoundaryCurveList->getSubInnerBoundaryIDIncluded(iBounInner,iBounSubInnerIncluded);

				OutputFiles::m_logFile << "# Insert nodes of the sub-inner boundary " << iBounSubInner << std::endl;

				//const int offset = m_nodeList.getTotalNumberOfNode();

				const int numNodes = ptrBoundaryCurveList->getNumNodeSubInnerBoundary( iBounSubInner );
		
				int locType = Triangle::INSIDE_TRIANGLE;
				//int iNodeFirst = 0;
				//for( int iNode = 0; iNode < numNodes; ++iNode ){
				//	const int nodeIDBoundCurveList = ( ptrBoundaryCurveList->getPointerToSubInnerBoundary( iBounSubInner ) )->getNodeID(iNode);
				//	insertNewNodeAndFlip( ptrBoundaryCurveList->getPointCoordSubInnerBoundary( iBounSubInner, 0 ), locType, nodeIDBoundCurveList );
				//	if( locType != Triangle::ON_NODE ){
				//		iNodeFirst = iNode;
				//		break;
				//	}
				//}
				const int nodeIDBoundCurveList = ( ptrBoundaryCurveList->getPointerToSubInnerBoundary( iBounSubInner ) )->getNodeID(0);
				const int nodeIDTriangleList = insertNewNodeAndFlip( ptrBoundaryCurveList->getPointCoordSubInnerBoundary( iBounSubInner, 0 ), locType, nodeIDBoundCurveList );
				if( nodeIDTriangleList < 0 ){
					OutputFiles::m_logFile << " Error : Node " << 0 << " of the boundary curve ( " << nodeIDBoundCurveList << " ) can't be inserted." << std::endl;
					exit(1);
				}

				const int iNodeFirst(nodeIDTriangleList);
				int iNodePre(nodeIDTriangleList);

//#ifdef _DEBUG_WRITE
//				std::ostringstream oss;
//				oss << "temp_triangle.subinner." << iBounInner << "_" << 0 << ".vtk";
//				writeTrianglesToVTK( oss.str().c_str() );
//#endif

				for( int iNode = 1; iNode < numNodes; ++iNode ){

#ifdef _DEBUG_WRITE
				std::cout << "** iNode = " << iNode << std::endl;
				std::cout << "** coord = " << ptrBoundaryCurveList->getPointCoordSubInnerBoundary( iBounSubInner, iNode ).X << " " << ptrBoundaryCurveList->getPointCoordSubInnerBoundary( iBounSubInner, iNode ).Y << std::endl;
#endif

					const int nodeIDBoundCurveList = ( ptrBoundaryCurveList->getPointerToSubInnerBoundary( iBounSubInner ) )->getNodeID(iNode);
					const int nodeIDTriangleList = insertNewNodeAndFlip( ptrBoundaryCurveList->getPointCoordSubInnerBoundary( iBounSubInner, iNode ), locType, nodeIDBoundCurveList );
					if( nodeIDTriangleList < 0 ){
						OutputFiles::m_logFile << " Error : Node " << iNode << " of the boundary curve ( " << nodeIDBoundCurveList << " ) can't be inserted." << std::endl;
						exit(1);
					}

#ifdef _DEBUG_WRITE
					std::cout << "** locType = " << locType << std::endl;
#endif

					if( locType == Triangle::ON_NODE ){
						iNodePre = nodeIDTriangleList;
						continue;
					}

//#ifdef _DEBUG_WRITE
//					std::ostringstream oss;
//					oss << "temp_triangle.subinner." << iBounInner << "_" << iNode << ".vtk";
//					writeTrianglesToVTK( oss.str().c_str() );
//#endif

					std::vector<int> triangleIDs;
					searchElemBetweenNodes( iNodePre, nodeIDTriangleList, triangleIDs );
					if( !triangleIDs.empty() ){
						// Subdivide polygon to keep boundary
						subdividePolygon( triangleIDs, iNodePre, nodeIDTriangleList );
					}
					iNodePre = nodeIDTriangleList;
				}

				std::vector<int> triangleIDs;
				searchElemBetweenNodes( iNodePre, iNodeFirst, triangleIDs );
				if( !triangleIDs.empty() ){
					// Subdivide polygon to keep boundary
					subdividePolygon( triangleIDs, iNodePre, iNodeFirst );
				}

			}
			//-----------------------------------------------------------------------------------------------------------------

		}
		//=================================================================================================================

	}


}

//----- Do not delete for future use >>>>>
//// Remove out-of-domain triangle 
//void TriangleList::removeOutOfDomainTriangles(){
//
//	OutputFiles::m_logFile << "# Remove out-of-domain triangle." << std::endl;
//
//	int iTriangle(0);
//	for( std::vector<Triangle>::iterator itr = m_triangles.begin(); itr != m_triangles.end(); ++itr, ++iTriangle ){
//
//		std::vector<int> stack;
//		for( int iNode = 0; iNode < 3; ++iNode ){
//			if( 0 <= itr->getNodeID(iNode) || itr->getNodeID(iNode) <= 2 ){
//				stack.push_back(iNode);
//			}
//		}
//
//		const int numCount = static_cast<int>( stack.size() );
//		if( numCount == 1 ){
//
//			const int edgeID = ( stack[0] + 1 ) % 3;
//
//			const int triangleIDAdj = itr->getAdjacentTriangleID(edgeID);
//
//			if( triangleIDAdj < 0 ){
//				OutputFiles::m_logFile << " Error : Triangle face to a node of super triangle is -1." << std::endl;
//				exit(1);
//			}
//
//			const int edgeIDAdj = itr->getEdgeIDOfAdjacentTriangle(edgeID);
//			m_triangles[triangleIDAdj].setAdjacentTriangleID( edgeIDAdj, -1 );
//			m_triangles[triangleIDAdj].setEdgeIDOfAdjacentTriangle( edgeIDAdj, -1 );
//
//			for( int iNode = 0; iNode < 3; ++iNode ){
//				if( iNode != stack[0] ){
//					m_nodeList.getPointerToNode( itr->getNodeID(iNode) )->eraseTriangleIDs( iTriangle );
//				}
//			}
//
//			itr = m_triangles.erase( itr );
//
//		}
//		else if( numCount == 2 ){
//
//			for( int iNode = 0; iNode < 3; ++iNode ){
//				if( iNode != stack[0] && iNode != stack[1] ){
//					m_nodeList.getPointerToNode( itr->getNodeID(iNode) )->eraseTriangleIDs( iTriangle );
//				}
//			}
//			itr = m_triangles.erase( itr );
//
//		}
//		else if( numCount == 3 ){
//			OutputFiles::m_logFile << " Error : Triangle " << iTriangle << " share all nodes with the super triangle." << std::endl;
//			exit(1);
//		}
//
//	}
//
//	m_nodeList.removeNode(0);
//	m_nodeList.removeNode(1);
//	m_nodeList.removeNode(2);
//
//};
//----- Do not delete for future use <<<<<

// Assign domain type to triangles
void TriangleList::assignDomainTypeToTriangles(){

	OutputFiles::m_logFile << "# Assign domain type to triangles" << std::endl;

	const BoundaryCurveList* const ptrBoundaryCurveList = &m_boundaryCurveList;
	const int numOuterBoundary = ptrBoundaryCurveList->getTotalNumberOuterBoundaries();
	//const NodeList* const ptrNode2DList = ptrBoundaryCurveList->getPointerToNodeList();

	for( std::vector<Triangle>::iterator itr = m_triangles.begin(); itr != m_triangles.end(); ++itr ){
		itr->setDomainType( CommonParameters::OUTSIDE_OF_DOMAIN );
	}

	for( int iBounOuter = 0; iBounOuter < numOuterBoundary; ++iBounOuter ){

#ifdef _DEBUG_WRITE
		std::cout << "iBounOuter = " << iBounOuter << std::endl;
#endif

		assignDomainTypeToTrianglesWithinBoundaryCurve( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ), true );

		const int numBounInnerIncluded = ptrBoundaryCurveList->getNumInnerBoundaryIncluded( iBounOuter );
		for( int iBounInnerIncluded = 0; iBounInnerIncluded < numBounInnerIncluded; ++iBounInnerIncluded ){
			const int iBounInner = ptrBoundaryCurveList->getInnerBoundaryIDIncluded(iBounOuter,iBounInnerIncluded);

#ifdef _DEBUG_WRITE
			std::cout << "iBounInner = " << iBounInner << std::endl;
#endif
			assignDomainTypeToTrianglesWithinBoundaryCurve( ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner ), false );

			const int numBounSubInnerIncluded = ptrBoundaryCurveList->getNumSubInnerBoundaryIncluded( iBounInner );
			for( int iBounSubInnerIncluded = 0; iBounSubInnerIncluded < numBounSubInnerIncluded; ++iBounSubInnerIncluded ){
				const int iBounSubInner = ptrBoundaryCurveList->getSubInnerBoundaryIDIncluded(iBounInner,iBounSubInnerIncluded);

#ifdef _DEBUG_WRITE
				std::cout << "iBounSubInner = " << iBounSubInner << std::endl;
#endif
				assignDomainTypeToTrianglesWithinBoundaryCurve( ptrBoundaryCurveList->getPointerToSubInnerBoundary( iBounSubInner ), false );

			}

		}

	}

}

// Assign domain type to triangles within the specified boundary curve
void TriangleList::assignDomainTypeToTrianglesWithinBoundaryCurve( const BoundaryCurve* const ptrBoundaryCurve, const bool clockwise ){

#ifdef _MOD_FOR_NMT
	if( ptrBoundaryCurve->getGeologicalType() == CommonParameters::NMT_DIPOLE ){
		// Exclude NMT dipole because NMT dipole is not closed
		return;
	}
#endif

	std::vector<int> triangles;
	searchTrianglesWithinBoundaryCurve( ptrBoundaryCurve, clockwise, triangles );

	const int iDomainType = ptrBoundaryCurve->getGeologicalType();
	for( std::vector<int>::iterator itr = triangles.begin(); itr != triangles.end(); ++itr ){
		m_triangles[*itr].setDomainType( iDomainType );
	}

}

// Search triangles within the specified boundary curve
void TriangleList::searchTrianglesWithinBoundaryCurve( const BoundaryCurve* const ptrBoundaryCurve, const bool clockwise, std::vector<int>& triangles ){

	const int numNodeBoundCurve = ptrBoundaryCurve->getNumOfPoints();
	
	std::map<int,int>::const_iterator itr = m_nodeIDBoundCurve2Triangle.find( ptrBoundaryCurve->getNodeID(numNodeBoundCurve-1) );
	if( itr == m_nodeIDBoundCurve2Triangle.end() ){
		OutputFiles::m_logFile << " Error : A node of the boundary curve ( " << ptrBoundaryCurve->getNodeID(numNodeBoundCurve-1)<< " ) cannot be found in the node list of the surface mesh" << std::endl;
		exit(1);
	}
	int nodeIDPre = itr->second;
#ifdef _MOD_FOR_NMT
	int iNodeFirst = 0;

	if( ptrBoundaryCurve->getGeologicalType() == CommonParameters::NMT_DIPOLE ){
		// Because NMT dipole is not closed
		iNodeFirst = 1;
		itr = m_nodeIDBoundCurve2Triangle.find( ptrBoundaryCurve->getNodeID(0) );
		nodeIDPre = itr->second;
	}
#endif

	std::set<int> adjTriangleToBeSearched;
	std::vector<int> stackTriangleOuter;
	std::vector<int> stackTriangleInner;
#ifdef _MOD_FOR_NMT
	for( int iNode = iNodeFirst; iNode < numNodeBoundCurve; ++iNode ){
#else
	for( int iNode = 0; iNode < numNodeBoundCurve; ++iNode ){
#endif

#ifdef _DEBUG_WRITE
		std::cout << "iNode = " << iNode << std::endl;
#endif

		itr = m_nodeIDBoundCurve2Triangle.find( ptrBoundaryCurve->getNodeID(iNode) );
		if( itr == m_nodeIDBoundCurve2Triangle.end() ){
			OutputFiles::m_logFile << " Error : A node of the boundary curve cannot be found in the node list of the surface mesh" << std::endl;
			exit(1);
		}
		int nodeIDCur = itr->second;

#ifdef _DEBUG_WRITE
		std::cout << "nodeIDPre nodeIDCur : " << nodeIDPre << " " << nodeIDCur << std::endl;
#endif

		int triLeft(-1);
		int triRight(-1);
		//int edgeLeft(-1);
		//int edgeRight(-1);
		//getTriangleShareTwoNodes( nodeIDPre, nodeIDCur, triLeft, triRight, edgeLeft, edgeRight );
		getTriangleShareTwoNodes( nodeIDPre, nodeIDCur, triLeft, triRight );
		const int triInner = clockwise ? triRight : triLeft;
		const int triOuter = clockwise ? triLeft  : triRight;
		//const int edgeInner = clockwise ? edgeRight : edgeLeft;
		//const int edgeOuter = clockwise ? edgeLeft  : edgeRight;
		if( triInner < 0 ){
			OutputFiles::m_logFile << " Error : ID of an inner triangle is nagative !!" << std::endl;
			exit(1);
		}
		stackTriangleInner.push_back( triInner );
		stackTriangleOuter.push_back( triOuter );
		triangles.push_back(triInner);

#ifdef _DEBUG_WRITE
		std::cout << "triInner triOuter : " << stackTriangleInner.back() << " " << stackTriangleOuter.back() << std::endl;
		//std::cout << "triInner triOuter : " << triInner << " " << triOuter << std::endl;
		//std::cout << "edgeInner edgeOuter : " << edgeInner << " " << edgeOuter << std::endl;
		//std::cout << "triInner  : " << triInner << std::endl;
		//std::cout << "edgeInner : " << edgeInner << std::endl;
#endif
		nodeIDPre = nodeIDCur;
	}

	for( std::vector<int>::iterator itrTriInner = stackTriangleInner.begin(); itrTriInner != stackTriangleInner.end(); ++itrTriInner ){
		for( int iAdj = 0; iAdj < 3; ++iAdj ){
			const int triangleIDAdj = m_triangles[*itrTriInner].getAdjacentTriangleID( iAdj );

			if( find( stackTriangleOuter.begin(), stackTriangleOuter.end(), triangleIDAdj ) != stackTriangleOuter.end() ){
				// Triangle locate out of the domain
				continue;
			}
			if( triangleIDAdj < 0 || find( triangles.begin(), triangles.end(), triangleIDAdj ) != triangles.end() ){
				continue;
			}

#ifdef _DEBUG_WRITE
			std::cout << "triangleIDAdj : " << triangleIDAdj << std::endl;
#endif

			triangles.push_back(triangleIDAdj);
			adjTriangleToBeSearched.insert( triangleIDAdj );
		}
	}

	const int maxCount = static_cast<int>( m_triangles.size() );

	int icount(0);
	while( !adjTriangleToBeSearched.empty() ){

		if( ++icount > maxCount ){
			OutputFiles::m_logFile << " Error : Reach maximum count number " << maxCount << std::endl;
			exit(1);
		}

		std::set<int> stack;

		for( std::set<int>::iterator itr = adjTriangleToBeSearched.begin(); itr != adjTriangleToBeSearched.end(); ++itr ){

#ifdef _DEBUG_WRITE
			std::cout << "*tri : " << *itr << std::endl;
#endif
			for( int iAdj = 0; iAdj < 3; ++iAdj ){

				const int triangleIDAdj = m_triangles[*itr].getAdjacentTriangleID( iAdj );
				if( find( stackTriangleOuter.begin(), stackTriangleOuter.end(), triangleIDAdj ) != stackTriangleOuter.end() ){
					// Triangle locate out of the domain
					continue;
				}
				if( triangleIDAdj < 0 || find( triangles.begin(), triangles.end(), triangleIDAdj ) != triangles.end() ){
					continue;
				}

#ifdef _DEBUG_WRITE
				std::cout << "triangleIDAdj : " << triangleIDAdj << std::endl;
#endif

				triangles.push_back(triangleIDAdj);
				stack.insert( triangleIDAdj );
			}
		}

		adjTriangleToBeSearched.swap( stack );
	}

	std::sort(triangles.begin(), triangles.end());
	triangles.erase( std::unique( triangles.begin(), triangles.end() ), triangles.end() );

}

//// Search boundary curve containing the segment containing specified two nodes as ends
//const BoundaryCurve* TriangleList::searchBoundaryCurve( const int nodeID0, const int nodeID1, int& nodeIDBoundCurve0, int& nodeIDBoundCurve1 ){
//
//	std::map<int,int>::iterator itr = m_nodeIDTriangle2BoundCurve.find( nodeID0 );
//	if( itr == m_nodeIDTriangle2BoundCurve.end() ){
//		return NULL;
//	}
//	nodeIDBoundCurve0 = itr->second;
//
//	itr = m_nodeIDTriangle2BoundCurve.find( nodeID1 );
//	if( itr == m_nodeIDTriangle2BoundCurve.end() ){
//		return NULL;
//	}
//	nodeIDBoundCurve1 = itr->second;
//
//	const BoundaryCurveList* const ptrBoundaryCurveList = BoundaryCurveList::getInstance();
//	const int numOuterBoundary = ptrBoundaryCurveList->getTotalNumberOuterBoundaries();
//	const NodeList* const ptrNode2DList = ptrBoundaryCurveList->getPointerToNodeList();
//
//	for( int iBounOuter = 0; iBounOuter < numOuterBoundary; ++iBounOuter ){
//
//		const BoundaryCurveOuter* ptrBoundaryCurveOuter = ptrBoundaryCurveList->getPointerToOuterBoundary(iBounOuter);
//		if( ptrBoundaryCurveOuter->nodePairOfBoundaryCurve( nodeIDBoundCurve0, nodeIDBoundCurve1 ) ){
//			return ptrBoundaryCurveOuter;
//		}
//		
//		const int numBounInnerIncluded = ptrBoundaryCurveList->getNumInnerBoundaryIncluded( iBounOuter );
//		for( int iBounInner = 0; iBounInner < numBounInnerIncluded; ++iBounInner ){
//
//			const BoundaryCurveInner* ptrBoundaryCurveInner = ptrBoundaryCurveList->getPointerToInnerBoundary(iBounInner);
//			if( ptrBoundaryCurveInner->nodePairOfBoundaryCurve( nodeIDBoundCurve0, nodeIDBoundCurve1 ) ){
//				return ptrBoundaryCurveInner;
//			}
//
//			const int numBounSubInnerIncluded = ptrBoundaryCurveList->getNumSubInnerBoundaryIncluded( iBounInner );
//			for( int iBounSubInner = 0; iBounSubInner < numBounSubInnerIncluded; ++iBounSubInner ){
//
//				const BoundaryCurveSubInner* ptrBoundaryCurveSubInner = ptrBoundaryCurveList->getPointerToSubInnerBoundary(iBounSubInner);
//				if( ptrBoundaryCurveSubInner->nodePairOfBoundaryCurve( nodeIDBoundCurve0, nodeIDBoundCurve1 ) ){
//					return ptrBoundaryCurveSubInner;
//				}
//
//			}
//
//		}
//
//	}
//
//}

// Get triangle sharing the specified two nodes
//void TriangleList::getTriangleShareTwoNodes( const int node0, const int node1, int& triLeft, int& triRight, int& edgeLeft, int& edgeRight ){
void TriangleList::getTriangleShareTwoNodes( const int node0, const int node1, int& triLeft, int& triRight ){
	
	int icount(0);

	std::vector<int> triangles = m_nodeList.getPointerToNode( node0 )->getTriangleIDs();

#ifdef _DEBUG_WRITE
	std::cout << "triangles" << std::endl;
	for( std::vector<int>::iterator itrTriangle = triangles.begin(); itrTriangle != triangles.end(); ++itrTriangle ){
		std::cout << *itrTriangle << std::endl;
	}
#endif

	for( std::vector<int>::iterator itrTriangle = triangles.begin(); itrTriangle != triangles.end(); ++itrTriangle ){

#ifdef _DEBUG_WRITE
		std::cout << "trianglesID : " << *itrTriangle << std::endl;
#endif

		int iNode0 = -1;
		int iNode1 = -1;
		int nodeOther = -1;
		for( int iNode = 0; iNode < 3; ++iNode ){
			const int nodeIDGlobal = m_triangles[*itrTriangle].getNodeID( iNode );
			if( node0 == nodeIDGlobal ){
				iNode0 = iNode;
			}
			else if( node1 == nodeIDGlobal ){
				iNode1 = iNode;
			}
			else{
				nodeOther = nodeIDGlobal;
			}
		}

#ifdef _DEBUG_WRITE
		std::cout << "iNode0 iNode1 nodeOther : " << iNode0 << " " << iNode1 << " " << nodeOther << std::endl;
#endif

		if( iNode0 >= 0 && iNode1 >= 0 ){

			int edgeID(0);
			switch( iNode0 + iNode1 ){
				case 1:
					edgeID = 0;
					break;
				case 2:
					edgeID = 2;
					break;
				case 3:
					edgeID = 1;
					break;
				default:
					OutputFiles::m_logFile << " Error : Wrong local node pair ( " << iNode0 << ", " << iNode1 << " ) ." << std::endl;
					exit(1);
					break;
			}

			if( Util::locateRightHandSide( m_nodeList.getCoordXYOfPoints( node0 ), m_nodeList.getCoordXYOfPoints( node1 ), m_nodeList.getCoordXYOfPoints( nodeOther ) ) ){
				//edgeRight = edgeID;
				triRight = *itrTriangle;
			}
			else{
				//edgeLeft = edgeID;
				triLeft = *itrTriangle;
			}

#ifdef _DEBUG_WRITE
			std::cout << " triLeft  triRight : " << triLeft << " " << triRight << std::endl;
			//std::cout << "edgeLeft edgeRight : " << edgeLeft << " " << edgeRight << std::endl;
#endif

			if( triRight >= 0 && triLeft >= 0 && triLeft != triRight ){
				return;
			}
			else if( m_triangles[*itrTriangle].getAdjacentTriangleID(edgeID) < 0 && ( triRight >= 0 || triLeft >= 0 ) ){
				return;
			}
		}

	}

	OutputFiles::m_logFile << " Error : Triangle pair sharing the specified two nodes ( " << node0 << ", " << node1 << " ) cannot be found." << std::endl;
	exit(1);

}

// Refine triangles
void TriangleList::refineTrianglesHavingLongEdges(){

	OutputFiles::m_logFile << "# Refine triangle having long edges" << std::endl;

	bool continueRefinement = true;

	std::vector<int> triangleToBeRefined;

	const double alpha = ( Control::getInstance() )->getAlpha();
	const double beta = ( Control::getInstance() )->getBeta();

	for( std::vector<Triangle>::iterator itr = m_triangles.begin(); itr != m_triangles.end(); ++itr ){
		itr->setSatisfyCriteria( false );// Initialize
	}

	int icount(0);
	int maxCount = ( Control::getInstance() )->getMaxIterNumRefiningTriangleWithLongEdge();
	while(continueRefinement){

		if( ++icount > maxCount ){
			//OutputFiles::m_logFile << " Error : Loop count reach maximum value in refining triangles with long edge " << maxCount << "." << std::endl;
			//exit(1);
			OutputFiles::m_logFile << "# Reach maximum interation number ( " << maxCount << " )." << std::endl;
			return;
		}

		//std::map<int,int> stackTriangle;
		std::vector<CommonParameters::XY> nodesInserted;

		int iTriangle(0);
		for( std::vector<Triangle>::iterator itr = m_triangles.begin(); itr != m_triangles.end(); ++itr, ++iTriangle ){

			if( itr->getSatisfyCriteria() ){
				continue;
			}

			if( m_triangles[iTriangle].getDomainType() == CommonParameters::OUTSIDE_OF_DOMAIN ||
				m_triangles[iTriangle].getDomainType() == CommonParameters::UNKNOWN ){
				continue;
			}

			const int nodes[3] = { itr->getNodeID(0), itr->getNodeID(1), itr->getNodeID(2) };
			const double disMax[3] = {
				calcMaximumEdgeLength(nodes[0]),
				calcMaximumEdgeLength(nodes[1]),
				calcMaximumEdgeLength(nodes[2]),
			};

			bool found(false);
			for( int iEdge = 0; iEdge < 3; ++iEdge ){
				const double edgeLength = Util::calcEdgeLength(m_nodeList.getCoordXYOfPoints( nodes[iEdge] ), m_nodeList.getCoordXYOfPoints( nodes[(iEdge+1)%3] ));
				if( edgeLength > alpha * std::min( disMax[iEdge], disMax[(iEdge+1)%3] ) ){
					found = true;
					break;
				}
			}
			if( !found ){
				itr->setSatisfyCriteria( true );
				continue;// Do not refine
			}

			bool tooClose(false);
			const CommonParameters::XY centerCoord = itr->getCoordGravCenter( &m_nodeList );
			const double maxLeng = m_locationOfPlane == CommonParameters::LAYER ? ( Control::getInstance() )->calcMaximumEdgeLength( centerCoord, m_locationOfPlane, m_indexLayerInterface ) : ( Control::getInstance() )->calcMaximumEdgeLength( centerCoord, m_locationOfPlane );
			for( std::vector<CommonParameters::XY>::iterator itr = nodesInserted.begin(); itr != nodesInserted.end(); ++itr ){
				if( Util::calcEdgeLength( centerCoord, *itr ) < beta * maxLeng ){// too close to another node to be inserted
					tooClose = true;
					break;// Do not refine
				}
			}
			for( int iNode = 0; iNode < 3; ++iNode ){
				if( Util::calcEdgeLength( centerCoord, m_nodeList.getCoordXYOfPoints( nodes[iNode] ) ) < beta * maxLeng ){// too close to another node
					tooClose = true;
					break;// Do not refine
				}
			}
			if( tooClose ){
				continue;// Do not refine
			}

			nodesInserted.push_back( centerCoord );

		}


		if( nodesInserted.empty() ){
			continueRefinement = false;
		}

		OutputFiles::m_logFile << "# Iteration : " << std::setw(10) << icount << ", Number of nodes to be inserted : " << std::setw(10) << nodesInserted.size() << std::endl;

		// Add new nodes
		for( std::vector<CommonParameters::XY>::iterator itr = nodesInserted.begin(); itr != nodesInserted.end(); ++itr ){

			int locType = Triangle::INSIDE_TRIANGLE;
			insertNewNodeAndFlip( *itr, locType, -1 );

		}

		//std::ostringstream oss;
		//oss << "surface_triangle.refine1_iter" << icount << ".vtk";
		//writeTrianglesToVTK( oss.str().c_str() );

	}
	
}

// Calculate maximum edge length of the specified node
double TriangleList::calcMaximumEdgeLength( const int iNode ){

	const CommonParameters::XY coordXY = m_nodeList.getPointerToNode( iNode )->getCoordXY();
	if( m_locationOfPlane == CommonParameters::LAYER ){
		return ( Control::getInstance() )->calcMaximumEdgeLength( coordXY, m_locationOfPlane, m_indexLayerInterface );
	}
	else{
		return ( Control::getInstance() )->calcMaximumEdgeLength( coordXY, m_locationOfPlane );
	}

}


// Insert node pair of element having large angle to container
bool TriangleList::insertNodePairToContainer( const int triangleID, ArrayType& container ){

	const double maxAngle = CommonParameters::DEG2RAD * ( Control::getInstance() )->getThresholdAngle();
	const double beta = ( Control::getInstance() )->getBeta();

	const int domainType = m_triangles[triangleID].getDomainType();
	if( domainType == CommonParameters::OUTSIDE_OF_DOMAIN ||
		domainType == CommonParameters::UNKNOWN ){
		return false;
	}

	const int nodes[3] = {
		m_triangles[triangleID].getNodeID(0),
		m_triangles[triangleID].getNodeID(1),
		m_triangles[triangleID].getNodeID(2)
	};

	const CommonParameters::XY coord[3] = {
		m_nodeList.getCoordXYOfPoints( nodes[0] ),
		m_nodeList.getCoordXYOfPoints( nodes[1] ),
		m_nodeList.getCoordXYOfPoints( nodes[2] )
	};

	double angle[3] = { 0.0, 0.0, 0.0 };

	//int iNodeAngleMin(0);
	//double angleMin( CommonParameters::PI );
	int iNodeAngleMax(0);
	double angleMax( 0.0 );
	for( int iNode = 0; iNode < 3; ++iNode ){
		angle[iNode] = fabs( Util::calcAngle( coord[(iNode+2)%3], coord[iNode], coord[(iNode+1)%3] ) );
		//if( angle[iNode] < angleMin ){
		//	angleMin = angle[iNode];
		//	iNodeAngleMin = iNode;
		//}
		if( angle[iNode] > angleMax ){
			angleMax = angle[iNode];
			iNodeAngleMax = iNode;
		}
	}

	if( angleMax <= maxAngle ){
		return false;
	}

	const int iNode0 = (iNodeAngleMax+1)%3;
	const int iNode1 = (iNodeAngleMax+2)%3;
	const CommonParameters::XY coordInsert = {
		0.5*( coord[iNode0].X + coord[iNode1].X ),
		0.5*( coord[iNode0].Y + coord[iNode1].Y )
	};

	const double maxLeng = ( Control::getInstance() )->calcMaximumEdgeLength( coordInsert, m_locationOfPlane );
	for( int iNode = 0; iNode < 3; ++iNode ){
#ifdef _DEBUG_WRITE
		std::cout << "maxLeng length : " << maxLeng << " " << Util::calcEdgeLength( coordInsert, coord[iNode] ) << std::endl;
#endif
		if( Util::calcEdgeLength( coordInsert, coord[iNode] ) < beta * maxLeng ){// too close to another node
#ifdef _DEBUG_WRITE
			std::cout << "Not add to container" << std::endl;
#endif
			return false;// Do not refine
		}
	}
	const int edgeID = ( iNodeAngleMax + 1 ) % 3;
	const int adjTriangleID = m_triangles[triangleID].getAdjacentTriangleID( edgeID );
	if( adjTriangleID >= 0 ){
		for( int iNode = 0; iNode < 3; ++iNode ){
			const int nodeIDAdj = m_triangles[adjTriangleID].getNodeID(iNode);
#ifdef _DEBUG_WRITE
			std::cout << "maxLeng length : " << maxLeng << " " << Util::calcEdgeLength( coordInsert, m_nodeList.getCoordXYOfPoints( nodeIDAdj ) ) << std::endl;
#endif
			if( Util::calcEdgeLength( coordInsert, m_nodeList.getCoordXYOfPoints( nodeIDAdj ) ) < beta * maxLeng ){// too close to another node
#ifdef _DEBUG_WRITE
				std::cout << "Not add to container" << std::endl;
#endif
				return false;// Do not refine
			}
		}
	}

#ifdef _DEBUG_WRITE
	std::cout << "Add to container" << std::endl;
#endif

	NodePair pairOfNode = std::make_pair( m_triangles[triangleID].getNodeID(iNode0), m_triangles[triangleID].getNodeID(iNode1) );
	container.insert( std::make_pair( angleMax, pairOfNode ) );

	return true;

}

// Refine triangles with large anle
void TriangleList::refineTrianglesHavingLargeAngle( const bool notInsertNodesOnBoundary ){

	OutputFiles::m_logFile << "# Refine triangle with large angle" << std::endl;

	//std::map<int,int> stackTriangle;
	//std::vector<CommonParameters::XY> nodesInserted;
	ArrayType toBeInserted;

	int iTriangle(0);
	for( std::vector<Triangle>::iterator itr = m_triangles.begin(); itr != m_triangles.end(); ++itr, ++iTriangle ){
		insertNodePairToContainer( iTriangle, toBeInserted );
	}
	
	int icount(0);
	int maxCount = ( Control::getInstance() )->getMaxIterNumRefiningTriangleWithLargeAngle();

	bool continueRefinement = true;
		
	while(continueRefinement){

		if( toBeInserted.empty() ){
			OutputFiles::m_logFile << "# Finish refinement" << std::endl;
			return;
		}

		if( ++icount > maxCount ){
			OutputFiles::m_logFile << "# Reach maximum interation number ( " << maxCount << " )." << std::endl;
			return;
		}

		ArrayType::iterator itrInserted = toBeInserted.end();
		advance( itrInserted, -1 );

		OutputFiles::m_logFile << "# Largest angle[deg] : "	<< std::setw(15) << itrInserted->first * CommonParameters::RAD2DEG << std::endl;

		const int nodeIDBoundCurve0 = getNodeIDBoundaryCurveList( (itrInserted->second).first );
		const int nodeIDBoundCurve1 = getNodeIDBoundaryCurveList( (itrInserted->second).second );
			
		BoundaryCurveList* const ptrBoundaryCurveList = &m_boundaryCurveList;

		int nodeIDInsertedBoundCurveList = -1;
		if( nodeIDBoundCurve0 >= 0 && nodeIDBoundCurve1 >= 0 ){
			if( notInsertNodesOnBoundary ){
				toBeInserted.erase( itrInserted );
				continue;
			}
#ifdef _DEBUG_WRITE
			std::cout << "nodeIDBoundCurve0 nodeIDBoundCurve1 : " << nodeIDBoundCurve0 << " " << nodeIDBoundCurve1 << std::endl;
			std::cout << "Insert node to boudary curve" << std::endl;
#endif
			nodeIDInsertedBoundCurveList = ptrBoundaryCurveList->addNewNodeBetweenSegment( nodeIDBoundCurve0, nodeIDBoundCurve1 );
		}

		const CommonParameters::XY coord[2] = {
			m_nodeList.getCoordXYOfPoints( (itrInserted->second).first ),
			m_nodeList.getCoordXYOfPoints( (itrInserted->second).second )
		};
		const CommonParameters::XY coordInsert = {
			0.5*( coord[0].X + coord[1].X ),
			0.5*( coord[0].Y + coord[1].Y )
		};

#ifdef _DEBUG_WRITE
		std::cout << "coord[0] : " << coord[0].X << " " << coord[0].Y << std::endl;
		std::cout << "coord[1] : " << coord[1].X << " " << coord[1].Y << std::endl;
		std::cout << "coordInsert : " << coordInsert.X << " " << coordInsert.Y << std::endl;
#endif

		int locType = Triangle::ON_EDGE;
		const int nodeIDInsertedTriangleList = insertNewNodeAndFlip( coordInsert, locType, nodeIDInsertedBoundCurveList );

#ifdef _DEBUG_WRITE
		std::cout << "nodeIDInsertedTriangleList : " << nodeIDInsertedTriangleList << std::endl;
		std::cout << "nodeIDInsertedBoundCurveList : " << nodeIDInsertedBoundCurveList << std::endl;
#endif

		//if( locType != Triangle::ON_EDGE ){
		//	OutputFiles::m_logFile << " Error : Specified point ( " << coordInsert.X << ", " << coordInsert.Y << " ) was not inserted on edge." << std::endl;
		//	exit(1);
		//}

		toBeInserted.erase( itrInserted );

//		if( nodeIDInsertedTriangleList < 0 ){
//			if( nodeIDInsertedBoundCurveList >= 0 && ( nodeIDInsertedTriangleList < 0 || locType != Triangle::ON_EDGE ) ){
//				ptrBoundaryCurveList->removeNode( nodeIDInsertedBoundCurveList );
//#ifdef _DEBUG_WRITE
//				std::cout << "Remove node from boudary curve" << std::endl;
//#endif
//			}
//			continue;
//		}
		if( nodeIDInsertedTriangleList < 0 || locType != Triangle::ON_EDGE ){
			if( nodeIDInsertedBoundCurveList >= 0 ){
				ptrBoundaryCurveList->removeNode( nodeIDInsertedBoundCurveList );
#ifdef _DEBUG_WRITE
				std::cout << "Remove node from boudary curve" << std::endl;
#endif
			}
			continue;
		}

#ifdef _DEBUG_WRITE
		if( icount % 100 == 0 ){
			std::ostringstream oss;
			oss << "temp_triangle.refined2_" << icount << ".vtk";
			writeTrianglesToVTK( oss.str().c_str() );
		}
#endif

		std::vector<int> triangles = m_nodeList.getTriangleIDs( nodeIDInsertedTriangleList );
		for( std::vector<int>::iterator itr = triangles.begin(); itr != triangles.end(); ++itr ){
#ifdef _DEBUG_WRITE
			std::cout << "*itr, : " << *itr << std::endl;
#endif
			insertNodePairToContainer( *itr, toBeInserted );
		}

	}
	
}

// Insert node to the triangles all of nodes of which locate on the coast line
void TriangleList::insertNodeToTriangleWithAllNodesOnCoast(){
	
	OutputFiles::m_logFile << "# Insert node to the triangles all nodes of which locate on the coast line" << std::endl;

	bool continueRefinement = true;

	int icount(0);
	const int maxCount = 1000;
	//while(true){
	while(continueRefinement){

		if( ++icount > maxCount ){
			OutputFiles::m_logFile << "# Reach maximum interation number ( " << maxCount << " )." << std::endl;
			return;
		}

		//--------------------------------------------------------------------------------------------------
		//----- Insert nodes on the center of the triangle all of which nodes locate on the coast line -----
		//--------------------------------------------------------------------------------------------------
		std::vector<CommonParameters::XY> nodesInsertedOnCenter;

		for( std::vector<Triangle>::iterator itr = m_triangles.begin(); itr != m_triangles.end(); ++itr ){

			if( itr->getDomainType() != CommonParameters::SEA && itr->getDomainType() != CommonParameters::LAKE ){
				continue; //Consider only the sea or lakes
			}

			//std::vector<int> nodeIDs;
			int icount(0);
			for( int iNode = 0; iNode < 3; ++iNode ){
				const int nodeID = itr->getNodeID(iNode);
				if( m_nodeList.getPointerToNode( nodeID )->getLocation() == Node::COAST_LINE || m_nodeList.getPointerToNode( nodeID )->getLocation() == Node::LAKE_LINE ){
					++icount;
				}
			}

			//if( static_cast<int>( nodeIDs.size() ) == 3 ){
			if( icount >= 3 ){

				//locTypes.push_back( Triangle::INSIDE_TRIANGLE );

				const CommonParameters::XY centerCoord = itr->getCoordGravCenter( &m_nodeList );
				nodesInsertedOnCenter.push_back( centerCoord );

				//const CommonParameters::XY centerCoord = itr->getCoordGravCenter( &m_nodeList );
				//int locType = Triangle::INSIDE_TRIANGLE;
				//insertNewNodeAndFlip( centerCoord, locType, -1 );
			}

		}

		OutputFiles::m_logFile << "# Iteration : " << std::setw(10) << icount << ", Number of nodes to be inserted on center : " << std::setw(10) << nodesInsertedOnCenter.size() << std::endl;

		// Add new nodes
		for( std::vector<CommonParameters::XY>::iterator itr = nodesInsertedOnCenter.begin(); itr != nodesInsertedOnCenter.end(); ++itr ){
			int locType = Triangle::INSIDE_TRIANGLE;
			insertNewNodeAndFlip( *itr, locType, -1 );
		}
		//--------------------------------------------------------------------------------------------------

		std::vector< std::pair<int,int> > nodePair;

		//--------------------------------------------------------------------------------------------------
		//----- Insert nodes on the edges of the triangle two of which nodes locate on the coast line ------
		//--------------------------------------------------------------------------------------------------
		for( std::vector<Triangle>::iterator itr = m_triangles.begin(); itr != m_triangles.end(); ++itr ){
			if( itr->getDomainType() != CommonParameters::SEA ){
				continue; //Consider only the sea
			}
			std::vector<int> nodeIDs;
			for( int iNode = 0; iNode < 3; ++iNode ){
				const int nodeID = itr->getNodeID(iNode);
				if( m_nodeList.getPointerToNode( nodeID )->getLocation() == Node::COAST_LINE ){
					nodeIDs.push_back( nodeID );
				}
			}
			if( static_cast<int>( nodeIDs.size() ) == 2 ){
				const std::map<int,int>::iterator itrTemp1 = m_nodeIDTriangle2BoundCurve.find( nodeIDs[0] );
				const std::map<int,int>::iterator itrTemp2 = m_nodeIDTriangle2BoundCurve.find( nodeIDs[1] );
				if( itrTemp1 != m_nodeIDTriangle2BoundCurve.end() && itrTemp2 != m_nodeIDTriangle2BoundCurve.end() &&
					m_boundaryCurveList.nodePairOfBoundaryCurve( itrTemp1->second, itrTemp2->second ) ){
					// Node pair is located on a boundary curve
					continue;
				}
				const int nd1 = nodeIDs[0] < nodeIDs[1] ? nodeIDs[0] : nodeIDs[1];
				const int nd2 = nodeIDs[0] < nodeIDs[1] ? nodeIDs[1] : nodeIDs[0];
				nodePair.push_back( std::make_pair( nd1, nd2 ) );
			}
		}

		//-----------------------------------------------------------------------------------------------------
		//----- Insert nodes on the edges of the triangle two of which nodes locate on the same lake line -----
		//-----------------------------------------------------------------------------------------------------
		for( std::vector<Triangle>::iterator itr = m_triangles.begin(); itr != m_triangles.end(); ++itr ){
			if( itr->getDomainType() != CommonParameters::LAKE ){
				continue; //Consider only the sea
			}
			std::vector<int> nodeIDs;
			for( int iNode = 0; iNode < 3; ++iNode ){
				const int nodeID = itr->getNodeID(iNode);
				if( m_nodeList.getPointerToNode( nodeID )->getLocation() == Node::LAKE_LINE ){
					nodeIDs.push_back( nodeID );
				}
			}
			if( static_cast<int>( nodeIDs.size() ) == 2 ){
				const std::map<int,int>::iterator itrTemp1 = m_nodeIDTriangle2BoundCurve.find( nodeIDs[0] );
				const std::map<int,int>::iterator itrTemp2 = m_nodeIDTriangle2BoundCurve.find( nodeIDs[1] );
				if( itrTemp1 != m_nodeIDTriangle2BoundCurve.end() && itrTemp2 != m_nodeIDTriangle2BoundCurve.end() &&
					m_boundaryCurveList.nodePairOfBoundaryCurve( itrTemp1->second, itrTemp2->second ) ){
					// Node pair is located on a boundary curve
					continue;
				}
				const int nd1 = nodeIDs[0] < nodeIDs[1] ? nodeIDs[0] : nodeIDs[1];
				const int nd2 = nodeIDs[0] < nodeIDs[1] ? nodeIDs[1] : nodeIDs[0];
				nodePair.push_back( std::make_pair( nd1, nd2 ) );
			}
		}

		//------------------------------------------------------------------------------------------------------ Not delete for future use >>>>>
		//const BoundaryCurveList* const ptrBoundaryCurveList = &m_boundaryCurveList;
		////------------------------------------------------------------------------------------------------------------------------------------
		////----- Search node pairs located on outer boundary curves except the nodes of the segments belonging to the same boundary curve -----
		////------------------------------------------------------------------------------------------------------------------------------------
		//const int numOuterBoundary = ptrBoundaryCurveList->getTotalNumberOuterBoundaries();
		//for( int iBounOuter1 = 0; iBounOuter1 < numOuterBoundary; ++iBounOuter1 ){// Loop of outer boundary 1 >>>>>>
		//	const int numNodes1 = ptrBoundaryCurveList->getNumNodeOuterBoundary( iBounOuter1 );
		//	if( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter1 )->getGeologicalType() == CommonParameters::SEA ||
		//		ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter1 )->getGeologicalType() == CommonParameters::LAKE ){
		//		for( int iNode1 = 0; iNode1 < numNodes1; ++iNode1 ){
		//			const int numInnerBoundaryOfItsOwn = ptrBoundaryCurveList->getNumInnerBoundaryIncluded(iBounOuter1);
		//			for( int iBounInnerItsOwn = 0; iBounInnerItsOwn < numInnerBoundaryOfItsOwn; ++iBounInnerItsOwn ){// Loop of inner boundary of its own >>>>>>

		//				const int iBounInner = ptrBoundaryCurveList->getInnerBoundaryIDIncluded( iBounOuter1, iBounInnerItsOwn );

		//				if( ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner )->getGeologicalType() != CommonParameters::LAND ){
		//					continue;// Consider only boundary curves surrounding land
		//				}

		//				const int numNodes2 = ptrBoundaryCurveList->getNumNodeInnerBoundary( iBounInner );
		//				for( int iNode2 = 0; iNode2 < numNodes2; ++iNode2 ){

		//					const int nodeIDBoundCurve1 = ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter1 ) )->getNodeID(iNode1);
		//					const int nodeID1 = convertNodeIDBoundCurve2Triangle(nodeIDBoundCurve1);
		//					const int nodeIDBoundCurve2 = ( ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner ) )->getNodeID(iNode2);
		//					const int nodeID2 = convertNodeIDBoundCurve2Triangle(nodeIDBoundCurve2);

		//					const int iTriangle = shareSameTriangleExceptDomainBoundary( nodeID1, nodeID2 );

		//					if( iTriangle >= 0 && ( m_triangles[iTriangle].getDomainType() == CommonParameters::SEA || m_triangles[iTriangle].getDomainType() == CommonParameters::LAKE ) ){
		//						const int nd1 = nodeID1 < nodeID2 ? nodeID1 : nodeID2;
		//						const int nd2 = nodeID1 < nodeID2 ? nodeID2 : nodeID1;
		//						nodePair.push_back( std::make_pair( nd1, nd2 ) );
		//					}
		//				}
		//			}//------------------------------------------------------------------------- // Loop of inner boundary of its own <<<<<
		//		}//------------------------------------------------------------------------- // Loop of outer boundary 1 <<<<<
		//	}
		//	else if( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter1 )->getGeologicalType() == CommonParameters::LAND ){
		//		for( int iNode1 = 0; iNode1 < numNodes1; ++iNode1 ){
		//			for( int iBounOuter2 = 0; iBounOuter2 < numOuterBoundary; ++iBounOuter2 ){// Loop of outer boundary 2 >>>>>>
		//				if( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter2 )->getGeologicalType() != CommonParameters::LAND ){
		//					continue;// Consider only boundary curves surrounding land
		//				}

		//				const int numNodes2 = ptrBoundaryCurveList->getNumNodeOuterBoundary( iBounOuter2 );
		//				for( int iNode2 = 0; iNode2 < numNodes2; ++iNode2 ){

		//					// Except the nodes of the segments belonging to the boundary curve for the nodes on the same boundary curve
		//					if( iBounOuter1 == iBounOuter2 ){
		//						if( iNode1 == iNode2 || iNode2 == (iNode1+1)%numNodes1 || iNode1 == (iNode2+1)%numNodes1 ){
		//							continue;
		//						}
		//					}

		//					const int nodeIDBoundCurve1 = ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter1 ) )->getNodeID(iNode1);
		//					const int nodeID1 = convertNodeIDBoundCurve2Triangle(nodeIDBoundCurve1);
		//					const int nodeIDBoundCurve2 = ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter2 ) )->getNodeID(iNode2);
		//					const int nodeID2 = convertNodeIDBoundCurve2Triangle(nodeIDBoundCurve2);

		//					const int iTriangle = shareSameTriangleExceptDomainBoundary( nodeID1, nodeID2 );

		//					if( iTriangle >= 0 && ( m_triangles[iTriangle].getDomainType() == CommonParameters::SEA || m_triangles[iTriangle].getDomainType() == CommonParameters::LAKE ) ){
		//						const int nd1 = nodeID1 < nodeID2 ? nodeID1 : nodeID2;
		//						const int nd2 = nodeID1 < nodeID2 ? nodeID2 : nodeID1;
		//						nodePair.push_back( std::make_pair( nd1, nd2 ) );
		//					}

		//				}
		//			}//------------------------------------------------------------------------- // Loop of outer boundary 2 <<<<<
		//		}
		//	}
		//}//------------------------------------------------------------------------- // Loop of outer boundary 1 <<<<<


		////------------------------------------------------------------------------------------------------------------------------------------
		////----- Search node pairs located on inner boundary curves except the nodes of the segments belonging to the same boundary curve -----
		////------------------------------------------------------------------------------------------------------------------------------------
		//const int numInnerBoundary = ptrBoundaryCurveList->getTotalNumberInnerBoundaries();
		//for( int iBounInner1 = 0; iBounInner1 < numInnerBoundary; ++iBounInner1 ){// Loop of inner boundary 1 >>>>>>
		//	const int numNodes1 = ptrBoundaryCurveList->getNumNodeInnerBoundary( iBounInner1 );
		//	if( ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner1 )->getGeologicalType() == CommonParameters::SEA ||
		//		ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner1 )->getGeologicalType() == CommonParameters::LAKE ){
		//		for( int iNode1 = 0; iNode1 < numNodes1; ++iNode1 ){
		//			const int numSubInnerBoundaryOfItsOwn = ptrBoundaryCurveList->getNumSubInnerBoundaryIncluded(iBounInner1);
		//			for( int iBounSubInnerOfItsOwn = 0; iBounSubInnerOfItsOwn < numSubInnerBoundaryOfItsOwn; ++iBounSubInnerOfItsOwn ){// Loop of sub-inner boundary of its own >>>>>>

		//				const int iBounSubInner = ptrBoundaryCurveList->getSubInnerBoundaryIDIncluded( iBounInner1, iBounSubInnerOfItsOwn );

		//				if( ptrBoundaryCurveList->getPointerToSubInnerBoundary( iBounSubInner )->getGeologicalType() != CommonParameters::LAND ){
		//					continue;// Consider only boundary curves surrounding land
		//				}

		//				const int numNodes2 = ptrBoundaryCurveList->getNumNodeSubInnerBoundary( iBounSubInner );
		//				for( int iNode2 = 0; iNode2 < numNodes2; ++iNode2 ){

		//					const int nodeIDBoundCurve1 = ( ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner1 ) )->getNodeID(iNode1);
		//					const int nodeID1 = convertNodeIDBoundCurve2Triangle(nodeIDBoundCurve1);
		//					const int nodeIDBoundCurve2 = ( ptrBoundaryCurveList->getPointerToSubInnerBoundary( iBounSubInner ) )->getNodeID(iNode2);
		//					const int nodeID2 = convertNodeIDBoundCurve2Triangle(nodeIDBoundCurve2);

		//					const int iTriangle = shareSameTriangleExceptDomainBoundary( nodeID1, nodeID2 );

		//					if( iTriangle >= 0 && ( m_triangles[iTriangle].getDomainType() == CommonParameters::SEA || m_triangles[iTriangle].getDomainType() == CommonParameters::LAKE ) ){
		//						const int nd1 = nodeID1 < nodeID2 ? nodeID1 : nodeID2;
		//						const int nd2 = nodeID1 < nodeID2 ? nodeID2 : nodeID1;
		//						nodePair.push_back( std::make_pair( nd1, nd2 ) );
		//					}

		//				}
		//			}//------------------------------------------------------------------------- // Loop of sub-inner boundary of its own <<<<<
		//		}
		//	}
		//	else if( ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner1 )->getGeologicalType() == CommonParameters::LAND ){
		//		for( int iNode1 = 0; iNode1 < numNodes1; ++iNode1 ){
		//			for( int iBounInner2 = 0; iBounInner2 < numInnerBoundary; ++iBounInner2 ){// Loop of inner boundary 2 >>>>>>
		//				if( ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner2 )->getGeologicalType() != CommonParameters::LAND ){
		//					continue;// Consider only boundary curves surrounding land
		//				}

		//				const int numNodes2 = ptrBoundaryCurveList->getNumNodeInnerBoundary( iBounInner2 );
		//				for( int iNode2 = 0; iNode2 < numNodes2; ++iNode2 ){

		//					// Except the nodes of the segments belonging to the boundary curve for the nodes on the same boundary curve
		//					if( iBounInner1 == iBounInner2 ){
		//						if( iNode1 == iNode2 || iNode2 == (iNode1+1)%numNodes1 || iNode1 == (iNode2+1)%numNodes1 ){
		//							continue;
		//						}
		//					}

		//					const int nodeIDBoundCurve1 = ( ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner1 ) )->getNodeID(iNode1);
		//					const int nodeID1 = convertNodeIDBoundCurve2Triangle(nodeIDBoundCurve1);
		//					const int nodeIDBoundCurve2 = ( ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner2 ) )->getNodeID(iNode2);
		//					const int nodeID2 = convertNodeIDBoundCurve2Triangle(nodeIDBoundCurve2);

		//					const int iTriangle = shareSameTriangleExceptDomainBoundary( nodeID1, nodeID2 );

		//					if( iTriangle >= 0 && ( m_triangles[iTriangle].getDomainType() == CommonParameters::SEA || m_triangles[iTriangle].getDomainType() == CommonParameters::LAKE ) ){
		//						const int nd1 = nodeID1 < nodeID2 ? nodeID1 : nodeID2;
		//						const int nd2 = nodeID1 < nodeID2 ? nodeID2 : nodeID1;
		//						nodePair.push_back( std::make_pair( nd1, nd2 ) );
		//					}

		//				}
		//			}//------------------------------------------------------------------------- // Loop of inner boundary 2 <<<<<
		//		}
		//	}
		//}//------------------------------------------------------------------------- // Loop of inner boundary 1 <<<<<

		////----------------------------------------------------------------------------------------------------------------------------------------
		////----- Search node pairs located on sub-inner boundary curves except the nodes of the segments belonging to the same boundary curve -----
		////----------------------------------------------------------------------------------------------------------------------------------------
		//const int numSubInnerBoundary = ptrBoundaryCurveList->getTotalNumberSubInnerBoundaries();
		//for( int iBounSubInner1 = 0; iBounSubInner1 < numSubInnerBoundary; ++iBounSubInner1 ){// Loop of sub-inner boundary 1 >>>>>>
		//	if( ptrBoundaryCurveList->getPointerToSubInnerBoundary( iBounSubInner1 )->getGeologicalType() != CommonParameters::LAND ){
		//		continue;// Consider only boundary curves surrounding land
		//	}

		//	const int numNodes1 = ptrBoundaryCurveList->getNumNodeSubInnerBoundary( iBounSubInner1 );
		//	for( int iNode1 = 0; iNode1 < numNodes1; ++iNode1 ){

		//		for( int iBounSubInner2 = 0; iBounSubInner2 < numSubInnerBoundary; ++iBounSubInner2 ){// Loop of sub-inner boundary 2 >>>>>>
		//			if( ptrBoundaryCurveList->getPointerToSubInnerBoundary( iBounSubInner2 )->getGeologicalType() != CommonParameters::LAND ){
		//				continue;// Consider only boundary curves surrounding land
		//			}

		//			const int numNodes2 = ptrBoundaryCurveList->getNumNodeSubInnerBoundary( iBounSubInner2 );
		//			for( int iNode2 = 0; iNode2 < numNodes2; ++iNode2 ){

		//				// Except the nodes of the segments belonging to the boundary curve for the nodes on the same boundary curve
		//				if( iBounSubInner1 == iBounSubInner2 ){
		//					if( iNode1 == iNode2 || iNode2 == (iNode1+1)%numNodes1 || iNode1 == (iNode2+1)%numNodes1 ){
		//						continue;
		//					}
		//				}

		//				const int nodeIDBoundCurve1 = ( ptrBoundaryCurveList->getPointerToSubInnerBoundary( iBounSubInner1 ) )->getNodeID(iNode1);
		//				const int nodeIDBoundCurve2 = ( ptrBoundaryCurveList->getPointerToSubInnerBoundary( iBounSubInner2 ) )->getNodeID(iNode2);
		//				const int nodeID1 = convertNodeIDBoundCurve2Triangle(nodeIDBoundCurve1);
		//				const int nodeID2 = convertNodeIDBoundCurve2Triangle(nodeIDBoundCurve2);

		//				const int iTriangle = shareSameTriangleExceptDomainBoundary( nodeID1, nodeID2 );

		//				if( iTriangle >= 0 && ( m_triangles[iTriangle].getDomainType() == CommonParameters::SEA || m_triangles[iTriangle].getDomainType() == CommonParameters::LAKE ) ){
		//					const int nd1 = nodeID1 < nodeID2 ? nodeID1 : nodeID2;
		//					const int nd2 = nodeID1 < nodeID2 ? nodeID2 : nodeID1;
		//					nodePair.push_back( std::make_pair( nd1, nd2 ) );
		//				}

		//			}
		//		}//------------------------------------------------------------------------- // Loop of sub-inner boundary 2 <<<<<

		//	}
		//}//------------------------------------------------------------------------- // Loop of sub-inner boundary 1 <<<<<
		//------------------------------------------------------------------------------------------------------ Not delete for future use <<<<<

		std::sort(nodePair.begin(), nodePair.end());
		nodePair.erase( std::unique( nodePair.begin(), nodePair.end() ), nodePair.end() );

		std::vector<CommonParameters::XY> nodesInsertedOnEdge;
		for( std::vector< std::pair<int,int> >::iterator itrNodePair = nodePair.begin(); itrNodePair != nodePair.end(); ++itrNodePair ){

				//std::cout << "Node Pair : " <<  itrNodePair->first << " " << itrNodePair->second << std::endl;

				const CommonParameters::XY coord[2] = {
					m_nodeList.getCoordXYOfPoints( itrNodePair->first ),
					m_nodeList.getCoordXYOfPoints( itrNodePair->second )
				};
				const CommonParameters::XY coordInsert = {
					0.5*( coord[0].X + coord[1].X ),
					0.5*( coord[0].Y + coord[1].Y )
				};
				nodesInsertedOnEdge.push_back( coordInsert );
	
		}

		OutputFiles::m_logFile << "# Iteration : " << std::setw(10) << icount << ", Number of nodes to be inserted on edge : " << std::setw(10) << nodesInsertedOnEdge.size() << std::endl;

		// Add new nodes
		for( std::vector<CommonParameters::XY>::iterator itr = nodesInsertedOnEdge.begin(); itr != nodesInsertedOnEdge.end(); ++itr ){
			int locType = Triangle::ON_EDGE;
			insertNewNodeAndFlip( *itr, locType, -1 );
		}
		//--------------------------------------------------------------------------------------------------

		if( nodesInsertedOnCenter.empty() && nodesInsertedOnEdge.empty() ){
			continueRefinement = false;
		}
		
		//if( nodesInsertedOnCenter.empty() ){
		//	break;
		//}

	}

}

// Perform laplacian method
void TriangleList::laplacian(){

	OutputFiles::m_logFile << "# Perform laplacian method" << std::endl;

	const double EPS = 1.0e-6;

	int maxCount = Control::getInstance()->getIterNumLaplacianMethod();

	for( int iter = 1; iter <= maxCount; ++iter ){

		OutputFiles::m_logFile << "# Iteration : " << std::setw(10) << iter << std::endl;

		const int numNode = m_nodeList.getTotalNumberOfNode();

		for( int iNode = 0; iNode < numNode; ++iNode ){

			if( getNodeIDBoundaryCurveList( iNode ) >= 0 ){// Node belong to boundary curve
				continue;
			}
		
			std::vector<int> triangles = m_nodeList.getTriangleIDs( iNode );
			std::set<int> surroundingNodes;
			for( std::vector<int>::iterator itr = triangles.begin(); itr != triangles.end(); ++itr ){
				for( int i = 0; i < 3; ++i ){
					const int anotherNode = m_triangles[*itr].getNodeID( i );
					if( anotherNode != iNode ){
						surroundingNodes.insert( anotherNode );
					}
				}
			}

			const CommonParameters::XY coordOrg = m_nodeList.getCoordXYOfPoints( iNode );

			CommonParameters::XY coordNew = { 0.0, 0.0 };
			for( std::set<int>::iterator itr = surroundingNodes.begin(); itr != surroundingNodes.end(); ++itr ){
				const CommonParameters::XY coord = m_nodeList.getCoordXYOfPoints( *itr );
				coordNew.X += coord.X;
				coordNew.Y += coord.Y;
			}
			coordNew.X /= static_cast<int>( surroundingNodes.size() );
			coordNew.Y /= static_cast<int>( surroundingNodes.size() );

			m_nodeList.setCoordOfPoints( iNode, coordNew );

			for( std::vector<int>::iterator itr = triangles.begin(); itr != triangles.end(); ++itr ){
				if( m_triangles[*itr].calcDeterminant( &m_nodeList ) < EPS ){
					m_nodeList.setCoordOfPoints( iNode, coordOrg );
					break;
				}
			}

		}

		//std::ostringstream oss;
		//oss << "surface_triangle.laplacian" << iter << ".vtk";
		//writeTrianglesToVTK( oss.str().c_str() );

	}

}

// Get node ID of boundary curve list from node ID of triangle list
int TriangleList::getNodeIDBoundaryCurveList( const int nodeIDTriangleList ) const{

	std::map<int,int>::const_iterator itr = m_nodeIDTriangle2BoundCurve.find( nodeIDTriangleList );
	if( itr == m_nodeIDTriangle2BoundCurve.end() ){
		return -1;
	}

	return itr->second;

}

// Make PLCs
void TriangleList::makePLCs( std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList ){

	for( int i = 0; i < 12; ++i ){
		m_additionalCornerNodes[i] = -1;
	}

	makeSurfMeshOfEarth( nodeList, facetList );

	if((Control::getInstance())->getFlagWhetherMakeSurfMesOnAllBoundary()){

		if( (Control::getInstance())->getIncludeSedimentLayer() ){
			makeSurfMeshForSedimentLayer( nodeList, facetList );
		}

		std::map<int,int>* outerBounNode2MeshNode(NULL);
		std::map<int,int>* innerBounNode2MeshNode(NULL);
		std::map<int,int>* subInnerBounNode2MeshNode(NULL);
		makeBoundaryCurveForSeaAndLakeMeshes( nodeList, outerBounNode2MeshNode, innerBounNode2MeshNode, subInnerBounNode2MeshNode );

		makeSurfMeshOfSea( outerBounNode2MeshNode, innerBounNode2MeshNode, subInnerBounNode2MeshNode, nodeList, facetList );

		makeFacetOnTopAndBot( nodeList, facetList );

		makeFacetOnSideOfSea( outerBounNode2MeshNode, nodeList, facetList );

		makeFacetOnSideOfAirAndLand( outerBounNode2MeshNode, CommonParameters::YZ_MINUS, nodeList, facetList );

		makeFacetOnSideOfAirAndLand( outerBounNode2MeshNode, CommonParameters::YZ_PLUS, nodeList, facetList );

		makeFacetOnSideOfAirAndLand( outerBounNode2MeshNode, CommonParameters::ZX_MINUS, nodeList, facetList );

		makeFacetOnSideOfAirAndLand( outerBounNode2MeshNode, CommonParameters::ZX_PLUS, nodeList, facetList );

		makeFacetOfLowerLayers( nodeList, facetList );

	}
	else{
		makeFacetOnInnnerSurface( nodeList, facetList );

		makeFacetOnTopAndBot( nodeList, facetList );

		makeFacetOnSideOfSea( NULL, nodeList, facetList );

		makeFacetOnSeaSurface( nodeList, facetList );

		makeFacetOnSideOfAirAndLand( NULL, CommonParameters::YZ_MINUS, nodeList, facetList );

		makeFacetOnSideOfAirAndLand( NULL, CommonParameters::YZ_PLUS, nodeList, facetList );

		makeFacetOnSideOfAirAndLand( NULL, CommonParameters::ZX_MINUS, nodeList, facetList );

		makeFacetOnSideOfAirAndLand( NULL, CommonParameters::ZX_PLUS, nodeList, facetList );

		makeFacetOfExtendedRegion( nodeList, facetList );
	}

}

// Make surface mesh of the earth
void TriangleList::makeSurfMeshOfEarth( std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList ){

	const int numNodes = m_nodeList.getTotalNumberOfNode();
	for( int iNode = 0; iNode < numNodes; ++iNode ){
		nodeList.push_back( m_nodeList.getCoordXYZOfPoints( iNode ) );
	}

	// Surface of the earth
	//const int numTriangles = static_cast<int>( m_triangles.size() );
	for( std::vector< Triangle >::const_iterator itr = m_triangles.begin(); itr != m_triangles.end(); ++itr ){

		//if( itr->getDomainType() == CommonParameters::SEA || itr->getDomainType() == CommonParameters::LAKE ){
		//	bool allLocateCoastLine(true);
		//	for( int iNode = 0; iNode < 3; ++iNode ){
		//		if( m_nodeList.getPointerToNode( itr->getNodeID(iNode) )->getLocation() != Node::COAST_LINE ){
		//			allLocateCoastLine = false;
		//			break;
		//		}
		//	}
		//	if( allLocateCoastLine ){
		//		continue;
		//	}
		//}

		TriangleList::PolygonNodes nodes;
		for( int iNode = 0; iNode < 3; ++iNode ){
			nodes.push_back( itr->getNodeID(iNode) );
		}
		TriangleList::Facet facetTmp;
		facetTmp.polygons.push_back( nodes );
		facetList.push_back( facetTmp );
	}

}

// Make PLCs for the inner surface inside of land
void TriangleList::makeFacetOnInnnerSurface( std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList ){

	const BoundaryCurveList* const ptrBoundaryCurveList = &m_boundaryCurveList;

	const int numOuterBoundary = ptrBoundaryCurveList->getTotalNumberOuterBoundaries();
	for( int iBounOuter = 0; iBounOuter < numOuterBoundary; ++iBounOuter ){

		const int numBounInnerIncluded = ptrBoundaryCurveList->getNumInnerBoundaryIncluded( iBounOuter );
		for( int iBounInnerIncluded = 0; iBounInnerIncluded < numBounInnerIncluded; ++iBounInnerIncluded ){
			const int iBounInner = ptrBoundaryCurveList->getInnerBoundaryIDIncluded(iBounOuter,iBounInnerIncluded);

			if( ( ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner ) )->getGeologicalType() != CommonParameters::SEA &&
				( ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner ) )->getGeologicalType() != CommonParameters::LAKE ){
				continue;// Treat only the sea or lakes
			}

			TriangleList::PolygonNodes nodesInnerBound;
			const int numNodesInnerBound = ptrBoundaryCurveList->getNumNodeInnerBoundary( iBounInner );
			for( int iNodeInnerBound = 0; iNodeInnerBound < numNodesInnerBound; ++iNodeInnerBound ){

				const int nodeIDInnerBound = ( ptrBoundaryCurveList->getPointerToInnerBoundary(iBounInner) )->getNodeID( iNodeInnerBound );
				nodesInnerBound.push_back( convertNodeIDBoundCurve2Triangle(nodeIDInnerBound) );

			}

			TriangleList::Facet facetTmpInnerBound;
			facetTmpInnerBound.polygons.push_back( nodesInnerBound );
			facetList.push_back( facetTmpInnerBound );

			const int numBounSubInnerIncluded = ptrBoundaryCurveList->getNumSubInnerBoundaryIncluded( iBounInner );
			for( int iBounSubInnerIncluded = 0; iBounSubInnerIncluded < numBounSubInnerIncluded; ++iBounSubInnerIncluded ){
				const int iBounSubInner = ptrBoundaryCurveList->getSubInnerBoundaryIDIncluded(iBounInner,iBounSubInnerIncluded);

				if( ( ptrBoundaryCurveList->getPointerToSubInnerBoundary( iBounSubInner ) )->getGeologicalType() != CommonParameters::SEA &&
					( ptrBoundaryCurveList->getPointerToSubInnerBoundary( iBounSubInner ) )->getGeologicalType() != CommonParameters::LAKE ){
					continue;// Treat only the sea or lakes
				}

				TriangleList::PolygonNodes nodesSubInnerBound;
				const int numNodesSubInnerBound = ptrBoundaryCurveList->getNumNodeSubInnerBoundary( iBounSubInner );
				for( int iNodeSubInnerBound = 0; iNodeSubInnerBound < numNodesSubInnerBound; ++iNodeSubInnerBound ){

					const int nodeIDSubInnerBound = ( ptrBoundaryCurveList->getPointerToSubInnerBoundary(iBounSubInner) )->getNodeID( iNodeSubInnerBound );
					nodesSubInnerBound.push_back( convertNodeIDBoundCurve2Triangle(nodeIDSubInnerBound) );

				}
				TriangleList::Facet facetTmpSubInnerBound;
				facetTmpSubInnerBound.polygons.push_back( nodesSubInnerBound );
				facetList.push_back( facetTmpSubInnerBound );

			}

		}

	}
	
}

// Make surface mesh of the sea
void TriangleList::makeSurfMeshOfSea(  const std::map<int,int>* const outerBounNode2MeshNode, const std::map<int,int>* const innerBounNode2MeshNode, const std::map<int,int>* const subInnerBounNode2MeshNode,
	std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList ){

	const BoundaryCurveList* const ptrBoundaryCurveList = &m_boundaryCurveList;

	std::map<int,int> nodesAlreadyInserted;
	//const int numNodes = m_nodeList.getTotalNumberOfNode();
	//for( int iNode = 0; iNode < numNodes; ++iNode ){
	//	nodesAlreadyInserted.insert(std::make_pair(iNode,iNode));
	//}

	//---- [Note] Not delete for future use >>>>>
	//std::set<int> triangleAreadyInserted;

	//// Outer boundary
	//const int numOuterBoundary = ptrBoundaryCurveList->getTotalNumberOuterBoundaries();
	//for( int iBoun = 0; iBoun < numOuterBoundary; ++iBoun ){
	//	const BoundaryCurveOuter* const ptrBound = ptrBoundaryCurveList->getPointerToOuterBoundary(iBoun);
	//	if( ptrBound->getGeologicalType() != CommonParameters::SEA && ptrBound->getGeologicalType() != CommonParameters::LAKE ){
	//		continue;// Treat only the sea or lakes
	//	}

	//	std::vector<int> triangles;
	//	searchTrianglesWithinBoundaryCurve( ptrBound, true, triangles );
	//	for( std::vector<int>::iterator itr = triangles.begin(); itr != triangles.end(); ++itr ){

	//		if( m_triangles[*itr].getDomainType() != CommonParameters::SEA && m_triangles[*itr].getDomainType() != CommonParameters::LAKE ){
	//			continue;
	//		}
	//		std::set<int>::const_iterator itrTriangleAlreadyInserted = triangleAreadyInserted.find(*itr);
	//		if( itrTriangleAlreadyInserted != triangleAreadyInserted.end() ){// Found
	//			continue;
	//		}

	//		TriangleList::PolygonNodes nodes;
	//		for( int iNode = 0; iNode < 3; ++iNode ){
	//			const int nodeIDTriangle = m_triangles[*itr].getNodeID(iNode);
	//			//if( m_nodeList.getPointerToNode(nodeIDTriangle)->getLocation() == Node::COAST_LINE ){
	//			//	const int nodeIDBoun = convertNodeIDTriangle2BoundCurve(nodeIDTriangle);
	//			//	std::map<int,int>::const_iterator itrOuterBounNode2MeshNode = outerBounNode2MeshNode[iBoun].find(nodeIDBoun);
	//			//	if( itrOuterBounNode2MeshNode != outerBounNode2MeshNode[iBoun].end() ){
	//			//		nodes.push_back(itrOuterBounNode2MeshNode->second);
	//			//	}
	//			//	else{
	//			//		OutputFiles::m_logFile << " Error : Node " << nodeIDBoun << " cannot be found in outerBounNode2MeshNode[" << iBoun << "] ." << std::endl;
	//			//		exit(1);
	//			//	}
	//			//}				
	//			const int nodeIDBoun = convertNodeIDTriangle2BoundCurve(nodeIDTriangle);
	//			if( nodeIDBoun >= 0 ){
	//				std::map<int,int>::const_iterator itrOuterBounNode2MeshNode = outerBounNode2MeshNode[iBoun].find(nodeIDBoun);
	//				if( itrOuterBounNode2MeshNode != outerBounNode2MeshNode[iBoun].end() ){
	//					nodes.push_back(itrOuterBounNode2MeshNode->second);
	//				}				
	//				else{
	//					OutputFiles::m_logFile << " Error : Node " << nodeIDBoun << " cannot be found in outerBounNode2MeshNode[" << iBoun << "] ." << std::endl;
	//					exit(1);
	//				}
	//			}
	//			else{
	//				std::map<int,int>::const_iterator itrNodesAlreadyInserted = nodesAlreadyInserted.find(nodeIDTriangle);
	//				if( itrNodesAlreadyInserted != nodesAlreadyInserted.end() ){// Found
	//					nodes.push_back( itrNodesAlreadyInserted->second );
	//				}
	//				else{// Not found
	//					const CommonParameters::XY coordXY = m_nodeList.getPointerToNode(nodeIDTriangle)->getCoordXY();
	//					const CommonParameters::XYZ coordXYZ = { coordXY.X, coordXY.Y, 0.0 };
	//					nodeList.push_back( coordXYZ );
	//					nodes.push_back(static_cast<int>(nodeList.size())-1);
	//					nodesAlreadyInserted.insert(std::make_pair(nodeIDTriangle,static_cast<int>(nodeList.size())-1));
	//				}
	//			}
	//		}
	//		TriangleList::Facet facetTmp;
	//		facetTmp.polygons.push_back( nodes );
	//		facetList.push_back( facetTmp );
	//		triangleAreadyInserted.insert(*itr);
	//	}
	//}

	//// Inner boundary
	//const int numInnerBoundary = ptrBoundaryCurveList->getTotalNumberInnerBoundaries();
	//for( int iBoun = 0; iBoun < numInnerBoundary; ++iBoun ){
	//	const BoundaryCurveInner* const ptrBound = ptrBoundaryCurveList->getPointerToInnerBoundary(iBoun);
	//	if( ptrBound->getGeologicalType() != CommonParameters::SEA && ptrBound->getGeologicalType() != CommonParameters::LAKE ){
	//		continue;// Treat only the sea or lakes
	//	}

	//	std::vector<int> triangles;
	//	searchTrianglesWithinBoundaryCurve( ptrBound, false, triangles );
	//	for( std::vector<int>::iterator itr = triangles.begin(); itr != triangles.end(); ++itr ){

	//		if( m_triangles[*itr].getDomainType() != CommonParameters::SEA && m_triangles[*itr].getDomainType() != CommonParameters::LAKE ){
	//			continue;
	//		}
	//		std::set<int>::const_iterator itrTriangleAlreadyInserted = triangleAreadyInserted.find(*itr);
	//		if( itrTriangleAlreadyInserted != triangleAreadyInserted.end() ){// Found
	//			continue;
	//		}

	//		TriangleList::PolygonNodes nodes;
	//		for( int iNode = 0; iNode < 3; ++iNode ){
	//			const int nodeIDTriangle = m_triangles[*itr].getNodeID(iNode);
	//			//if( m_nodeList.getPointerToNode(nodeIDTriangle)->getLocation() == Node::COAST_LINE ){
	//			//	const int nodeIDBoun = convertNodeIDTriangle2BoundCurve(nodeIDTriangle);
	//			//	std::map<int,int>::const_iterator itrInnerBounNode2MeshNode = innerBounNode2MeshNode[iBoun].find(nodeIDBoun);
	//			//	if( itrInnerBounNode2MeshNode != innerBounNode2MeshNode[iBoun].end() ){
	//			//		nodes.push_back(itrInnerBounNode2MeshNode->second);
	//			//	}
	//			//	else{
	//			//		OutputFiles::m_logFile << " Error : Node " << nodeIDBoun << " cannot be found in innerBounNode2MeshNode[" << iBoun << "] ." << std::endl;
	//			//		exit(1);
	//			//	}
	//			//}
	//			const int nodeIDBoun = convertNodeIDTriangle2BoundCurve(nodeIDTriangle);
	//			if( nodeIDBoun >= 0 ){
	//				std::map<int,int>::const_iterator itrInnerBounNode2MeshNode = innerBounNode2MeshNode[iBoun].find(nodeIDBoun);
	//				if( itrInnerBounNode2MeshNode != innerBounNode2MeshNode[iBoun].end() ){
	//					nodes.push_back(itrInnerBounNode2MeshNode->second);
	//				}
	//				else{
	//					OutputFiles::m_logFile << " Error : Node " << nodeIDBoun << " cannot be found in innerBounNode2MeshNode[" << iBoun << "] ." << std::endl;
	//					exit(1);
	//				}
	//			}
	//			else{
	//				std::map<int,int>::const_iterator itrNodesAlreadyInserted = nodesAlreadyInserted.find(nodeIDTriangle);
	//				if( itrNodesAlreadyInserted != nodesAlreadyInserted.end() ){// Found
	//					nodes.push_back( itrNodesAlreadyInserted->second );
	//				}
	//				else{// Not found
	//					const CommonParameters::XY coordXY = m_nodeList.getPointerToNode(nodeIDTriangle)->getCoordXY();
	//					const CommonParameters::XYZ coordXYZ = { coordXY.X, coordXY.Y, 0.0 };
	//					nodeList.push_back( coordXYZ );
	//					nodes.push_back(static_cast<int>(nodeList.size())-1);
	//					nodesAlreadyInserted.insert(std::make_pair(nodeIDTriangle,static_cast<int>(nodeList.size())-1));
	//				}
	//			}
	//		}
	//		TriangleList::Facet facetTmp;
	//		facetTmp.polygons.push_back( nodes );
	//		facetList.push_back( facetTmp );
	//		triangleAreadyInserted.insert(*itr);
	//	}
	//}

	//// Sub-inner boundary
	//const int numSubInnerBoundary = ptrBoundaryCurveList->getTotalNumberSubInnerBoundaries();
	//for( int iBoun = 0; iBoun < numSubInnerBoundary; ++iBoun ){
	//	const BoundaryCurveSubInner* const ptrBound = ptrBoundaryCurveList->getPointerToSubInnerBoundary(iBoun);
	//	if( ptrBound->getGeologicalType() != CommonParameters::SEA && ptrBound->getGeologicalType() != CommonParameters::LAKE ){
	//		continue;// Treat only the sea or lakes
	//	}

	//	std::vector<int> triangles;
	//	searchTrianglesWithinBoundaryCurve( ptrBound, false, triangles );
	//	for( std::vector<int>::iterator itr = triangles.begin(); itr != triangles.end(); ++itr ){

	//		if( m_triangles[*itr].getDomainType() != CommonParameters::SEA && m_triangles[*itr].getDomainType() != CommonParameters::LAKE ){
	//			continue;
	//		}
	//		std::set<int>::const_iterator itrTriangleAlreadyInserted = triangleAreadyInserted.find(*itr);
	//		if( itrTriangleAlreadyInserted != triangleAreadyInserted.end() ){// Found
	//			continue;
	//		}

	//		TriangleList::PolygonNodes nodes;
	//		for( int iNode = 0; iNode < 3; ++iNode ){
	//			const int nodeIDTriangle = m_triangles[*itr].getNodeID(iNode);
	//			//if( m_nodeList.getPointerToNode(nodeIDTriangle)->getLocation() == Node::COAST_LINE ){
	//			//	const int nodeIDBoun = convertNodeIDTriangle2BoundCurve(nodeIDTriangle);
	//			//	std::map<int,int>::const_iterator itrSubInnerBounNode2MeshNode = subInnerBounNode2MeshNode[iBoun].find(nodeIDBoun);
	//			//	if( itrSubInnerBounNode2MeshNode != subInnerBounNode2MeshNode[iBoun].end() ){
	//			//		nodes.push_back(itrSubInnerBounNode2MeshNode->second);
	//			//	}
	//			//	else{
	//			//		OutputFiles::m_logFile << " Error : Node " << nodeIDBoun << " cannot be found in subInnerBounNode2MeshNode[" << iBoun << "] ." << std::endl;
	//			//		exit(1);
	//			//	}
	//			//}
	//			const int nodeIDBoun = convertNodeIDTriangle2BoundCurve(nodeIDTriangle);
	//			if( nodeIDBoun >= 0 ){
	//				std::map<int,int>::const_iterator itrSubInnerBounNode2MeshNode = subInnerBounNode2MeshNode[iBoun].find(nodeIDBoun);
	//				if( itrSubInnerBounNode2MeshNode != subInnerBounNode2MeshNode[iBoun].end() ){
	//					nodes.push_back(itrSubInnerBounNode2MeshNode->second);
	//				}
	//				else{
	//					OutputFiles::m_logFile << " Error : Node " << nodeIDBoun << " cannot be found in subInnerBounNode2MeshNode[" << iBoun << "] ." << std::endl;
	//					exit(1);
	//				}
	//			}
	//			else{
	//				std::map<int,int>::const_iterator itrNodesAlreadyInserted = nodesAlreadyInserted.find(nodeIDTriangle);
	//				if( itrNodesAlreadyInserted != nodesAlreadyInserted.end() ){// Found
	//					nodes.push_back( itrNodesAlreadyInserted->second );
	//				}
	//				else{// Not found
	//					const CommonParameters::XY coordXY = m_nodeList.getPointerToNode(nodeIDTriangle)->getCoordXY();
	//					const CommonParameters::XYZ coordXYZ = { coordXY.X, coordXY.Y, 0.0 };
	//					nodeList.push_back( coordXYZ );
	//					nodes.push_back(static_cast<int>(nodeList.size())-1);
	//					nodesAlreadyInserted.insert(std::make_pair(nodeIDTriangle,static_cast<int>(nodeList.size())-1));
	//				}
	//			}
	//		}
	//		TriangleList::Facet facetTmp;
	//		facetTmp.polygons.push_back( nodes );
	//		facetList.push_back( facetTmp );
	//		triangleAreadyInserted.insert(*itr);
	//	}
	//}
	//---- [Note] Not delete for future use <<<<<

	const LakeList* const ptrLakeList = LakeList::getInstance(); 

	for( std::vector< Triangle >::const_iterator itr = m_triangles.begin(); itr != m_triangles.end(); ++itr ){

		const int domainType = itr->getDomainType();
		if( domainType != CommonParameters::SEA && domainType != CommonParameters::LAKE ){
			continue;
		}

		TriangleList::PolygonNodes nodes;
		for( int iNode = 0; iNode < 3; ++iNode ){
			const int nodeIDTriangle = itr->getNodeID(iNode);
			const int nodeIDBoun = convertNodeIDTriangle2BoundCurve(nodeIDTriangle);
			if( nodeIDBoun >= 0 ){

				bool found(false);

				// Outer boundary
				const int numOuterBoundary = ptrBoundaryCurveList->getTotalNumberOuterBoundaries();
				for( int iBoun = 0; iBoun < numOuterBoundary; ++iBoun ){
					const BoundaryCurveOuter* const ptrBoundOuter = ptrBoundaryCurveList->getPointerToOuterBoundary(iBoun);
					std::map<int,int>::const_iterator itrOuterBounNode2MeshNode = outerBounNode2MeshNode[iBoun].find(nodeIDBoun);
					if( itrOuterBounNode2MeshNode != outerBounNode2MeshNode[iBoun].end() ){
						nodes.push_back(itrOuterBounNode2MeshNode->second);
						found = true;
						break;
					}		
				}
				if(found){
					continue;
				}

				// Inner boundary
				const int numInnerBoundary = ptrBoundaryCurveList->getTotalNumberInnerBoundaries();
				for( int iBoun = 0; iBoun < numInnerBoundary; ++iBoun ){
					const BoundaryCurveInner* const ptrBound = ptrBoundaryCurveList->getPointerToInnerBoundary(iBoun);
					std::map<int,int>::const_iterator itrInnerBounNode2MeshNode = innerBounNode2MeshNode[iBoun].find(nodeIDBoun);
					if( itrInnerBounNode2MeshNode != innerBounNode2MeshNode[iBoun].end() ){
						nodes.push_back(itrInnerBounNode2MeshNode->second);
						found = true;
						break;
					}
				}
				if(found){
					continue;
				}

				// Sub-inner boundary
				const int numSubInnerBoundary = ptrBoundaryCurveList->getTotalNumberSubInnerBoundaries();
				for( int iBoun = 0; iBoun < numSubInnerBoundary; ++iBoun ){
					const BoundaryCurveSubInner* const ptrBound = ptrBoundaryCurveList->getPointerToSubInnerBoundary(iBoun);
					std::map<int,int>::const_iterator itrSubInnerBounNode2MeshNode = subInnerBounNode2MeshNode[iBoun].find(nodeIDBoun);
					if( itrSubInnerBounNode2MeshNode != subInnerBounNode2MeshNode[iBoun].end() ){
						nodes.push_back(itrSubInnerBounNode2MeshNode->second);
						found = true;
						break;
					}
				}
				if(found){
					continue;
				}

				OutputFiles::m_logFile << " Error : Node " << nodeIDBoun << " cannot be found in arrays relating nodes of boudary and the ones of mesh." << std::endl;
				exit(1);
			}
			else{
				std::map<int,int>::const_iterator itrNodesAlreadyInserted = nodesAlreadyInserted.find(nodeIDTriangle);
				if( itrNodesAlreadyInserted != nodesAlreadyInserted.end() ){// Found
					nodes.push_back( itrNodesAlreadyInserted->second );
				}
				else{// Not found
					const CommonParameters::XY coordXY = m_nodeList.getPointerToNode(nodeIDTriangle)->getCoordXY();
					double coordZ(0.0);
					if( domainType == CommonParameters::LAKE ){
						const int lakeIndex = ptrBoundaryCurveList->getLakeIndexFromCoordinate( coordXY );
						coordZ = - ptrLakeList->getLakeHeight( lakeIndex );
					}
					const CommonParameters::XYZ coordXYZ = { coordXY.X, coordXY.Y, coordZ };
					nodeList.push_back( coordXYZ );
					nodes.push_back(static_cast<int>(nodeList.size())-1);
					nodesAlreadyInserted.insert(std::make_pair(nodeIDTriangle,static_cast<int>(nodeList.size())-1));
				}
			}
		}
		TriangleList::Facet facetTmp;
		facetTmp.polygons.push_back( nodes );
		facetList.push_back( facetTmp );
	}

}

// Make surface mesh for sediment layer
void TriangleList::makeSurfMeshForSedimentLayer( std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList ){

	const bool makeSurfMesOnAllBoundary = (Control::getInstance())->getFlagWhetherMakeSurfMesOnAllBoundary();
	if( !makeSurfMesOnAllBoundary ){
		OutputFiles::m_logFile << " Error : Sediment layer can be included if and only if surface meshes are made on all boundary faces." << std::endl;
		exit(1);
	}

	const double thicknessSediment = (Control::getInstance())->getThicknessSediment();

	std::vector<TriangleList::Facet>::const_iterator itrLandSurf;
	for( std::vector< Triangle >::const_iterator itr = m_triangles.begin(); itr != m_triangles.end(); ++itr ){
		if( itr->getDomainType() != CommonParameters::SEA && itr->getDomainType() != CommonParameters::LAKE ){
			continue;
		}

		TriangleList::PolygonNodes nodes;
		for( int iNode = 0; iNode < 3; ++iNode ){
			const int nodeIDTriangle = itr->getNodeID(iNode);

			std::map<int,int>::const_iterator itrNodesAlreadyInserted = m_nodeLandSurfToSediment.find(nodeIDTriangle);
			if( itrNodesAlreadyInserted != m_nodeLandSurfToSediment.end() ){// Found
				nodes.push_back( itrNodesAlreadyInserted->second );
			}
			else{// Not found
				CommonParameters::XYZ coordXYZ = m_nodeList.getPointerToNode(nodeIDTriangle)->getCoordXYZ();
				coordXYZ.Z += thicknessSediment;
				nodeList.push_back( coordXYZ );
				const int nodeIDSediment = static_cast<int>( nodeList.size() ) - 1;
				nodes.push_back(nodeIDSediment);
				m_nodeLandSurfToSediment.insert( std::make_pair( nodeIDTriangle, nodeIDSediment ) );
			}
		}

		TriangleList::Facet facetTmp;
		facetTmp.polygons.push_back( nodes );
		facetList.push_back( facetTmp );
	}

	//-----------------------------------
	//--- Side of the sediment layers ---
	//-----------------------------------
	const BoundaryCurveList* const ptrBoundaryCurveList = &m_boundaryCurveList;

	//--- Outer boundary ---
	const int numOuterBoundary = ptrBoundaryCurveList->getTotalNumberOuterBoundaries();
	for( int iBounOuter = 0; iBounOuter < numOuterBoundary; ++iBounOuter ){

		if( ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getGeologicalType() != CommonParameters::SEA &&
			( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getGeologicalType() != CommonParameters::LAKE ){
			continue;// Treat only the sea or lakes
		}

		const int numNodes = ptrBoundaryCurveList->getNumNodeOuterBoundary( iBounOuter );
		
		const int nodeIDBoundCurveListLast = ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getNodeID(numNodes - 1);
		int nodeIDPre = convertNodeIDBoundCurve2Triangle( nodeIDBoundCurveListLast );
		for( int iNode = 0; iNode < numNodes; ++iNode ){
			const int nodeIDBoundCurveList = ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getNodeID(iNode);
			const int nodeID = convertNodeIDBoundCurve2Triangle( nodeIDBoundCurveList );			
			{// First triangle
				TriangleList::PolygonNodes nodes1;
				nodes1.push_back( nodeIDPre );
				nodes1.push_back( nodeID );
				nodes1.push_back( m_nodeLandSurfToSediment[nodeID] );
				TriangleList::Facet facetTmp1;
				facetTmp1.polygons.push_back( nodes1 );
				facetList.push_back( facetTmp1 );
			}
			{// Second triangle
				TriangleList::PolygonNodes nodes2;
				nodes2.push_back( nodeIDPre );
				nodes2.push_back( m_nodeLandSurfToSediment[nodeIDPre] );
				nodes2.push_back( m_nodeLandSurfToSediment[nodeID] );
				TriangleList::Facet facetTmp2;
				facetTmp2.polygons.push_back( nodes2 );
				facetList.push_back( facetTmp2 );
			}
			nodeIDPre = nodeID;
		}
	}

	//--- Inner boundary ---
	const int numInnerBoundary = ptrBoundaryCurveList->getTotalNumberInnerBoundaries();
	for( int iBounInner = 0; iBounInner < numInnerBoundary; ++iBounInner ){

		if( ( ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner ) )->getGeologicalType() != CommonParameters::LAND ){
			continue;// Treat only the inner land boundary
		}

		const int numNodes = ptrBoundaryCurveList->getNumNodeInnerBoundary( iBounInner );

		const int nodeIDBoundCurveListLast = ( ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner ) )->getNodeID(numNodes - 1);
		int nodeIDPre = convertNodeIDBoundCurve2Triangle( nodeIDBoundCurveListLast );
		for( int iNode = 0; iNode < numNodes; ++iNode ){
			const int nodeIDBoundCurveList = ( ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner ) )->getNodeID(iNode);
			const int nodeID = convertNodeIDBoundCurve2Triangle( nodeIDBoundCurveList );			
			{// First triangle
				TriangleList::PolygonNodes nodes1;
				nodes1.push_back( nodeIDPre );
				nodes1.push_back( nodeID );
				nodes1.push_back( m_nodeLandSurfToSediment[nodeID] );
				TriangleList::Facet facetTmp1;
				facetTmp1.polygons.push_back( nodes1 );
				facetList.push_back( facetTmp1 );
			}
			{// Second triangle
				TriangleList::PolygonNodes nodes2;
				nodes2.push_back( nodeIDPre );
				nodes2.push_back( m_nodeLandSurfToSediment[nodeIDPre] );
				nodes2.push_back( m_nodeLandSurfToSediment[nodeID] );
				TriangleList::Facet facetTmp2;
				facetTmp2.polygons.push_back( nodes2 );
				facetList.push_back( facetTmp2 );
			}
			nodeIDPre = nodeID;
		}
	}

	//--- Subinner boundary ---
	const int numSubInnerBoundary = ptrBoundaryCurveList->getTotalNumberSubInnerBoundaries();
	for( int iBounSubInner = 0; iBounSubInner < numSubInnerBoundary; ++iBounSubInner ){

		if( ( ptrBoundaryCurveList->getPointerToSubInnerBoundary( iBounSubInner ) )->getGeologicalType() != CommonParameters::SEA &&
			( ptrBoundaryCurveList->getPointerToSubInnerBoundary( iBounSubInner ) )->getGeologicalType() != CommonParameters::LAKE ){
			continue;// Treat only the sea or lakes
		}

		const int numNodes = ptrBoundaryCurveList->getNumNodeSubInnerBoundary( iBounSubInner );

		const int nodeIDBoundCurveListLast = ( ptrBoundaryCurveList->getPointerToSubInnerBoundary( iBounSubInner ) )->getNodeID(numNodes - 1);
		int nodeIDPre = convertNodeIDBoundCurve2Triangle( nodeIDBoundCurveListLast );
		for( int iNode = 0; iNode < numNodes; ++iNode ){
			const int nodeIDBoundCurveList = ( ptrBoundaryCurveList->getPointerToSubInnerBoundary( iBounSubInner ) )->getNodeID(iNode);
			const int nodeID = convertNodeIDBoundCurve2Triangle( nodeIDBoundCurveList );			
			{// First triangle
				TriangleList::PolygonNodes nodes1;
				nodes1.push_back( nodeIDPre );
				nodes1.push_back( nodeID );
				nodes1.push_back( m_nodeLandSurfToSediment[nodeID] );
				TriangleList::Facet facetTmp1;
				facetTmp1.polygons.push_back( nodes1 );
				facetList.push_back( facetTmp1 );
			}
			{// Second triangle
				TriangleList::PolygonNodes nodes2;
				nodes2.push_back( nodeIDPre );
				nodes2.push_back( m_nodeLandSurfToSediment[nodeIDPre] );
				nodes2.push_back( m_nodeLandSurfToSediment[nodeID] );
				TriangleList::Facet facetTmp2;
				facetTmp2.polygons.push_back( nodes2 );
				facetList.push_back( facetTmp2 );
			}
			nodeIDPre = nodeID;
		}
	}

}

// Make PLCs for the top and bottom
void TriangleList::makeFacetOnTopAndBot( std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList ){

	const AnalysisDomain* const ptrAnalysisDomain = AnalysisDomain::getInstance();
	const bool isLayers = (Control::getInstance())->getIncludeLayers();

	const double xMin = ptrAnalysisDomain->getMinCoordX();
	const double xMax = ptrAnalysisDomain->getMaxCoordX(); 
	const double yMin = ptrAnalysisDomain->getMinCoordY(); 
	const double yMax = ptrAnalysisDomain->getMaxCoordY(); 
	const double zMin = ptrAnalysisDomain->getMinCoordZ(); 
	const double zMax = isLayers ? (Control::getInstance())->getDepthLayerInterfaces(0) : ptrAnalysisDomain->getMaxCoordZ(); 

	const CommonParameters::XYZ coordXPlusYMinusTop  = { xMax, yMin, zMin };
	const CommonParameters::XYZ coordXPlusYPlusTop   = { xMax, yMax, zMin };
	const CommonParameters::XYZ coordXMinusYPlusTop  = { xMin, yMax, zMin };
	const CommonParameters::XYZ coordXMinusYMinusTop = { xMin, yMin, zMin };

	const CommonParameters::XYZ coordXPlusYMinusBot  = { xMax, yMin, zMax };
	const CommonParameters::XYZ coordXPlusYPlusBot   = { xMax, yMax, zMax };
	const CommonParameters::XYZ coordXMinusYPlusBot  = { xMin, yMax, zMax };
	const CommonParameters::XYZ coordXMinusYMinusBot = { xMin, yMin, zMax };

	int nodeIDCur = static_cast<int>( nodeList.size() );

	if((Control::getInstance())->getFlagWhetherMakeSurfMesOnAllBoundary()){

		nodeList.push_back( coordXPlusYMinusTop );
		m_additionalCornerNodes[TriangleList::XPLUS_YMINUS_TOP]  = nodeIDCur++; 
		nodeList.push_back( coordXPlusYPlusTop );
		m_additionalCornerNodes[TriangleList::XPLUS_YPLUS_TOP]   = nodeIDCur++; 
		nodeList.push_back( coordXMinusYPlusTop );
		m_additionalCornerNodes[TriangleList::XMINUS_YPLUS_TOP]  = nodeIDCur++; 
		nodeList.push_back( coordXMinusYMinusTop );
		m_additionalCornerNodes[TriangleList::XMINUS_YMINUS_TOP] = nodeIDCur++;

		// X Minus
		std::vector<double> terms;
		Util::calculateTermsOfGeometricProgressionExceptLast(
			(Control::getInstance())->calcMaximumEdgeLength(coordXMinusYMinusTop),
			(Control::getInstance())->calcMaximumEdgeLength(coordXMinusYPlusTop),
			yMax - yMin, terms );
		for( std::vector<double>::const_iterator itr = terms.begin(); itr != terms.end(); ++itr ){
			const CommonParameters::XYZ coord = { xMin, yMin + *itr, zMin };
			nodeList.push_back( coord );
			m_additionalSideNodes[TriangleList::XMINUS_TOP].push_back(nodeIDCur++);
		}

		// X Plus
		terms.clear();
		Util::calculateTermsOfGeometricProgressionExceptLast(
			(Control::getInstance())->calcMaximumEdgeLength(coordXPlusYMinusTop),
			(Control::getInstance())->calcMaximumEdgeLength(coordXPlusYPlusTop),
			yMax - yMin, terms );
		for( std::vector<double>::const_iterator itr = terms.begin(); itr != terms.end(); ++itr ){
			const CommonParameters::XYZ coord = { xMax, yMin + *itr, zMin };
			nodeList.push_back( coord );
			m_additionalSideNodes[TriangleList::XPLUS_TOP].push_back(nodeIDCur++);
		}

		// Y Minus
		terms.clear();
		Util::calculateTermsOfGeometricProgressionExceptLast(
			(Control::getInstance())->calcMaximumEdgeLength(coordXMinusYMinusTop),
			(Control::getInstance())->calcMaximumEdgeLength(coordXPlusYMinusTop),
			xMax - xMin, terms );
		for( std::vector<double>::const_iterator itr = terms.begin(); itr != terms.end(); ++itr ){
			const CommonParameters::XYZ coord = { xMin + *itr, yMin, zMin };
			nodeList.push_back( coord );
			m_additionalSideNodes[TriangleList::YMINUS_TOP].push_back(nodeIDCur++);
		}

		// Y Plus
		terms.clear();
		Util::calculateTermsOfGeometricProgressionExceptLast(
			(Control::getInstance())->calcMaximumEdgeLength(coordXMinusYPlusTop),
			(Control::getInstance())->calcMaximumEdgeLength(coordXPlusYPlusTop),
			xMax - xMin, terms );
		for( std::vector<double>::const_iterator itr = terms.begin(); itr != terms.end(); ++itr ){
			const CommonParameters::XYZ coord = { xMin + *itr, yMax, zMin };
			nodeList.push_back( coord );
			m_additionalSideNodes[TriangleList::YPLUS_TOP].push_back(nodeIDCur++);
		}

		TriangleList::PolygonNodes nodesTop;
		nodesTop.push_back(m_additionalCornerNodes[TriangleList::XMINUS_YMINUS_TOP]);
		for( std::vector<int>::const_iterator itr = m_additionalSideNodes[TriangleList::YMINUS_TOP].begin(); itr != m_additionalSideNodes[TriangleList::YMINUS_TOP].end(); ++itr ){
			nodesTop.push_back(*itr);
		}
		nodesTop.push_back(m_additionalCornerNodes[TriangleList::XPLUS_YMINUS_TOP]);
		for( std::vector<int>::const_iterator itr = m_additionalSideNodes[TriangleList::XPLUS_TOP].begin(); itr != m_additionalSideNodes[TriangleList::XPLUS_TOP].end(); ++itr ){
			nodesTop.push_back(*itr);
		}
		nodesTop.push_back(m_additionalCornerNodes[TriangleList::XPLUS_YPLUS_TOP]);
		for( std::vector<int>::const_reverse_iterator itr = m_additionalSideNodes[TriangleList::YPLUS_TOP].rbegin(); itr != m_additionalSideNodes[TriangleList::YPLUS_TOP].rend(); ++itr ){
			nodesTop.push_back(*itr);
		}
		nodesTop.push_back(m_additionalCornerNodes[TriangleList::XMINUS_YPLUS_TOP]);
		for( std::vector<int>::const_reverse_iterator itr = m_additionalSideNodes[TriangleList::XMINUS_TOP].rbegin(); itr != m_additionalSideNodes[TriangleList::XMINUS_TOP].rend(); ++itr ){
			nodesTop.push_back(*itr);
		}
	
		nodeList.push_back( coordXPlusYMinusBot );
		m_additionalCornerNodes[TriangleList::XPLUS_YMINUS_BOTTOM]  = nodeIDCur++; 
		nodeList.push_back( coordXPlusYPlusBot );
		m_additionalCornerNodes[TriangleList::XPLUS_YPLUS_BOTTOM]   = nodeIDCur++; 
		nodeList.push_back( coordXMinusYPlusBot );
		m_additionalCornerNodes[TriangleList::XMINUS_YPLUS_BOTTOM]  = nodeIDCur++; 
		nodeList.push_back( coordXMinusYMinusBot );
		m_additionalCornerNodes[TriangleList::XMINUS_YMINUS_BOTTOM] = nodeIDCur++; 

		// X Minus
		terms.clear();
		Util::calculateTermsOfGeometricProgressionExceptLast(
			(Control::getInstance())->calcMaximumEdgeLength(coordXMinusYMinusTop),
			(Control::getInstance())->calcMaximumEdgeLength(coordXMinusYPlusTop),
			yMax - yMin, terms );
		for( std::vector<double>::const_iterator itr = terms.begin(); itr != terms.end(); ++itr ){
			const CommonParameters::XYZ coord = { xMin, yMin + *itr, zMax };
			nodeList.push_back( coord );
			m_additionalSideNodes[TriangleList::XMINUS_BOTTOM].push_back(nodeIDCur++);
		}

		// X Plus
		terms.clear();
		Util::calculateTermsOfGeometricProgressionExceptLast(
			(Control::getInstance())->calcMaximumEdgeLength(coordXPlusYMinusTop),
			(Control::getInstance())->calcMaximumEdgeLength(coordXPlusYPlusTop),
			yMax - yMin, terms );
		for( std::vector<double>::const_iterator itr = terms.begin(); itr != terms.end(); ++itr ){
			const CommonParameters::XYZ coord = { xMax, yMin + *itr, zMax };
			nodeList.push_back( coord );
			m_additionalSideNodes[TriangleList::XPLUS_BOTTOM].push_back(nodeIDCur++);
		}

		// Y Minus
		terms.clear();
		Util::calculateTermsOfGeometricProgressionExceptLast(
			(Control::getInstance())->calcMaximumEdgeLength(coordXMinusYMinusTop),
			(Control::getInstance())->calcMaximumEdgeLength(coordXPlusYMinusTop),
			xMax - xMin, terms );
		for( std::vector<double>::const_iterator itr = terms.begin(); itr != terms.end(); ++itr ){
			const CommonParameters::XYZ coord = { xMin + *itr, yMin, zMax };
			nodeList.push_back( coord );
			m_additionalSideNodes[TriangleList::YMINUS_BOTTOM].push_back(nodeIDCur++);
		}

		// Y Plus
		terms.clear();
		Util::calculateTermsOfGeometricProgressionExceptLast(
			(Control::getInstance())->calcMaximumEdgeLength(coordXMinusYPlusTop),
			(Control::getInstance())->calcMaximumEdgeLength(coordXPlusYPlusTop),
			xMax - xMin, terms );
		for( std::vector<double>::const_iterator itr = terms.begin(); itr != terms.end(); ++itr ){
			const CommonParameters::XYZ coord = { xMin + *itr, yMax, zMax };
			nodeList.push_back( coord );
			m_additionalSideNodes[TriangleList::YPLUS_BOTTOM].push_back(nodeIDCur++);
		}

		TriangleList::PolygonNodes nodesBot;
		nodesBot.push_back(m_additionalCornerNodes[TriangleList::XMINUS_YMINUS_BOTTOM]);
		for( std::vector<int>::const_iterator itr = m_additionalSideNodes[TriangleList::YMINUS_BOTTOM].begin(); itr != m_additionalSideNodes[TriangleList::YMINUS_BOTTOM].end(); ++itr ){
			nodesBot.push_back(*itr);
		}
		nodesBot.push_back(m_additionalCornerNodes[TriangleList::XPLUS_YMINUS_BOTTOM]);
		for( std::vector<int>::const_iterator itr = m_additionalSideNodes[TriangleList::XPLUS_BOTTOM].begin(); itr != m_additionalSideNodes[TriangleList::XPLUS_BOTTOM].end(); ++itr ){
			nodesBot.push_back(*itr);
		}
		nodesBot.push_back(m_additionalCornerNodes[TriangleList::XPLUS_YPLUS_BOTTOM]);
		for( std::vector<int>::const_reverse_iterator itr = m_additionalSideNodes[TriangleList::YPLUS_BOTTOM].rbegin(); itr != m_additionalSideNodes[TriangleList::YPLUS_BOTTOM].rend(); ++itr ){
			nodesBot.push_back(*itr);
		}
		nodesBot.push_back(m_additionalCornerNodes[TriangleList::XMINUS_YPLUS_BOTTOM]);
		for( std::vector<int>::const_reverse_iterator itr = m_additionalSideNodes[TriangleList::XMINUS_BOTTOM].rbegin(); itr != m_additionalSideNodes[TriangleList::XMINUS_BOTTOM].rend(); ++itr ){
			nodesBot.push_back(*itr);
		}
	
		makeFacetFromSurroundingNodes( nodesTop, false, CommonParameters::TOP, CommonParameters::LAND, nodeList, facetList );
		if( isLayers ){
			makeFacetFromSurroundingNodes( nodesBot, false, CommonParameters::LAYER, CommonParameters::LAND, nodeList, facetList, 0 );
		}else{
			makeFacetFromSurroundingNodes( nodesBot, false, CommonParameters::BOT, CommonParameters::LAND, nodeList, facetList );
		}
	}
	else{
		nodeList.push_back( coordXPlusYMinusTop );
		m_additionalCornerNodes[TriangleList::XPLUS_YMINUS_TOP]  = nodeIDCur; 
		nodeList.push_back( coordXPlusYPlusTop );
		m_additionalCornerNodes[TriangleList::XPLUS_YPLUS_TOP]   = nodeIDCur + 1; 
		nodeList.push_back( coordXMinusYPlusTop );
		m_additionalCornerNodes[TriangleList::XMINUS_YPLUS_TOP]  = nodeIDCur + 2; 
		nodeList.push_back( coordXMinusYMinusTop );
		m_additionalCornerNodes[TriangleList::XMINUS_YMINUS_TOP] = nodeIDCur + 3;

		nodeList.push_back( coordXPlusYMinusBot );
		m_additionalCornerNodes[TriangleList::XPLUS_YMINUS_BOTTOM]  = nodeIDCur + 4; 
		nodeList.push_back( coordXPlusYPlusBot );
		m_additionalCornerNodes[TriangleList::XPLUS_YPLUS_BOTTOM]   = nodeIDCur + 5; 
		nodeList.push_back( coordXMinusYPlusBot );
		m_additionalCornerNodes[TriangleList::XMINUS_YPLUS_BOTTOM]  = nodeIDCur + 6; 
		nodeList.push_back( coordXMinusYMinusBot );
		m_additionalCornerNodes[TriangleList::XMINUS_YMINUS_BOTTOM] = nodeIDCur + 7; 

		TriangleList::PolygonNodes nodesTop;
		for( int i = 0; i < 4; ++i ){
			nodesTop.push_back( nodeIDCur + i );
		}

		TriangleList::Facet facetTmpTop;
		facetTmpTop.polygons.push_back( nodesTop );
		facetList.push_back( facetTmpTop );

		TriangleList::PolygonNodes nodesBot;
		for( int i = 4; i < 8; ++i ){
			nodesBot.push_back( nodeIDCur + i );
		}

		TriangleList::Facet facetTmpBot;
		facetTmpBot.polygons.push_back( nodesBot );
		facetList.push_back( facetTmpBot );
	}
}

// Make PLCs for the sides of the sea
void TriangleList::makeFacetOnSideOfSea( const std::map<int,int>* const outerBounNode2MeshNode, std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList ){

	const BoundaryCurveList* const ptrBoundaryCurveList = &m_boundaryCurveList;

	const AnalysisDomain* const ptrAnalysisDomain = AnalysisDomain::getInstance();

	const bool makeSurfMesOnAllBoundary = (Control::getInstance())->getFlagWhetherMakeSurfMesOnAllBoundary();

	const int numOuterBoundary = ptrBoundaryCurveList->getTotalNumberOuterBoundaries();
	for( int iBounOuter = 0; iBounOuter < numOuterBoundary; ++iBounOuter ){//==========================================================

		if( ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getGeologicalType() != CommonParameters::SEA &&
			( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getGeologicalType() != CommonParameters::LAKE ){
			continue;// Treat only the sea or lakes
		}

		//std::cout << "iBounOuter : " << iBounOuter << std::endl;

		TriangleList::PolygonNodes nodesSide;
		std::vector<PolygonNodes> nodeSidesArray;

		bool isIntersectPre = false;
		const int numNodes = ptrBoundaryCurveList->getNumNodeOuterBoundary( iBounOuter );
		for( int iNode = 0; iNode < numNodes; ++iNode ){//-------------------------------------------------------------------------------

			const CommonParameters::XY coordXY = ptrBoundaryCurveList->getPointCoordOuterBoundary(iBounOuter, iNode);

			if( ptrAnalysisDomain->doesIntersectWithBoundary( coordXY ) ){

				const int nodeIDBoundCurveList = ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getNodeID(iNode);
				const int nodeID = convertNodeIDBoundCurve2Triangle( nodeIDBoundCurveList );
				
				nodesSide.push_back( nodeID );

				TriangleList::CornerPoint cornerPointType = TriangleList::UNDEFINED_CORNER;
				if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, AnalysisDomain::XPLUS_YMINUS ) ){
					cornerPointType = TriangleList::XPLUS_YMINUS_SURFACE;
				}
				else if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, AnalysisDomain::XPLUS_YPLUS ) ){
					cornerPointType = TriangleList::XPLUS_YPLUS_SURFACE;
				}
				else if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, AnalysisDomain::XMINUS_YPLUS ) ){
					cornerPointType = TriangleList::XMINUS_YPLUS_SURFACE;
				}
				else if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, AnalysisDomain::XMINUS_YMINUS ) ){
					cornerPointType = TriangleList::XMINUS_YMINUS_SURFACE;
				}

				if( cornerPointType != TriangleList::UNDEFINED_CORNER ){
					if( makeSurfMesOnAllBoundary ){
						const TriangleList::PolygonNodes nodesSideTmp = nodesSide;
						for( std::vector<int>::const_reverse_iterator itrNodeSide = nodesSideTmp.rbegin(); itrNodeSide != nodesSideTmp.rend(); ++itrNodeSide ){
							const int nodeIDBoun = convertNodeIDTriangle2BoundCurve(*itrNodeSide);
							std::map<int,int>::const_iterator itrOuterBounNode2MeshNode = outerBounNode2MeshNode[iBounOuter].find(nodeIDBoun);
							if( itrOuterBounNode2MeshNode == outerBounNode2MeshNode[iBounOuter].end() ){// Not found
								OutputFiles::m_logFile << " Error : Node " << *itrNodeSide << " cannot be found in outerBounNode2MeshNode[" << iBounOuter << "] ." << std::endl;
								exit(1);
							}
							if( m_additionalCornerNodes[cornerPointType] < 0 ){
								m_additionalCornerNodes[cornerPointType] = itrOuterBounNode2MeshNode->second;
							}
							if( std::find( nodesSide.begin(), nodesSide.end(), itrOuterBounNode2MeshNode->second ) != nodesSide.end() ){// Already inserted
								break;
							}
							nodesSide.push_back(itrOuterBounNode2MeshNode->second);
						}
						if( !nodeSidesArray.empty() ){
							const CommonParameters::Boundary boundaryType = searchTypeOfSideBoundary(nodeList, nodesSide);
							TriangleList::PolygonNodes nodesSideReversed = nodesSide;
							std::reverse( nodesSideReversed.begin(), nodesSideReversed.end() );
							switch( boundaryType ){
								case CommonParameters::YZ_PLUS:
									makeFacetFromSurroundingNodes( nodesSide, true, CommonParameters::YZ_PLUS, CommonParameters::SEA, nodeList, facetList );
									break;
								case CommonParameters::ZX_PLUS:
									makeFacetFromSurroundingNodes( nodesSide, true, CommonParameters::ZX_PLUS, CommonParameters::SEA, nodeList, facetList );
									break;
								case CommonParameters::YZ_MINUS:
									makeFacetFromSurroundingNodes( nodesSide, false, CommonParameters::YZ_MINUS, CommonParameters::SEA, nodeList, facetList );
									break;
								case CommonParameters::ZX_MINUS:
									makeFacetFromSurroundingNodes( nodesSide, false, CommonParameters::ZX_MINUS, CommonParameters::SEA, nodeList, facetList );
									break;
								default:
									OutputFiles::m_logFile << " Error : Corner point type is wrong : " << boundaryType << std::endl;
									exit(1);
									break;
							}
						}
					}
					else{
						if( m_additionalCornerNodes[cornerPointType] < 0 ){
							const CommonParameters::XYZ cornerCoord = { coordXY.X, coordXY.Y, 0.0 };
							nodeList.push_back( cornerCoord );
							m_additionalCornerNodes[cornerPointType] = static_cast<int>( nodeList.size() ) - 1;
						}
						nodesSide.push_back( m_additionalCornerNodes[cornerPointType] );
						if( !nodeSidesArray.empty() ){
							TriangleList::Facet facetTmp;
							facetTmp.polygons.push_back( nodesSide );
							facetList.push_back( facetTmp );
						}
					}
					nodeSidesArray.push_back( nodesSide );
					nodesSide.clear();
					nodesSide.push_back( m_additionalCornerNodes[cornerPointType] );
					nodesSide.push_back( nodeID );
				}

				isIntersectPre= true; 

			}
			else{

				if( isIntersectPre ){
					if(makeSurfMesOnAllBoundary ){
						const TriangleList::PolygonNodes nodesSideTmp = nodesSide;
						std::vector<int>::const_reverse_iterator itrNodeSide = nodesSideTmp.rbegin();
						++itrNodeSide;
						for( ;itrNodeSide != nodesSideTmp.rend(); ++itrNodeSide ){
							//std::map<int,int>::const_iterator itrOuterBounNode2MeshNode = outerBounNode2MeshNode[iBounOuter].find(*itrNodeSide);
							const int nodeIDBoun = convertNodeIDTriangle2BoundCurve(*itrNodeSide);
							std::map<int,int>::const_iterator itrOuterBounNode2MeshNode = outerBounNode2MeshNode[iBounOuter].find(nodeIDBoun);
							if( itrOuterBounNode2MeshNode == outerBounNode2MeshNode[iBounOuter].end() ){// Not found
								OutputFiles::m_logFile << " Error : Node " << *itrNodeSide << " cannot be found in outerBounNode2MeshNode[" << iBounOuter << "] ." << std::endl;
								exit(1);
							}
							if( std::find( nodesSide.begin(), nodesSide.end(), itrOuterBounNode2MeshNode->second ) != nodesSide.end() ){// Already inserted
								break;
							}
							nodesSide.push_back(itrOuterBounNode2MeshNode->second);
						}
						if( !nodeSidesArray.empty() ){
							const CommonParameters::Boundary boundaryType = searchTypeOfSideBoundary(nodeList, nodesSide);
							switch( boundaryType ){
								case CommonParameters::YZ_PLUS:
									makeFacetFromSurroundingNodes( nodesSide, true, CommonParameters::YZ_PLUS, CommonParameters::SEA, nodeList, facetList );
									break;
								case CommonParameters::ZX_PLUS:
									makeFacetFromSurroundingNodes( nodesSide, true, CommonParameters::ZX_PLUS, CommonParameters::SEA, nodeList, facetList );
									break;
								case CommonParameters::YZ_MINUS:
									makeFacetFromSurroundingNodes( nodesSide, false, CommonParameters::YZ_MINUS, CommonParameters::SEA, nodeList, facetList );
									break;
								case CommonParameters::ZX_MINUS:
									makeFacetFromSurroundingNodes( nodesSide, false, CommonParameters::ZX_MINUS, CommonParameters::SEA, nodeList, facetList );
									break;
								default:
									OutputFiles::m_logFile << " Error : Corner point type is wrong : " << boundaryType << std::endl;
									exit(1);
									break;
							}
						}
					}
					else{
						if( !nodeSidesArray.empty() ){
							TriangleList::Facet facetTmp;
							facetTmp.polygons.push_back( nodesSide );
							facetList.push_back( facetTmp );
						}
					}
					nodeSidesArray.push_back( nodesSide );
					nodesSide.clear();
				}

				isIntersectPre = false;
			}

		}//-------------------------------------------------------------------------------------------------------------

		if( isIntersectPre ){//-------------------------------------------------------------------------------------------------------------
			if(makeSurfMesOnAllBoundary){
				TriangleList::PolygonNodes nodesSideTmp = nodesSide;
				TriangleList::PolygonNodes nodesAdded;
				for( std::vector<int>::const_reverse_iterator itrNodeSide = nodesSideTmp.rbegin() ;itrNodeSide != nodesSideTmp.rend(); ++itrNodeSide ){
					const int nodeIDBoun = convertNodeIDTriangle2BoundCurve(*itrNodeSide);
					std::map<int,int>::const_iterator itrOuterBounNode2MeshNode = outerBounNode2MeshNode[iBounOuter].find(nodeIDBoun);
					if( itrOuterBounNode2MeshNode == outerBounNode2MeshNode[iBounOuter].end() ){// Not found
						OutputFiles::m_logFile << " Error : Node " << *itrNodeSide << " cannot be found in outerBounNode2MeshNode[" << iBounOuter << "] ." << std::endl;
						exit(1);
					}
					if( std::find( nodesSide.begin(), nodesSide.end(), itrOuterBounNode2MeshNode->second ) != nodesSide.end() ){// Already inserted
						break;
					}
					nodesAdded.push_back(itrOuterBounNode2MeshNode->second);
				}
				nodesSide.insert( nodesSide.begin(), nodesAdded.begin(), nodesAdded.end() );
			}

			for( std::vector<int>::const_iterator itrNodeSide = nodeSidesArray[0].begin(); itrNodeSide != nodeSidesArray[0].end(); ++itrNodeSide ){
				if( std::find( nodesSide.begin(), nodesSide.end(), *itrNodeSide ) != nodesSide.end() ){// Already inserted
					break;
				}
				nodesSide.push_back(*itrNodeSide);
			}
			if(makeSurfMesOnAllBoundary){
				const CommonParameters::Boundary boundaryType = searchTypeOfSideBoundary(nodeList, nodesSide);
				switch( boundaryType ){
					case CommonParameters::YZ_PLUS:
						makeFacetFromSurroundingNodes( nodesSide, true, CommonParameters::YZ_PLUS, CommonParameters::SEA, nodeList, facetList );
						break;
					case CommonParameters::ZX_PLUS:
						makeFacetFromSurroundingNodes( nodesSide, true, CommonParameters::ZX_PLUS, CommonParameters::SEA, nodeList, facetList );
						break;
					case CommonParameters::YZ_MINUS:
						makeFacetFromSurroundingNodes( nodesSide, false, CommonParameters::YZ_MINUS, CommonParameters::SEA, nodeList, facetList );
						break;
					case CommonParameters::ZX_MINUS:
						makeFacetFromSurroundingNodes( nodesSide, false, CommonParameters::ZX_MINUS, CommonParameters::SEA, nodeList, facetList );
						break;
					default:
						OutputFiles::m_logFile << " Error : Corner point type is wrong : " << boundaryType << std::endl;
						exit(1);
						break;
				}
			}
			else{
				TriangleList::Facet facetTmp;
				facetTmp.polygons.push_back( nodesSide );
				facetList.push_back( facetTmp );
			}
		}//-------------------------------------------------------------------------------------------------------------

	}//=================================================================================================================

}

// Make boundary curve of the sides of the sea
void TriangleList::makeBoundaryCurveOfSideOfSea( std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList ){

	const BoundaryCurveList* const ptrBoundaryCurveList = &m_boundaryCurveList;

	const AnalysisDomain* const ptrAnalysisDomain = AnalysisDomain::getInstance();

	const bool makeSurfMesOnAllBoundary = (Control::getInstance())->getFlagWhetherMakeSurfMesOnAllBoundary();

	const int numOuterBoundary = ptrBoundaryCurveList->getTotalNumberOuterBoundaries();
	for( int iBounOuter = 0; iBounOuter < numOuterBoundary; ++iBounOuter ){

		if( ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getGeologicalType() != CommonParameters::SEA &&
			( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getGeologicalType() != CommonParameters::LAKE ){
			continue;// Treat only the sea or lakes
		}

		//std::cout << "iBounOuter : " << iBounOuter << std::endl;

		TriangleList::PolygonNodes nodesSide;

		bool isIntersectPre = false;
		const int listSizeInit = static_cast<int>( facetList.size() );
		const int numNodes = ptrBoundaryCurveList->getNumNodeOuterBoundary( iBounOuter );
		for( int iNode = 0; iNode < numNodes; ++iNode ){

			//std::cout << "iNode : " << iNode << std::endl;
	
			const int nodeIDBoundCurveList = ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getNodeID(iNode);	

			//std::cout << "nodeIDBoundCurveList : " << nodeIDBoundCurveList << std::endl;

			const CommonParameters::XY coordXY = ptrBoundaryCurveList->getPointCoordOuterBoundary(iBounOuter, iNode);

			const bool isIntersect = ptrAnalysisDomain->doesIntersectWithBoundary( coordXY ); 

			if( isIntersect ){

				const int nodeID = convertNodeIDBoundCurve2Triangle( nodeIDBoundCurveList );

				//std::cout << "nodeID : " << nodeID << std::endl;
				
				nodesSide.push_back( nodeID );

				int cornerNode = -1;
				if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, AnalysisDomain::XPLUS_YMINUS ) ){

					//std::cout << "XPLUS_YMINUS iNode : " << iNode << std::endl;

					if( m_additionalCornerNodes[TriangleList::XPLUS_YMINUS_SURFACE] < 0 ){
						const CommonParameters::XYZ cornerCoord = { coordXY.X, coordXY.Y, 0.0 };
						nodeList.push_back( cornerCoord );
						m_additionalCornerNodes[TriangleList::XPLUS_YMINUS_SURFACE] = static_cast<int>( nodeList.size() ) - 1;
					}
					nodesSide.push_back( m_additionalCornerNodes[TriangleList::XPLUS_YMINUS_SURFACE] );
					cornerNode = m_additionalCornerNodes[TriangleList::XPLUS_YMINUS_SURFACE];
				}
				else if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, AnalysisDomain::XPLUS_YPLUS ) ){

					//std::cout << "XPLUS_YPLUS iNode : " << iNode << std::endl;

					if( m_additionalCornerNodes[TriangleList::XPLUS_YPLUS_SURFACE] < 0 ){
						const CommonParameters::XYZ cornerCoord = { coordXY.X, coordXY.Y, 0.0 };
						nodeList.push_back( cornerCoord );
						m_additionalCornerNodes[TriangleList::XPLUS_YPLUS_SURFACE] = static_cast<int>( nodeList.size() ) - 1;
					}
					nodesSide.push_back( m_additionalCornerNodes[TriangleList::XPLUS_YPLUS_SURFACE] );
					cornerNode = m_additionalCornerNodes[TriangleList::XPLUS_YPLUS_SURFACE];
				}
				else if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, AnalysisDomain::XMINUS_YPLUS ) ){

					//std::cout << "XMINUS_YPLUS iNode : " << iNode << std::endl;

					if( m_additionalCornerNodes[TriangleList::XMINUS_YPLUS_SURFACE] < 0 ){
						const CommonParameters::XYZ cornerCoord = { coordXY.X, coordXY.Y, 0.0 };
						nodeList.push_back( cornerCoord );
						m_additionalCornerNodes[TriangleList::XMINUS_YPLUS_SURFACE] = static_cast<int>( nodeList.size() ) - 1;
					}
					nodesSide.push_back( m_additionalCornerNodes[TriangleList::XMINUS_YPLUS_SURFACE] );
					cornerNode = m_additionalCornerNodes[TriangleList::XMINUS_YPLUS_SURFACE];
				}
				else if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, AnalysisDomain::XMINUS_YMINUS ) ){

					//std::cout << "XMINUS_YMINUS iNode : " << iNode << std::endl;

					if( m_additionalCornerNodes[TriangleList::XMINUS_YMINUS_SURFACE] < 0 ){
						const CommonParameters::XYZ cornerCoord = { coordXY.X, coordXY.Y, 0.0 };
						nodeList.push_back( cornerCoord );
						m_additionalCornerNodes[TriangleList::XMINUS_YMINUS_SURFACE] = static_cast<int>( nodeList.size() ) - 1;
					}
					nodesSide.push_back( m_additionalCornerNodes[TriangleList::XMINUS_YMINUS_SURFACE] );
					cornerNode = m_additionalCornerNodes[TriangleList::XMINUS_YMINUS_SURFACE];
				}

				if( cornerNode >= 0 ){
					TriangleList::Facet facetTmp;
					facetTmp.polygons.push_back( nodesSide );
					facetList.push_back( facetTmp );
					nodesSide.clear();
					nodesSide.push_back( cornerNode );
					nodesSide.push_back( nodeID );
				}

				isIntersectPre= true; 

			}
			else{

				if( isIntersectPre ){
					TriangleList::Facet facetTmp;
					facetTmp.polygons.push_back( nodesSide );
					facetList.push_back( facetTmp );
					nodesSide.clear();
				}

				isIntersectPre = false;
			}

		}

		if( listSizeInit == static_cast<int>( facetList.size() ) ){
			OutputFiles::m_logFile << " Error : Size of list was not increased at the function " << __FUNCTION__ << " ." << std::endl;
			exit(1);
		}

		if( isIntersectPre ){
			//facetList[listSizeInit - 1].polygons.insert( facetList[listSizeInit].polygons.begin(), nodesSide.begin(), nodesSide.end() );
			(facetList[listSizeInit].polygons.begin())->insert( (facetList[listSizeInit].polygons.begin())->begin(), nodesSide.begin(), nodesSide.end() );
		}

	}

}

// Make PLCs for the surface of the sea
void TriangleList::makeFacetOnSeaSurface( std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList ){

	const BoundaryCurveList* const ptrBoundaryCurveList = &m_boundaryCurveList;
	const AnalysisDomain* const ptrAnalysisDomain = AnalysisDomain::getInstance();

	const int numOuterBoundary = ptrBoundaryCurveList->getTotalNumberOuterBoundaries();
	for( int iBounOuter = 0; iBounOuter < numOuterBoundary; ++iBounOuter ){

		if( ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getGeologicalType() != CommonParameters::SEA &&
			( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getGeologicalType() != CommonParameters::LAKE ){
			continue;// Treat only the sea or lakes
		}

		TriangleList::PolygonNodes nodesSide;

		//bool isIntersectPre = false;
		const int numNodes = ptrBoundaryCurveList->getNumNodeOuterBoundary( iBounOuter );

		const CommonParameters::XY coordXYLast = ptrBoundaryCurveList->getPointCoordOuterBoundary(iBounOuter, numNodes-1);
		bool isIntersectPre = ptrAnalysisDomain->doesIntersectWithBoundary( coordXYLast ); 
		const int listSizeInit = static_cast<int>( facetList.size() );

		for( int iNode = 0; iNode < numNodes; ++iNode ){//-------------------------------------------------------------------------

			const int nodeIDBoundCurveList = ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getNodeID(iNode);	
			const int nodeID = convertNodeIDBoundCurve2Triangle( nodeIDBoundCurveList );

			const CommonParameters::XY coordXY = ptrBoundaryCurveList->getPointCoordOuterBoundary(iBounOuter, iNode);
			const bool isIntersect = ptrAnalysisDomain->doesIntersectWithBoundary( coordXY ); 

			if( isIntersect ){

				int cornerNode = -1;
				if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, AnalysisDomain::XPLUS_YMINUS ) ){
					if( m_additionalCornerNodes[TriangleList::XPLUS_YMINUS_SURFACE] < 0 ){
						OutputFiles::m_logFile << " Error : Additional corner node has not been set yet." << std::endl;
						exit(1);
					}
					nodesSide.push_back( m_additionalCornerNodes[TriangleList::XPLUS_YMINUS_SURFACE] );
					cornerNode = m_additionalCornerNodes[TriangleList::XPLUS_YMINUS_SURFACE];
				}
				else if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, AnalysisDomain::XPLUS_YPLUS ) ){
					if( m_additionalCornerNodes[TriangleList::XPLUS_YPLUS_SURFACE] < 0 ){
						OutputFiles::m_logFile << " Error : Additional corner node has not been set yet." << std::endl;
						exit(1);
					}
					nodesSide.push_back( m_additionalCornerNodes[TriangleList::XPLUS_YPLUS_SURFACE] );
					cornerNode = m_additionalCornerNodes[TriangleList::XPLUS_YPLUS_SURFACE];
				}
				else if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, AnalysisDomain::XMINUS_YPLUS ) ){
					if( m_additionalCornerNodes[TriangleList::XMINUS_YPLUS_SURFACE] < 0 ){
						OutputFiles::m_logFile << " Error : Additional corner node has not been set yet." << std::endl;
						exit(1);
					}
					nodesSide.push_back( m_additionalCornerNodes[TriangleList::XMINUS_YPLUS_SURFACE] );
					cornerNode = m_additionalCornerNodes[TriangleList::XMINUS_YPLUS_SURFACE];
				}
				else if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, AnalysisDomain::XMINUS_YMINUS ) ){
					if( m_additionalCornerNodes[TriangleList::XMINUS_YMINUS_SURFACE] < 0 ){
						OutputFiles::m_logFile << " Error : Additional corner node has not been set yet." << std::endl;
						exit(1);
					}
					nodesSide.push_back( m_additionalCornerNodes[TriangleList::XMINUS_YMINUS_SURFACE] );
					cornerNode = m_additionalCornerNodes[TriangleList::XMINUS_YMINUS_SURFACE];
				}
				else{
					const int iNodeNext = (iNode + 1) % numNodes;
					const CommonParameters::XY coordXYNext = ptrBoundaryCurveList->getPointCoordOuterBoundary(iBounOuter, iNodeNext);
					if( !isIntersectPre || !ptrAnalysisDomain->doesIntersectWithBoundary( coordXYNext ) ){
						// Next point does not intersect with the outer edge of the analysis domain
						nodesSide.push_back( nodeID );
					}
				}

				isIntersectPre= true; 

			}
			else{
				nodesSide.push_back( nodeID );

				isIntersectPre = false;
			}

		}//----------------------------------------------------------------------------------------------------------------------------

		TriangleList::Facet facetTmp;
		facetTmp.polygons.push_back( nodesSide );

		////----------------------------------------------------------------------------------------------------------------------------
		////----- Search node pairs located on the boundary curve except the nodes of the segments belonging to the boundary curve -----
		////----------------------------------------------------------------------------------------------------------------------------
		//std::vector< std::pair<int,int> > nodePair;
		//for( int iNode1 = 0; iNode1 < numNodes; ++iNode1 ){
		//	for( int i2 = 2; i2 < numNodes - 1; ++i2 ){// Except the nodes of the segments belonging to the boundary curve

		//		const int iNode2 = ( iNode1 + i2 ) % numNodes;

		//		const int nodeIDBoundCurve1 = ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getNodeID(iNode1);
		//		const int nodeIDBoundCurve2 = ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getNodeID(iNode2);
		//		const int nodeID1 = convertNodeIDBoundCurve2Triangle(nodeIDBoundCurve1);
		//		const int nodeID2 = convertNodeIDBoundCurve2Triangle(nodeIDBoundCurve2);

		//		if( shareSameTriangle( nodeID1, nodeID2 ) >= 0 ){
		//			const int nd1 = nodeID1 < nodeID2 ? nodeID1 : nodeID2;
		//			const int nd2 = nodeID1 < nodeID2 ? nodeID2 : nodeID1;
		//			nodePair.push_back( std::make_pair( nd1, nd2 ) );
		//		}

		//	}
		//}
		//std::sort(nodePair.begin(), nodePair.end());
		//nodePair.erase( std::unique( nodePair.begin(), nodePair.end() ), nodePair.end() );

		//for( std::vector< std::pair<int,int> >::iterator itrNodePair = nodePair.begin(); itrNodePair != nodePair.end(); ++itrNodePair ){

		//	std::cout << "Node Pair : " <<  itrNodePair->first << " " << itrNodePair->second << std::endl;

		//	TriangleList::PolygonNodes nodesSide;
		//	nodesSide.push_back( itrNodePair->first );
		//	nodesSide.push_back( itrNodePair->second );
		//	TriangleList::Facet facetTmp;
		//	facetTmp.polygons.push_back( nodesSide );
		//}

		//--------------------------------
		//----- For inner boundaries -----
		//--------------------------------
		const int numBounInnerIncluded = ptrBoundaryCurveList->getNumInnerBoundaryIncluded( iBounOuter );
		for( int iBounInner = 0; iBounInner < numBounInnerIncluded; ++iBounInner ){

			if( ( ptrBoundaryCurveList->getPointerToInnerBoundary( iBounInner ) )->getGeologicalType() != CommonParameters::LAND ){
				continue;// Treat only inner surface of land
			}

			TriangleList::PolygonNodes nodesInnerBound;
			const int numNodesInnerBound = ptrBoundaryCurveList->getNumNodeInnerBoundary( iBounInner );
			for( int iNodeInnerBound = 0; iNodeInnerBound < numNodesInnerBound; ++iNodeInnerBound ){

				const int nodeIDInnerBound = ( ptrBoundaryCurveList->getPointerToInnerBoundary(iBounInner) )->getNodeID( iNodeInnerBound );
				nodesInnerBound.push_back( convertNodeIDBoundCurve2Triangle(nodeIDInnerBound) );

			}

			facetTmp.polygons.push_back( nodesInnerBound );

			// Calculate coordinate of hole
			//bool found(false);
			std::vector<int> triangleIDs = ( m_nodeList.getPointerToNode( nodesInnerBound[0] ) )->getTriangleIDs();
			for( std::vector<int>::iterator itr = triangleIDs.begin(); itr != triangleIDs.end(); ++itr ){
				if( m_triangles[*itr].getDomainType() == CommonParameters::LAND ){
					const CommonParameters::XY holeCoordXY = m_triangles[*itr].getCoordGravCenter( &m_nodeList );
					const CommonParameters::XYZ holeCoord = { holeCoordXY.X, holeCoordXY.Y, 0.0 };
					facetTmp.holes.push_back( holeCoord );
					break;
				}
			}

			if( facetTmp.holes.empty() ){
				OutputFiles::m_logFile << " Error : Could't find trialgle of land surface at the function " << __FUNCTION__ << " ." << std::endl;
				exit(1);
			}

		}

		facetList.push_back( facetTmp );

		if( listSizeInit == static_cast<int>( facetList.size() ) ){
			OutputFiles::m_logFile << " Error : Size of list was not increased at the function " << __FUNCTION__ << " ." << std::endl;
			exit(1);
		}
		
	}

	//----- Do not delete for future use >>>>>
	//std::map<int,int> nodeIDEarth2SeaSurf;

	//const AnalysisDomain* const ptrAnalysisDomain = AnalysisDomain::getInstance();

	//// Surface of the sea
	//const int numTriangles = static_cast<int>( m_triangles.size() );
	//for( std::vector< Triangle >::const_iterator itr = m_triangles.begin(); itr != m_triangles.end(); ++itr ){

	//	if( itr->getDomainType() != CommonParameters::SEA && itr->getDomainType() != CommonParameters::LAKE ){
	//		continue;
	//	}

	//	//int nodeIDs[3] = { -1, -1, -1 };

	//	TriangleList::PolygonNodes nodes;

	//	for( int iNode = 0; iNode < 3; ++iNode ){

	//		const int nodeIDOrg = itr->getNodeID(iNode);
	//		int nodeID(-1);

	//		const CommonParameters::XY coordXY = m_nodeList.getCoordXYOfPoints( nodeIDOrg );

	//		if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, AnalysisDomain::XPLUS_YMINUS ) ){
	//			nodeID = additionalCornerNodes[TriangleList::XPLUS_YMINUS_SURFACE];
	//		}
	//		else if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, AnalysisDomain::XPLUS_YPLUS ) ){
	//			nodeID = additionalCornerNodes[TriangleList::XPLUS_YPLUS_SURFACE];
	//		}
	//		else if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, AnalysisDomain::XMINUS_YPLUS ) ){
	//			nodeID = additionalCornerNodes[TriangleList::XMINUS_YPLUS_SURFACE];
	//		}
	//		else if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, AnalysisDomain::XMINUS_YMINUS ) ){
	//			nodeID = additionalCornerNodes[TriangleList::XMINUS_YMINUS_SURFACE];
	//		}
	//		else{// No locate on a corner of the domain

	//			if( m_nodeList.getPointerToNode( nodeIDOrg )->getLocation() != Node::COAST_LINE ){
	//	
	//				std::map<int,int>::iterator itrMap = nodeIDEarth2SeaSurf.find( nodeIDOrg );

	//				if( itrMap == nodeIDEarth2SeaSurf.end() ){// Not found
	//					CommonParameters::XYZ coord = {coordXY.X, coordXY.Y, 0.0 };
	//					nodeList.push_back( coord );
	//					nodeID = static_cast<int>( nodeList.size() ) - 1;
	//					nodeIDEarth2SeaSurf.insert( std::make_pair( nodeIDOrg, nodeID ) );
	//				}
	//				else{
	//					nodeID = itrMap->second;
	//				}

	//			}
	//			else{

	//				nodeID = nodeIDOrg;

	//			}

	//		}

	//		nodes.push_back( nodeID );

	//	}

	//	TriangleList::Facet facetTmp;
	//	facetTmp.polygons.push_back( nodes );
	//	facetList.push_back( facetTmp );

	//}
	//----- Do not delete for future use <<<<<

}

// Make PLCs for the sides of the air and land
void TriangleList::makeFacetOnSideOfAirAndLand( const std::map<int,int>* const outerBounNode2MeshNode, const CommonParameters::Boundary& boundaryType, std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList ){

	const BoundaryCurveList* const ptrBoundaryCurveList = &m_boundaryCurveList;

	CommonParameters::XYZ cornerCoordsOfLandSurface[4] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
	calcCornerCoordsOfLandSurface( cornerCoordsOfLandSurface );

	const bool includeSedimentLayer = (Control::getInstance())->getIncludeSedimentLayer();
	if( includeSedimentLayer ){
		const double thickness = (Control::getInstance())->getThicknessSediment();
		for( int iCorner = 0; iCorner < 4; ++iCorner ){
			cornerCoordsOfLandSurface[iCorner].Z += thickness;
		}
	}

	const AnalysisDomain* const ptrAnalysisDomain = AnalysisDomain::getInstance();

	const double xMin = ptrAnalysisDomain->getMinCoordX();
	const double xMax = ptrAnalysisDomain->getMaxCoordX(); 
	const double yMin = ptrAnalysisDomain->getMinCoordY(); 
	const double yMax = ptrAnalysisDomain->getMaxCoordY(); 

	//const CommonParameters::XY coordXPlusYMinus  = { xMax, yMin };
	//const CommonParameters::XY coordXPlusYPlus   = { xMax, yMax };
	//const CommonParameters::XY coordXMinusYPlus  = { xMin, yMax };
	//const CommonParameters::XY coordXMinusYMinus = { xMin, yMin };

	CommonParameters::XY targerCoord = { -1.0, -1.0 };

	TriangleList::CornerPoint startConerTop = TriangleList::UNDEFINED_CORNER;
	TriangleList::CornerPoint startConerSurface = TriangleList::UNDEFINED_CORNER;
	TriangleList::CornerPoint startConerBottom = TriangleList::UNDEFINED_CORNER;
	TriangleList::Side startSideAir = UNDEFINED_SIDE;
	TriangleList::Side startSideLand = UNDEFINED_SIDE;
	TriangleList::CornerPoint endConerTop = TriangleList::UNDEFINED_CORNER;
	TriangleList::CornerPoint endConerSurface = TriangleList::UNDEFINED_CORNER;
	TriangleList::CornerPoint endConerBottom = TriangleList::UNDEFINED_CORNER;
	TriangleList::Side endSideAir = UNDEFINED_SIDE;
	TriangleList::Side endSideLand = UNDEFINED_SIDE;
	AnalysisDomain::FourCorner startCorner = AnalysisDomain::UNDEFINED_CORNER;
	AnalysisDomain::FourCorner endCorner = AnalysisDomain::UNDEFINED_CORNER;
	bool reverseAir = false;
	bool reverseLand = false;

	switch( boundaryType ){
		case CommonParameters::YZ_MINUS:
			startConerTop =     TriangleList::XMINUS_YPLUS_TOP;
			startConerSurface = TriangleList::XMINUS_YPLUS_SURFACE;
			startConerBottom =  TriangleList::XMINUS_YPLUS_BOTTOM;
			startSideAir = XMINUS_YPLUS_VERTICAL_AIR;
			startSideLand = XMINUS_YPLUS_VERTICAL_LAND;
			startCorner = AnalysisDomain::XMINUS_YPLUS;
			targerCoord.X = xMin;
			targerCoord.Y = yMax;

			endConerTop =     TriangleList::XMINUS_YMINUS_TOP;
			endConerSurface = TriangleList::XMINUS_YMINUS_SURFACE;
			endConerBottom =  TriangleList::XMINUS_YMINUS_BOTTOM;
			endSideAir = XMINUS_YMINUS_VERTICAL_AIR;
			endSideLand = XMINUS_YMINUS_VERTICAL_LAND;
			endCorner = AnalysisDomain::XMINUS_YMINUS;

			reverseAir = false;
			reverseLand = true;
			break;
		case CommonParameters::YZ_PLUS:
			startConerTop =     TriangleList::XPLUS_YMINUS_TOP;
			startConerSurface = TriangleList::XPLUS_YMINUS_SURFACE;
			startConerBottom =  TriangleList::XPLUS_YMINUS_BOTTOM;
			startSideAir = XPLUS_YMINUS_VERTICAL_AIR;
			startSideLand = XPLUS_YMINUS_VERTICAL_LAND;
			startCorner = AnalysisDomain::XPLUS_YMINUS;
			targerCoord.X = xMax;
			targerCoord.Y = yMin;

			endConerTop =     TriangleList::XPLUS_YPLUS_TOP;
			endConerSurface = TriangleList::XPLUS_YPLUS_SURFACE;
			endConerBottom =  TriangleList::XPLUS_YPLUS_BOTTOM;
			endSideAir = XPLUS_YPLUS_VERTICAL_AIR;
			endSideLand = XPLUS_YPLUS_VERTICAL_LAND;
			endCorner = AnalysisDomain::XPLUS_YPLUS;

			reverseAir = true;
			reverseLand = false;
			break;
		case CommonParameters::ZX_MINUS:
			startConerTop =     TriangleList::XMINUS_YMINUS_TOP;
			startConerSurface = TriangleList::XMINUS_YMINUS_SURFACE;
			startConerBottom =  TriangleList::XMINUS_YMINUS_BOTTOM;
			startSideAir = XMINUS_YMINUS_VERTICAL_AIR;
			startSideLand = XMINUS_YMINUS_VERTICAL_LAND;
			startCorner = AnalysisDomain::XMINUS_YMINUS;
			targerCoord.X = xMin;
			targerCoord.Y = yMin;

			endConerTop =     TriangleList::XPLUS_YMINUS_TOP;
			endConerSurface = TriangleList::XPLUS_YMINUS_SURFACE;
			endConerBottom =  TriangleList::XPLUS_YMINUS_BOTTOM;
			endSideAir = XPLUS_YMINUS_VERTICAL_AIR;
			endSideLand = XPLUS_YMINUS_VERTICAL_LAND;
			endCorner = AnalysisDomain::XPLUS_YMINUS;

			reverseAir = false;
			reverseLand = true;
			break;
		case CommonParameters::ZX_PLUS:
			startConerTop =     TriangleList::XPLUS_YPLUS_TOP;
			startConerSurface = TriangleList::XPLUS_YPLUS_SURFACE;
			startConerBottom =  TriangleList::XPLUS_YPLUS_BOTTOM;
			startSideAir = XPLUS_YPLUS_VERTICAL_AIR;
			startSideLand = XPLUS_YPLUS_VERTICAL_LAND;
			startCorner = AnalysisDomain::XPLUS_YPLUS;
			targerCoord.X = xMax;
			targerCoord.Y = yMax;

			endConerTop =     TriangleList::XMINUS_YPLUS_TOP;
			endConerSurface = TriangleList::XMINUS_YPLUS_SURFACE;
			endConerBottom =  TriangleList::XMINUS_YPLUS_BOTTOM;
			endSideAir = XMINUS_YPLUS_VERTICAL_AIR;
			endSideLand = XMINUS_YPLUS_VERTICAL_LAND;
			endCorner = AnalysisDomain::XMINUS_YPLUS;

			reverseAir = true;
			reverseLand = false;
			break;
		default:
			OutputFiles::m_logFile << " Error : Wrong boundary type : " << boundaryType << " ." << std::endl;
			exit(1);
			break;
	};

	const bool makeSurfMesOnAllBoundary = (Control::getInstance())->getFlagWhetherMakeSurfMesOnAllBoundary();
	int nodeIDCur = static_cast<int>( nodeList.size() );

	TriangleList::PolygonNodes nodesSideLand;
	TriangleList::PolygonNodes nodesSideAir;

	nodesSideLand.push_back( m_additionalCornerNodes[startConerBottom] );
	nodesSideAir.push_back( m_additionalCornerNodes[startConerTop] );

	if(makeSurfMesOnAllBoundary){
		if( m_additionalSideNodes[startSideLand].empty() ){
			const CommonParameters::XYZ coordBot = nodeList[m_additionalCornerNodes[startConerBottom]];
			CommonParameters::XYZ coordSurf = cornerCoordsOfLandSurface[startCorner];
			std::vector<double> terms;
			Util::calculateTermsOfGeometricProgressionExceptLast(
				(Control::getInstance())->calcMaximumEdgeLength( Util::transformXYZToXY(coordSurf), boundaryType ),
				(Control::getInstance())->calcMaximumEdgeLength( Util::transformXYZToXY(coordBot), boundaryType ),
				coordBot.Z - coordSurf.Z, terms );
			for( std::vector<double>::const_iterator itr = terms.begin(); itr != terms.end(); ++itr ){
				const CommonParameters::XYZ coord = { coordSurf.X, coordSurf.Y, coordSurf.Z + *itr };
				nodeList.push_back( coord );
				m_additionalSideNodes[startSideLand].push_back(nodeIDCur++);
			}
		}
		nodesSideLand.insert( nodesSideLand.end(), m_additionalSideNodes[startSideLand].rbegin(), m_additionalSideNodes[startSideLand].rend() );

		if( m_additionalSideNodes[startSideAir].empty() ){
			const CommonParameters::XYZ coordTop = nodeList[m_additionalCornerNodes[startConerTop]];
			CommonParameters::XYZ coordSurf = { coordTop.X, coordTop.Y, 0.0 };
			std::vector<double> terms;
			Util::calculateTermsOfGeometricProgressionExceptLast(
				(Control::getInstance())->calcMaximumEdgeLength( Util::transformXYZToXY(coordTop), boundaryType ),
				(Control::getInstance())->calcMaximumEdgeLength( Util::transformXYZToXY(coordSurf), boundaryType ),
				coordSurf.Z - coordTop.Z, terms );
			for( std::vector<double>::const_iterator itr = terms.begin(); itr != terms.end(); ++itr ){
				const CommonParameters::XYZ coord = { coordTop.X, coordTop.Y, coordTop.Z + *itr };
				nodeList.push_back( coord );
				m_additionalSideNodes[startSideAir].push_back(nodeIDCur++);
			}
		}
		nodesSideAir.insert( nodesSideAir.end(), m_additionalSideNodes[startSideAir].begin(), m_additionalSideNodes[startSideAir].end() );
	}

	int courveIDPre(-1);

	//------------------------------------------------------------------ LOOP >>>>>>>>>>
	while(true){

		int curveID(-1);
		int nodeID(-1);

		if( courveIDPre < 0 ){// The first iteration
			ptrBoundaryCurveList->searchPointOfOuterBoundaryByCoord( targerCoord, curveID, nodeID );
		}
		else{
			ptrBoundaryCurveList->searchPointOfOuterBoundaryByCoordExceptACurve( targerCoord, courveIDPre, curveID, nodeID );
		}

		const int nodeIDBoundCurve = ( ptrBoundaryCurveList->getPointerToOuterBoundary(curveID) )->getNodeID( nodeID );

		bool finished = false;
		if( ( ptrBoundaryCurveList->getPointerToOuterBoundary( curveID ) )->getGeologicalType() == CommonParameters::LAND ){// boundary curve surrounds land

			nodesSideLand.push_back( convertNodeIDBoundCurve2Triangle(nodeIDBoundCurve) );
			nodesSideAir.push_back( convertNodeIDBoundCurve2Triangle(nodeIDBoundCurve) );
	
			CommonParameters::XY coordXYPre = ptrBoundaryCurveList->getPointCoordOuterBoundary(curveID, nodeID);
			const int numNodes = ptrBoundaryCurveList->getNumNodeOuterBoundary( curveID );
			for( int i = 1; i < numNodes; ++i ){
				const int iNode = ( nodeID + i ) % numNodes;

				const CommonParameters::XY coordXY = ptrBoundaryCurveList->getPointCoordOuterBoundary(curveID, iNode);

				if( !ptrAnalysisDomain->doesIntersectWithBoundary( coordXY ) ){
					nodesSideLand.pop_back();
					nodesSideAir.pop_back();
					targerCoord = coordXYPre;
					courveIDPre = curveID;
					break;
				}

				const int nodeIDBoundCurve = ( ptrBoundaryCurveList->getPointerToOuterBoundary(curveID) )->getNodeID( iNode );
				nodesSideLand.push_back( convertNodeIDBoundCurve2Triangle(nodeIDBoundCurve) );
				nodesSideAir.push_back( convertNodeIDBoundCurve2Triangle(nodeIDBoundCurve) );

				if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, endCorner ) ){
					if(makeSurfMesOnAllBoundary){
						if( m_additionalSideNodes[endSideLand].empty() ){
							const CommonParameters::XYZ coordBot = nodeList[m_additionalCornerNodes[endConerBottom]];
							CommonParameters::XYZ coordSurf = cornerCoordsOfLandSurface[endCorner];
							std::vector<double> terms;
							Util::calculateTermsOfGeometricProgressionExceptLast(
								(Control::getInstance())->calcMaximumEdgeLength( Util::transformXYZToXY(coordSurf), boundaryType ),
								(Control::getInstance())->calcMaximumEdgeLength( Util::transformXYZToXY(coordBot), boundaryType ),
								coordBot.Z - coordSurf.Z, terms );
							for( std::vector<double>::const_iterator itr = terms.begin(); itr != terms.end(); ++itr ){
								const CommonParameters::XYZ coord = { coordSurf.X, coordSurf.Y, coordSurf.Z + *itr };
								nodeList.push_back( coord );
								m_additionalSideNodes[endSideLand].push_back(nodeIDCur++);
							}
						}
						nodesSideLand.insert( nodesSideLand.end(), m_additionalSideNodes[endSideLand].begin(), m_additionalSideNodes[endSideLand].end() );

						if( m_additionalSideNodes[endSideAir].empty() ){
							const CommonParameters::XYZ coordTop = nodeList[m_additionalCornerNodes[endConerTop]];
							CommonParameters::XYZ coordSurf = { coordTop.X, coordTop.Y, 0.0 };
							std::vector<double> terms;
							Util::calculateTermsOfGeometricProgressionExceptLast(
								(Control::getInstance())->calcMaximumEdgeLength( Util::transformXYZToXY(coordTop), boundaryType ),
								(Control::getInstance())->calcMaximumEdgeLength( Util::transformXYZToXY(coordSurf), boundaryType ),
								coordSurf.Z - coordTop.Z, terms );
							for( std::vector<double>::const_iterator itr = terms.begin(); itr != terms.end(); ++itr ){
								const CommonParameters::XYZ coord = { coordTop.X, coordTop.Y, coordTop.Z + *itr };
								nodeList.push_back( coord );
								m_additionalSideNodes[endSideAir].push_back(nodeIDCur++);
							}
						}
						nodesSideAir.insert( nodesSideAir.end(), m_additionalSideNodes[endSideAir].rbegin(), m_additionalSideNodes[endSideAir].rend() );
					}
					nodesSideLand.push_back( m_additionalCornerNodes[endConerBottom] );
					nodesSideAir.push_back( m_additionalCornerNodes[endConerTop] );
					finished = true;
					break;
				}

				coordXYPre = coordXY;
			}

		}
		else{// boundary curve surrounds the sea and lakes

			const int nodeIDTriangle = convertNodeIDBoundCurve2Triangle(nodeIDBoundCurve);
			if( includeSedimentLayer ){
				if( courveIDPre >= 0 ){
					nodesSideLand.push_back(nodeIDTriangle);
				} 
				std::map<int,int>::const_iterator itrNodeLandSurfToSediment = m_nodeLandSurfToSediment.find(nodeIDTriangle);
				if( itrNodeLandSurfToSediment == m_nodeLandSurfToSediment.end() ){// Not found
					OutputFiles::m_logFile << " Error : Node " << nodeIDTriangle << " cannot be found in m_nodeLandSurfToSediment." << std::endl;
					exit(1);
				}
				nodesSideLand.push_back( itrNodeLandSurfToSediment->second );
			}
			else{
				nodesSideLand.push_back(nodeIDTriangle);
			}

			if(makeSurfMesOnAllBoundary){
				std::map<int,int>::const_iterator itrOuterBounNode2MeshNode = outerBounNode2MeshNode[curveID].find(nodeIDBoundCurve);
				if( itrOuterBounNode2MeshNode == outerBounNode2MeshNode[curveID].end() ){// Not found
					OutputFiles::m_logFile << " Error : Node " << nodeIDBoundCurve << " cannot be found in outerBounNode2MeshNode[" << curveID << "] ." << std::endl;
					exit(1);
				}
				nodesSideAir.push_back(itrOuterBounNode2MeshNode->second);
			}
			else{
				if( courveIDPre < 0 ){
					nodesSideAir.push_back( m_additionalCornerNodes[startConerSurface] );
				}
				else{
					nodesSideAir.push_back( convertNodeIDBoundCurve2Triangle(nodeIDBoundCurve) );
				}
			}
		
			CommonParameters::XY coordXYPre = ptrBoundaryCurveList->getPointCoordOuterBoundary(curveID, nodeID); 
			const int numNodes = ptrBoundaryCurveList->getNumNodeOuterBoundary( curveID );
			for( int i = 1; i < numNodes; ++i ){
				const int iNode = ( nodeID + i ) % numNodes;

				const CommonParameters::XY coordXY = ptrBoundaryCurveList->getPointCoordOuterBoundary(curveID, iNode);

				if( !ptrAnalysisDomain->doesIntersectWithBoundary( coordXY ) ){
					if( !includeSedimentLayer ){					
						nodesSideLand.pop_back();
					}
					if(makeSurfMesOnAllBoundary){
						nodesSideAir.pop_back();
					}
					targerCoord = coordXYPre;
					courveIDPre = curveID;
					break;
				}

				const int nodeIDBoundCurve = ( ptrBoundaryCurveList->getPointerToOuterBoundary(curveID) )->getNodeID( iNode );
				const int nodeIDTriangle = convertNodeIDBoundCurve2Triangle(nodeIDBoundCurve);
				if( includeSedimentLayer ){
					std::map<int,int>::const_iterator itrNodeLandSurfToSediment = m_nodeLandSurfToSediment.find(nodeIDTriangle);
					if( itrNodeLandSurfToSediment == m_nodeLandSurfToSediment.end() ){// Not found
						OutputFiles::m_logFile << " Error : Node " << nodeIDTriangle << " cannot be found in m_nodeLandSurfToSediment." << std::endl;
						exit(1);
					}
					nodesSideLand.push_back(itrNodeLandSurfToSediment->second);
				}
				else{
					nodesSideLand.push_back(nodeIDTriangle);
				}
				if(makeSurfMesOnAllBoundary){
					std::map<int,int>::const_iterator itrOuterBounNode2MeshNode = outerBounNode2MeshNode[curveID].find(nodeIDBoundCurve);
					if( itrOuterBounNode2MeshNode == outerBounNode2MeshNode[curveID].end() ){// Not found
						OutputFiles::m_logFile << " Error : Node " << nodeIDBoundCurve << " cannot be found in outerBounNode2MeshNode[" << curveID << "] ." << std::endl;
						exit(1);
					}
					nodesSideAir.push_back(itrOuterBounNode2MeshNode->second);
				}

				if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, endCorner ) ){

					if(makeSurfMesOnAllBoundary){
						if( m_additionalSideNodes[endSideLand].empty() ){
							const CommonParameters::XYZ coordBot = nodeList[m_additionalCornerNodes[endConerBottom]];
							CommonParameters::XYZ coordSurf = cornerCoordsOfLandSurface[endCorner];
							std::vector<double> terms;
							Util::calculateTermsOfGeometricProgressionExceptLast(
								(Control::getInstance())->calcMaximumEdgeLength( Util::transformXYZToXY(coordSurf), boundaryType ),
								(Control::getInstance())->calcMaximumEdgeLength( Util::transformXYZToXY(coordBot), boundaryType ),
								coordBot.Z - coordSurf.Z, terms );
							for( std::vector<double>::const_iterator itr = terms.begin(); itr != terms.end(); ++itr ){
								const CommonParameters::XYZ coord = { coordSurf.X, coordSurf.Y, coordSurf.Z + *itr };
								nodeList.push_back( coord );
								m_additionalSideNodes[endSideLand].push_back(nodeIDCur++);
							}
						}
						nodesSideLand.insert( nodesSideLand.end(), m_additionalSideNodes[endSideLand].begin(), m_additionalSideNodes[endSideLand].end() );

						if( m_additionalSideNodes[endSideAir].empty() ){
							const CommonParameters::XYZ coordTop = nodeList[m_additionalCornerNodes[endConerTop]];
							CommonParameters::XYZ coordSurf = { coordTop.X, coordTop.Y, 0.0 };
							std::vector<double> terms;
							Util::calculateTermsOfGeometricProgressionExceptLast(
								(Control::getInstance())->calcMaximumEdgeLength( Util::transformXYZToXY(coordTop), boundaryType ),
								(Control::getInstance())->calcMaximumEdgeLength( Util::transformXYZToXY(coordSurf), boundaryType ),
								coordSurf.Z - coordTop.Z, terms );
							for( std::vector<double>::const_iterator itr = terms.begin(); itr != terms.end(); ++itr ){
								const CommonParameters::XYZ coord = { coordTop.X, coordTop.Y, coordTop.Z + *itr };
								nodeList.push_back( coord );
								m_additionalSideNodes[endSideAir].push_back(nodeIDCur++);
							}
						}
						nodesSideAir.insert( nodesSideAir.end(), m_additionalSideNodes[endSideAir].rbegin(), m_additionalSideNodes[endSideAir].rend() );
					}
					else{
						nodesSideAir.push_back( m_additionalCornerNodes[endConerSurface] );
					}

					nodesSideLand.push_back( m_additionalCornerNodes[endConerBottom] );
					nodesSideAir.push_back( m_additionalCornerNodes[endConerTop] );

					finished = true;
					break;
				}

				coordXYPre = coordXY;
			}

		}

		if( finished ){

			//std::cout << "break" << std::endl;

			if(makeSurfMesOnAllBoundary){
				switch( boundaryType ){
					case CommonParameters::YZ_MINUS:
						nodesSideLand.insert( nodesSideLand.end(), m_additionalSideNodes[XMINUS_BOTTOM].begin(), m_additionalSideNodes[XMINUS_BOTTOM].end() );
						nodesSideAir.insert( nodesSideAir.end(), m_additionalSideNodes[XMINUS_TOP].begin(), m_additionalSideNodes[XMINUS_TOP].end() );
						break;
					case CommonParameters::YZ_PLUS:
						nodesSideLand.insert( nodesSideLand.end(), m_additionalSideNodes[XPLUS_BOTTOM].rbegin(), m_additionalSideNodes[XPLUS_BOTTOM].rend() );
						nodesSideAir.insert( nodesSideAir.end(), m_additionalSideNodes[XPLUS_TOP].rbegin(), m_additionalSideNodes[XPLUS_TOP].rend() );
						break;
					case CommonParameters::ZX_MINUS:
						nodesSideLand.insert( nodesSideLand.end(), m_additionalSideNodes[YMINUS_BOTTOM].rbegin(), m_additionalSideNodes[YMINUS_BOTTOM].rend() );
						nodesSideAir.insert( nodesSideAir.end(), m_additionalSideNodes[YMINUS_TOP].rbegin(), m_additionalSideNodes[YMINUS_TOP].rend() );
						break;
					case CommonParameters::ZX_PLUS:
						nodesSideLand.insert( nodesSideLand.end(), m_additionalSideNodes[YPLUS_BOTTOM].begin(), m_additionalSideNodes[YPLUS_BOTTOM].end() );
						nodesSideAir.insert( nodesSideAir.end(), m_additionalSideNodes[YPLUS_TOP].begin(), m_additionalSideNodes[YPLUS_TOP].end() );
						break;
					default:
						OutputFiles::m_logFile << " Error : Wrong boundary type : " << boundaryType << " ." << std::endl;
						exit(1);
						break;
				};
				makeFacetFromSurroundingNodes( nodesSideAir, reverseAir, boundaryType, CommonParameters::AIR, nodeList, facetList );
				makeFacetFromSurroundingNodes( nodesSideLand, reverseLand, boundaryType, CommonParameters::LAND, nodeList, facetList );
			}
			else{
				TriangleList::Facet facetTmpLand;
				facetTmpLand.polygons.push_back( nodesSideLand );
				facetList.push_back( facetTmpLand );

				TriangleList::Facet facetTmpAir;
				facetTmpAir.polygons.push_back( nodesSideAir );
				facetList.push_back( facetTmpAir );
			}
			break;
		}

	}//------------------------------------------------------------------ LOOP <<<<<<<<<
	
}

// Make PLCs for lower layers
void TriangleList::makeFacetOfLowerLayers( std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList ){

	if( !(Control::getInstance())->getIncludeLayers() ){
		return;
	}

	const int numLayers = (Control::getInstance())->getNumLayers();

	// Side type
	const int static XMINUS = 0;
	const int static YMINUS = 1;
	const int static XPLUS = 2;
	const int static YPLUS = 3;

	// Corner point type
	const int static XPLUS_YMINUS = 0;
	const int static XPLUS_YPLUS = 1;
	const int static XMINUS_YPLUS = 2;
	const int static XMINUS_YMINUS = 3;

	const AnalysisDomain* const ptrAnalysisDomain = AnalysisDomain::getInstance();
	const double xMin = ptrAnalysisDomain->getMinCoordX();
	const double xMax = ptrAnalysisDomain->getMaxCoordX(); 
	const double yMin = ptrAnalysisDomain->getMinCoordY(); 
	const double yMax = ptrAnalysisDomain->getMaxCoordY(); 

	TriangleList::PolygonNodes nodesSide[4];
	TriangleList::PolygonNodes nodesBottom;

	int additionalCornerNodesPre[4] = { -1, -1, -1, -1 };
	additionalCornerNodesPre[XPLUS_YMINUS]  = m_additionalCornerNodes[TriangleList::XPLUS_YMINUS_BOTTOM];
	additionalCornerNodesPre[XPLUS_YPLUS]   = m_additionalCornerNodes[TriangleList::XPLUS_YPLUS_BOTTOM];
	additionalCornerNodesPre[XMINUS_YPLUS]  = m_additionalCornerNodes[TriangleList::XMINUS_YPLUS_BOTTOM];
	additionalCornerNodesPre[XMINUS_YMINUS] = m_additionalCornerNodes[TriangleList::XMINUS_YMINUS_BOTTOM];

	std::vector<int> additionalHorizontalNodesPre[4];
	additionalHorizontalNodesPre[YMINUS] = m_additionalSideNodes[TriangleList::YMINUS_BOTTOM];
	additionalHorizontalNodesPre[XPLUS]  = m_additionalSideNodes[TriangleList::XPLUS_BOTTOM];
	additionalHorizontalNodesPre[YPLUS]  = m_additionalSideNodes[TriangleList::YPLUS_BOTTOM];
	additionalHorizontalNodesPre[XMINUS] = m_additionalSideNodes[TriangleList::XMINUS_BOTTOM];

	CommonParameters::XYZ cornerCoordPre[4] = { -1, -1, -1, -1 };
	for( int i = 0; i < 4; ++i ){
		cornerCoordPre[i] = nodeList[additionalCornerNodesPre[i]];
	}

	for( int iz = 0; iz < numLayers; ++iz ){


		const int layerIndex = iz + 1;
		const bool isBottom = layerIndex == numLayers;

		const double coordZ = isBottom ? ptrAnalysisDomain->getMaxCoordZ() : (Control::getInstance())->getDepthLayerInterfaces(iz+1);

		CommonParameters::XYZ cornerCoord[4] = { { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 }, { -1, -1, -1 } };
		int additionalCornerNodes[4] = { -1, -1, -1, -1 };

		int nodeIDCur = static_cast<int>( nodeList.size() );
		for( int i = 0; i < 4; ++i ){
			const CommonParameters::XYZ coord = { cornerCoordPre[i].X, cornerCoordPre[i].Y, coordZ };
			nodeList.push_back( coord );
			cornerCoord[i] = coord;
			additionalCornerNodes[i]  = nodeIDCur++; 
		}

		std::vector<int> additionalHorizontalNodes[4];
		{// YZ Plane (X Minus)
			std::vector<double> terms;
			Util::calculateTermsOfGeometricProgressionExceptLast(
				(Control::getInstance())->calcMaximumEdgeLength( cornerCoord[XMINUS_YMINUS] ),
				(Control::getInstance())->calcMaximumEdgeLength( cornerCoord[XMINUS_YPLUS] ),
				yMax - yMin, terms );
			for( std::vector<double>::const_iterator itr = terms.begin(); itr != terms.end(); ++itr ){
				const CommonParameters::XYZ coord = { xMin, yMin + *itr, coordZ };
				nodeList.push_back( coord );
				additionalHorizontalNodes[XMINUS].push_back(nodeIDCur++);
			}
		}

		{// YZ Plane (X Plus)
			std::vector<double> terms;
			Util::calculateTermsOfGeometricProgressionExceptLast(
				(Control::getInstance())->calcMaximumEdgeLength( cornerCoord[XPLUS_YMINUS] ),
				(Control::getInstance())->calcMaximumEdgeLength( cornerCoord[XPLUS_YPLUS] ),
				yMax - yMin, terms );
			for( std::vector<double>::const_iterator itr = terms.begin(); itr != terms.end(); ++itr ){
				const CommonParameters::XYZ coord = { xMax, yMin + *itr, coordZ };
				nodeList.push_back( coord );
				additionalHorizontalNodes[XPLUS].push_back(nodeIDCur++);
			}
		}

		{// ZX Plane (Y Minus)
			std::vector<double> terms;
			Util::calculateTermsOfGeometricProgressionExceptLast(
				(Control::getInstance())->calcMaximumEdgeLength( cornerCoord[XMINUS_YMINUS] ),
				(Control::getInstance())->calcMaximumEdgeLength( cornerCoord[XPLUS_YMINUS] ),
				xMax - xMin, terms );
			for( std::vector<double>::const_iterator itr = terms.begin(); itr != terms.end(); ++itr ){
				const CommonParameters::XYZ coord = { xMin + *itr, yMin, coordZ };
				nodeList.push_back( coord );
				additionalHorizontalNodes[YMINUS].push_back(nodeIDCur++);
			}
		}

		{// ZX Plane (Y Plus)
			std::vector<double> terms;
			Util::calculateTermsOfGeometricProgressionExceptLast(
				(Control::getInstance())->calcMaximumEdgeLength( cornerCoord[XMINUS_YPLUS] ),
				(Control::getInstance())->calcMaximumEdgeLength( cornerCoord[XPLUS_YPLUS] ),
				xMax - xMin, terms );
			for( std::vector<double>::const_iterator itr = terms.begin(); itr != terms.end(); ++itr ){
				const CommonParameters::XYZ coord = { xMin + *itr, yMax, coordZ };
				nodeList.push_back( coord );
				additionalHorizontalNodes[YPLUS].push_back(nodeIDCur++);
			}
		}

		TriangleList::PolygonNodes nodesBot;
		nodesBot.push_back(additionalCornerNodes[XMINUS_YMINUS]);
		for( std::vector<int>::const_iterator itr = additionalHorizontalNodes[YMINUS].begin(); itr != additionalHorizontalNodes[YMINUS].end(); ++itr ){
			nodesBot.push_back(*itr);
		}
		nodesBot.push_back(additionalCornerNodes[XPLUS_YMINUS]);
		for( std::vector<int>::const_iterator itr = additionalHorizontalNodes[XPLUS].begin(); itr != additionalHorizontalNodes[XPLUS].end(); ++itr ){
			nodesBot.push_back(*itr);
		}
		nodesBot.push_back(additionalCornerNodes[XPLUS_YPLUS]);
		for( std::vector<int>::const_reverse_iterator itr = additionalHorizontalNodes[YPLUS].rbegin(); itr != additionalHorizontalNodes[YPLUS].rend(); ++itr ){
			nodesBot.push_back(*itr);
		}
		nodesBot.push_back(additionalCornerNodes[XMINUS_YPLUS]);
		for( std::vector<int>::const_reverse_iterator itr = additionalHorizontalNodes[XMINUS].rbegin(); itr != additionalHorizontalNodes[XMINUS].rend(); ++itr ){
			nodesBot.push_back(*itr);
		}

		std::vector<int> additionalVerticalNodes[4];
		for( int i = 0; i < 4; ++i ){
			const CommonParameters::XYZ coordBotPre  = nodeList[additionalCornerNodesPre[i]];
			const CommonParameters::XYZ coordBot = nodeList[additionalCornerNodes[i]];
			std::vector<double> terms;
			Util::calculateTermsOfGeometricProgressionExceptLast(
				(Control::getInstance())->calcMaximumEdgeLength( coordBotPre ),
				(Control::getInstance())->calcMaximumEdgeLength( coordBot ),
				coordBot.Z - coordBotPre.Z, terms );
			for( std::vector<double>::const_iterator itr = terms.begin(); itr != terms.end(); ++itr ){
				const CommonParameters::XYZ coord = { coordBotPre.X, coordBotPre.Y, coordBotPre.Z + *itr };
				nodeList.push_back( coord );
				additionalVerticalNodes[i].push_back(nodeIDCur++);
			}
		}
	
		if(isBottom){
			makeFacetFromSurroundingNodes( nodesBot, false, CommonParameters::BOT, CommonParameters::LAND, nodeList, facetList );
		}else{
			makeFacetFromSurroundingNodes( nodesBot, false, CommonParameters::LAYER, CommonParameters::LAND, nodeList, facetList, layerIndex );
		}
		
		{// YZ Plane (X Minus)
			TriangleList::PolygonNodes nodesSide;
			nodesSide.push_back( additionalCornerNodesPre[XMINUS_YMINUS] );
			nodesSide.insert( nodesSide.end(), additionalVerticalNodes[XMINUS_YMINUS].begin(), additionalVerticalNodes[XMINUS_YMINUS].end() );
			nodesSide.push_back( additionalCornerNodes[XMINUS_YMINUS] );
			nodesSide.insert( nodesSide.end(), additionalHorizontalNodes[XMINUS].begin(), additionalHorizontalNodes[XMINUS].end() );
			nodesSide.push_back( additionalCornerNodes[XMINUS_YPLUS] );
			nodesSide.insert( nodesSide.end(), additionalVerticalNodes[XMINUS_YPLUS].rbegin(), additionalVerticalNodes[XMINUS_YPLUS].rend() );
			nodesSide.push_back( additionalCornerNodesPre[XMINUS_YPLUS] );
			nodesSide.insert( nodesSide.end(), additionalHorizontalNodesPre[XMINUS].rbegin(), additionalHorizontalNodesPre[XMINUS].rend() );
			makeFacetFromSurroundingNodes( nodesSide, true, CommonParameters::YZ_MINUS, CommonParameters::LAND, nodeList, facetList );
		}

		{// YZ Plane (X Plus)
			TriangleList::PolygonNodes nodesSide;
			nodesSide.push_back( additionalCornerNodesPre[XPLUS_YMINUS] );
			nodesSide.insert( nodesSide.end(), additionalVerticalNodes[XPLUS_YMINUS].begin(), additionalVerticalNodes[XPLUS_YMINUS].end() );
			nodesSide.push_back( additionalCornerNodes[XPLUS_YMINUS] );
			nodesSide.insert( nodesSide.end(), additionalHorizontalNodes[XPLUS].begin(), additionalHorizontalNodes[XPLUS].end() );
			nodesSide.push_back( additionalCornerNodes[XPLUS_YPLUS] );
			nodesSide.insert( nodesSide.end(), additionalVerticalNodes[XPLUS_YPLUS].rbegin(), additionalVerticalNodes[XPLUS_YPLUS].rend() );
			nodesSide.push_back( additionalCornerNodesPre[XPLUS_YPLUS] );
			nodesSide.insert( nodesSide.end(), additionalHorizontalNodesPre[XPLUS].rbegin(), additionalHorizontalNodesPre[XPLUS].rend() );
			makeFacetFromSurroundingNodes( nodesSide, true, CommonParameters::YZ_PLUS, CommonParameters::LAND, nodeList, facetList );
		}

		{// ZX Plane (Y Minus)
			TriangleList::PolygonNodes nodesSide;
			nodesSide.push_back( additionalCornerNodesPre[XMINUS_YMINUS] );
			nodesSide.insert( nodesSide.end(), additionalVerticalNodes[XMINUS_YMINUS].begin(), additionalVerticalNodes[XMINUS_YMINUS].end() );
			nodesSide.push_back( additionalCornerNodes[XMINUS_YMINUS] );
			nodesSide.insert( nodesSide.end(), additionalHorizontalNodes[YMINUS].begin(), additionalHorizontalNodes[YMINUS].end() );
			nodesSide.push_back( additionalCornerNodes[XPLUS_YMINUS] );
			nodesSide.insert( nodesSide.end(), additionalVerticalNodes[XPLUS_YMINUS].rbegin(), additionalVerticalNodes[XPLUS_YMINUS].rend() );
			nodesSide.push_back( additionalCornerNodesPre[XPLUS_YMINUS] );
			nodesSide.insert( nodesSide.end(), additionalHorizontalNodesPre[YMINUS].rbegin(), additionalHorizontalNodesPre[YMINUS].rend() );
			makeFacetFromSurroundingNodes( nodesSide, false, CommonParameters::ZX_MINUS, CommonParameters::LAND, nodeList, facetList );
		}

		{// ZX Plane (Y Plus)
			TriangleList::PolygonNodes nodesSide;
			nodesSide.push_back( additionalCornerNodesPre[XMINUS_YPLUS] );
			nodesSide.insert( nodesSide.end(), additionalVerticalNodes[XMINUS_YPLUS].begin(), additionalVerticalNodes[XMINUS_YPLUS].end() );
			nodesSide.push_back( additionalCornerNodes[XMINUS_YPLUS] );
			nodesSide.insert( nodesSide.end(), additionalHorizontalNodes[YPLUS].begin(), additionalHorizontalNodes[YPLUS].end() );
			nodesSide.push_back( additionalCornerNodes[XPLUS_YPLUS] );
			nodesSide.insert( nodesSide.end(), additionalVerticalNodes[XPLUS_YPLUS].rbegin(), additionalVerticalNodes[XPLUS_YPLUS].rend() );
			nodesSide.push_back( additionalCornerNodesPre[XPLUS_YPLUS] );
			nodesSide.insert( nodesSide.end(), additionalHorizontalNodesPre[YPLUS].rbegin(), additionalHorizontalNodesPre[YPLUS].rend() );
			makeFacetFromSurroundingNodes( nodesSide, false, CommonParameters::ZX_PLUS, CommonParameters::LAND, nodeList, facetList );
		}

		// Cur => Pre
		for( int i = 0; i < 4; ++i ){
			additionalCornerNodesPre[i] = additionalCornerNodes[i];
			additionalHorizontalNodesPre[i] = additionalHorizontalNodes[i];
			cornerCoordPre[i] = cornerCoord[i];
		}
	}

}

// Make PLCs for extended region
void TriangleList::makeFacetOfExtendedRegion( std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList ){

	if( !( Control::getInstance() )->isExtendedRegionAdded() ){
		return;
	}

	//const BoundaryCurveList* const ptrBoundaryCurveList = BoundaryCurveList::getInstance();
	//const int numOuterBoundary = ptrBoundaryCurveList->getTotalNumberOuterBoundaries();
	//for( int iBounOuter = 0; iBounOuter < numOuterBoundary; ++iBounOuter ){
	//	if( ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getGeologicalType() == CommonParameters::LAND ){
	//		OutputFiles::m_logFile << " Error : Outermost boundary curves cannot include land if extended region was added. " << std::endl;
	//		exit(1);
	//	}
	//}

	const AnalysisDomain* const ptrAnalysisDomain = AnalysisDomain::getInstance();

	//const double xMin = ptrAnalysisDomain->getMinCoordX();
	//const double xMax = ptrAnalysisDomain->getMaxCoordX(); 
	//const double yMin = ptrAnalysisDomain->getMinCoordY(); 
	//const double yMax = ptrAnalysisDomain->getMaxCoordY(); 
	//const double zMin = ptrAnalysisDomain->getMinCoordZ(); 
	//const double zMax = ptrAnalysisDomain->getMaxCoordZ(); 

	//const double lengPlusX = ( Control::getInstance() )->getLengthOfExtendedRegionPlusX();
	//const double lengMinusX = ( Control::getInstance() )->getLengthOfExtendedRegionMinusX();
	//const double lengPlusY = ( Control::getInstance() )->getLengthOfExtendedRegionPlusY();
	//const double lengMinusY = ( Control::getInstance() )->getLengthOfExtendedRegionMinusY();

	const double xMinExtended = ptrAnalysisDomain->getMinCoordX() - ( Control::getInstance() )->getLengthOfExtendedRegionMinusX();
	const double xMaxExtended = ptrAnalysisDomain->getMaxCoordX() + ( Control::getInstance() )->getLengthOfExtendedRegionPlusX(); 
	const double yMinExtended = ptrAnalysisDomain->getMinCoordY() - ( Control::getInstance() )->getLengthOfExtendedRegionMinusY(); 
	const double yMaxExtended = ptrAnalysisDomain->getMaxCoordY() + ( Control::getInstance() )->getLengthOfExtendedRegionPlusY();
	const double zMin = ptrAnalysisDomain->getMinCoordZ(); 
	const double zMax = ptrAnalysisDomain->getMaxCoordZ(); 
	const double xMid = 0.5*(ptrAnalysisDomain->getMinCoordX()+ptrAnalysisDomain->getMaxCoordX());
	const double yMid = 0.5*(ptrAnalysisDomain->getMinCoordY()+ptrAnalysisDomain->getMaxCoordY());

	int nodeIDLast = static_cast<int>( nodeList.size() ) - 1;

	CommonParameters::XYZ coord = { 0.0, 0.0, 0.0 };
	coord.X = xMaxExtended;
	coord.Y = yMinExtended;
	coord.Z = zMin;
	nodeList.push_back( coord );
	const int nodeIDXPlusYMinusTopExtended = ++nodeIDLast;

	coord.X = xMaxExtended;
	coord.Y = yMaxExtended;
	nodeList.push_back( coord );
	const int nodeIDXPlusYPlusTopExtended = ++nodeIDLast;

	coord.X = xMinExtended;
	coord.Y = yMaxExtended;
	nodeList.push_back( coord );
	const int nodeIDXMinusYPlusTopExtended = ++nodeIDLast;

	coord.X = xMinExtended;
	coord.Y = yMinExtended;
	nodeList.push_back( coord );
	const int nodeIDXMinusYMinusTopExtended  = ++nodeIDLast;

	coord.X = xMaxExtended;
	coord.Y = yMinExtended;
	coord.Z = 0.0;
	nodeList.push_back( coord );
	const int nodeIDXPlusYMinusSurfExtended = ++nodeIDLast;

	coord.X = xMaxExtended;
	coord.Y = yMaxExtended;
	nodeList.push_back( coord );
	const int nodeIDXPlusYPlusSurfExtended = ++nodeIDLast;

	coord.X = xMinExtended;
	coord.Y = yMaxExtended;
	nodeList.push_back( coord );
	const int nodeIDXMinusYPlusSurfExtended = ++nodeIDLast;

	coord.X = xMinExtended;
	coord.Y = yMinExtended;
	nodeList.push_back( coord );
	const int nodeIDXMinusYMinusSurfExtended = ++nodeIDLast;

	coord.X = xMaxExtended;
	coord.Y = yMinExtended;
	coord.Z = zMax;
	nodeList.push_back( coord );
	const int nodeIDXPlusYMinusBotExtended = ++nodeIDLast;

	coord.X = xMaxExtended;
	coord.Y = yMaxExtended;
	nodeList.push_back( coord );
	const int nodeIDXPlusYPlusBotExtended = ++nodeIDLast;

	coord.X = xMinExtended;
	coord.Y = yMaxExtended;
	nodeList.push_back( coord );
	const int nodeIDXMinusYPlusBotExtended = ++nodeIDLast;

	coord.X = xMinExtended;
	coord.Y = yMinExtended;
	nodeList.push_back( coord );
	const int nodeIDXMinusYMinusBotExtended = ++nodeIDLast;

	{// Top
		TriangleList::PolygonNodes nodesOuter;
		nodesOuter.push_back(nodeIDXPlusYMinusTopExtended);
		nodesOuter.push_back(nodeIDXPlusYPlusTopExtended);
		nodesOuter.push_back(nodeIDXMinusYPlusTopExtended);
		nodesOuter.push_back(nodeIDXMinusYMinusTopExtended);

		TriangleList::PolygonNodes nodesInner;
		nodesInner.push_back(m_additionalCornerNodes[XPLUS_YMINUS_TOP]);
		nodesInner.push_back(m_additionalCornerNodes[XPLUS_YPLUS_TOP]);
		nodesInner.push_back(m_additionalCornerNodes[XMINUS_YPLUS_TOP]);
		nodesInner.push_back(m_additionalCornerNodes[XMINUS_YMINUS_TOP]);

		TriangleList::Facet facetTmp;
		facetTmp.polygons.push_back( nodesOuter );
		facetTmp.polygons.push_back( nodesInner );
		const CommonParameters::XYZ holeCoord = { xMid, yMid, zMin };
		facetTmp.holes.push_back( holeCoord );
		facetList.push_back( facetTmp );
	}
	
	{// Surface
		TriangleList::PolygonNodes nodesOuter;
		nodesOuter.push_back(nodeIDXPlusYMinusSurfExtended);
		nodesOuter.push_back(nodeIDXPlusYPlusSurfExtended);
		nodesOuter.push_back(nodeIDXMinusYPlusSurfExtended);
		nodesOuter.push_back(nodeIDXMinusYMinusSurfExtended);

		//TriangleList::PolygonNodes nodesInner;
		//nodesInner.push_back(additionalCornerNodes[XPLUS_YMINUS_SURFACE]);
		//nodesInner.push_back(additionalCornerNodes[XPLUS_YPLUS_SURFACE]);
		//nodesInner.push_back(additionalCornerNodes[XMINUS_YPLUS_SURFACE]);
		//nodesInner.push_back(additionalCornerNodes[XMINUS_YMINUS_SURFACE]);

		TriangleList::Facet facetTmp;
		facetTmp.polygons.push_back( nodesOuter );
		//facetTmp.polygons.push_back( nodesInner );
		std::vector<int> nodesSurroundigAnalysisDomain;
		getNodesSurroundingAnalysisDomain( m_additionalCornerNodes, nodeList, nodesSurroundigAnalysisDomain );
		facetTmp.polygons.push_back( nodesSurroundigAnalysisDomain );
		const CommonParameters::XYZ holeCoord = { xMid, yMid, 0.0 };
		facetTmp.holes.push_back( holeCoord );
		facetList.push_back( facetTmp );
	}

	{// Bottom
		TriangleList::PolygonNodes nodesOuter;
		nodesOuter.push_back(nodeIDXPlusYMinusBotExtended);
		nodesOuter.push_back(nodeIDXPlusYPlusBotExtended);
		nodesOuter.push_back(nodeIDXMinusYPlusBotExtended);
		nodesOuter.push_back(nodeIDXMinusYMinusBotExtended);

		TriangleList::PolygonNodes nodesInner;
		nodesInner.push_back(m_additionalCornerNodes[XPLUS_YMINUS_BOTTOM]);
		nodesInner.push_back(m_additionalCornerNodes[XPLUS_YPLUS_BOTTOM]);
		nodesInner.push_back(m_additionalCornerNodes[XMINUS_YPLUS_BOTTOM]);
		nodesInner.push_back(m_additionalCornerNodes[XMINUS_YMINUS_BOTTOM]);

		TriangleList::Facet facetTmp;
		facetTmp.polygons.push_back( nodesOuter );
		facetTmp.polygons.push_back( nodesInner );
		const CommonParameters::XYZ holeCoord = { xMid, yMid, zMax };
		facetTmp.holes.push_back( holeCoord );
		facetList.push_back( facetTmp );
	}

	{// Y-Z Plane Plus Side
		TriangleList::PolygonNodes nodes1;
		nodes1.push_back(nodeIDXPlusYMinusTopExtended);
		nodes1.push_back(nodeIDXPlusYPlusTopExtended);
		nodes1.push_back(nodeIDXPlusYPlusSurfExtended);
		nodes1.push_back(nodeIDXPlusYMinusSurfExtended);

		TriangleList::PolygonNodes nodes2;
		nodes2.push_back(nodeIDXPlusYMinusSurfExtended);
		nodes2.push_back(nodeIDXPlusYPlusSurfExtended);
		nodes2.push_back(nodeIDXPlusYPlusBotExtended);
		nodes2.push_back(nodeIDXPlusYMinusBotExtended);

		TriangleList::Facet facetTmp;
		facetTmp.polygons.push_back( nodes1 );
		facetTmp.polygons.push_back( nodes2 );
		facetList.push_back( facetTmp );
	}

	{// Z-X Plane Plus Side
		TriangleList::PolygonNodes nodes1;
		nodes1.push_back(nodeIDXPlusYPlusTopExtended);
		nodes1.push_back(nodeIDXMinusYPlusTopExtended);
		nodes1.push_back(nodeIDXMinusYPlusSurfExtended);
		nodes1.push_back(nodeIDXPlusYPlusSurfExtended);

		TriangleList::PolygonNodes nodes2;
		nodes2.push_back(nodeIDXPlusYPlusSurfExtended);
		nodes2.push_back(nodeIDXMinusYPlusSurfExtended);
		nodes2.push_back(nodeIDXMinusYPlusBotExtended);
		nodes2.push_back(nodeIDXPlusYPlusBotExtended);

		TriangleList::Facet facetTmp;
		facetTmp.polygons.push_back( nodes1 );
		facetTmp.polygons.push_back( nodes2 );
		facetList.push_back( facetTmp );
	}

	{// Y-Z Plane Minus Side
		TriangleList::PolygonNodes nodes1;
		nodes1.push_back(nodeIDXMinusYPlusTopExtended);
		nodes1.push_back(nodeIDXMinusYMinusTopExtended);
		nodes1.push_back(nodeIDXMinusYMinusSurfExtended);
		nodes1.push_back(nodeIDXMinusYPlusSurfExtended);

		TriangleList::PolygonNodes nodes2;
		nodes2.push_back(nodeIDXMinusYPlusSurfExtended);
		nodes2.push_back(nodeIDXMinusYMinusSurfExtended);
		nodes2.push_back(nodeIDXMinusYMinusBotExtended);
		nodes2.push_back(nodeIDXMinusYPlusBotExtended);

		TriangleList::Facet facetTmp;
		facetTmp.polygons.push_back( nodes1 );
		facetTmp.polygons.push_back( nodes2 );
		facetList.push_back( facetTmp );
	}

	{// Z-X Plane Minus Side
		TriangleList::PolygonNodes nodes1;
		nodes1.push_back(nodeIDXMinusYMinusTopExtended);
		nodes1.push_back(nodeIDXPlusYMinusTopExtended);
		nodes1.push_back(nodeIDXPlusYMinusSurfExtended);
		nodes1.push_back(nodeIDXMinusYMinusSurfExtended);

		TriangleList::PolygonNodes nodes2;
		nodes2.push_back(nodeIDXMinusYMinusSurfExtended);
		nodes2.push_back(nodeIDXPlusYMinusSurfExtended);
		nodes2.push_back(nodeIDXPlusYMinusBotExtended);
		nodes2.push_back(nodeIDXMinusYMinusBotExtended);

		TriangleList::Facet facetTmp;
		facetTmp.polygons.push_back( nodes1 );
		facetTmp.polygons.push_back( nodes2 );
		facetList.push_back( facetTmp );
	}

	//const CommonParameters::XYZ coordXPlusYMinusPlusXTop  = {  xMax + lengPlusX,  yMin, zMin };
	//const CommonParameters::XYZ coordXPlusYPlusPlusXTop   = {  xMax + lengPlusX,  yMax,  zMin };
	//const CommonParameters::XYZ coordXMinusYPlusMinusXTop  = { xMin - lengMinusX, yMax,  zMin };
	//const CommonParameters::XYZ coordXMinusYMinusMinusXTop = { xMin - lengMinusX, yMin, zMin };

	//const CommonParameters::XYZ coordXPlusYMinusSurfExtended  = { xMax + lengPlusX,  yMin - lengMinusY, 0.0 };
	//const CommonParameters::XYZ coordXPlusYPlusSurfExtended   = { xMax + lengPlusX,  yMax + lengPlusY,  0.0 };
	//const CommonParameters::XYZ coordXMinusYPlusSurfExtended  = { xMin - lengMinusX, yMax + lengPlusY,  0.0 };
	//const CommonParameters::XYZ coordXMinusYMinusSurfExtended = { xMin - lengMinusX, yMin - lengMinusY, 0.0 };

	//const CommonParameters::XYZ coordXPlusYMinusBotExtended  = { xMax + lengPlusX,  yMin - lengMinusY, zMax };
	//const CommonParameters::XYZ coordXPlusYPlusBotExtended   = { xMax + lengPlusX,  yMax + lengPlusY,  zMax };
	//const CommonParameters::XYZ coordXMinusYPlusBotExtended  = { xMin - lengMinusX, yMax + lengPlusY,  zMax };
	//const CommonParameters::XYZ coordXMinusYMinusBotExtended = { xMin - lengMinusX, yMin - lengMinusY, zMax };
	
}

// Write PLCs to poly file
void TriangleList::writePLCsToPolyFile( const std::string& fileName, const std::vector<CommonParameters::XYZ>& nodeList, const std::vector<TriangleList::Facet>& facetList ){

	// Open output vtk file -----
	std::ofstream ofsVTK( fileName.c_str() );
	if( !ofsVTK ) {
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << std::endl;
		exit(1);
	}

	const int numNodes = static_cast<int>( nodeList.size() );

	ofsVTK << std::setw(10) << numNodes << std::setw(10) << 3  << std::setw(10) << 0 << std::setw(10) << 0 << std::endl;
	ofsVTK << "# Part 1" << std::endl;

	ofsVTK.precision(15);
	ofsVTK << std::fixed;

	int icount(0);
	for( std::vector<CommonParameters::XYZ>::const_iterator itr = nodeList.begin(); itr != nodeList.end(); ++itr ){
		ofsVTK << std::setw(10) << ++icount
			   << std::setw(30) << std::scientific << itr->X
			   << std::setw(30) << std::scientific << itr->Y
			   << std::setw(30) << std::scientific << itr->Z << std::endl;
	}

	const int numFacet = static_cast<int>( facetList.size() );

	ofsVTK << "# Part 2" << std::endl;

	ofsVTK <<  std::setw(10) << numFacet <<  std::setw(10) << 0 << std::endl;

	int icountSurf = 0;
	int icountHole = 0;
	for( std::vector<TriangleList::Facet>::const_iterator itr = facetList.begin(); itr != facetList.end(); ++itr ){

		//ofsVTK << std::setw(10) << itr->polygons.size() << std::setw(10) << itr->holes.size() << std::setw(10) << -1 << std::endl;
		ofsVTK << std::setw(10) << itr->polygons.size() << std::setw(10) << itr->holes.size() << " #" << ++icountSurf << std::endl;
	
		for( std::vector<PolygonNodes>::const_iterator itrPoly = itr->polygons.begin(); itrPoly != itr->polygons.end(); ++itrPoly ){

			ofsVTK << std::setw(10) << static_cast<int>( itrPoly->size() );
			for( std::vector<int>::const_iterator itrNode = itrPoly->begin(); itrNode != itrPoly->end(); ++itrNode ){
				ofsVTK << std::setw(10) << *itrNode + 1;
			}
			ofsVTK << std::endl;
		}

		for( std::vector<CommonParameters::XYZ>::const_iterator itrHoles = itr->holes.begin(); itrHoles != itr->holes.end(); ++itrHoles ){

			ofsVTK << std::setw(10) << ++icountHole
				   << std::setw(30) << std::scientific << itrHoles->X
				   << std::setw(30) << std::scientific << itrHoles->Y
				   << std::setw(30) << std::scientific << itrHoles->Z << std::endl;

		}

	}

	ofsVTK << "# Part 3" << std::endl;
	ofsVTK << 0 << std::endl;

	ofsVTK << "# Part 4" << std::endl;
	ofsVTK << 0 << std::endl;

	ofsVTK.close();


}

// Convert node ID of boundary list to the one of triangle list
int TriangleList::convertNodeIDBoundCurve2Triangle( const int nodeIDBoundCurve ) const{

	std::map<int,int>::const_iterator itr = m_nodeIDBoundCurve2Triangle.find( nodeIDBoundCurve );

	if( itr == m_nodeIDBoundCurve2Triangle.end() ){
		OutputFiles::m_logFile << " Error : Node " << nodeIDBoundCurve << " of the boundary curve cannot be found in the node list of boundary curves ." << std::endl;
		exit(1);
	}
	
	return itr->second;

}

// Convert node ID of triangle to the one of boundary list
int TriangleList::convertNodeIDTriangle2BoundCurve( const int nodeIDTriangle ) const{

	std::map<int,int>::const_iterator itr = m_nodeIDTriangle2BoundCurve.find( nodeIDTriangle );

	//if( itr == m_nodeIDTriangle2BoundCurve.end() ){
	//	OutputFiles::m_logFile << " Error : Node " << nodeIDTriangle << " of the triangles cannot be found in the node list of triangles ." << std::endl;
	//	exit(1);
	//}
	if( itr == m_nodeIDTriangle2BoundCurve.end() ){
		return -1;
	}
	return itr->second;

}

//// Get flag whether the specified two node share the same triangle
//bool TriangleList::shareSameTriangle( const int nodeID1, const int nodeID2 ) const{
//
//	std::cout << "nodeID1 nodeID2 : " << nodeID1 << " " << nodeID2 << std::endl;
//
//	std::vector<int> triangles = m_nodeList.getTriangleIDs( nodeID1 );
//	for( std::vector<int>::iterator itr = triangles.begin(); itr != triangles.end(); ++itr ){
//		for( int iNode = 0; iNode < 3; ++iNode ){
//
//			std::cout << "triangle iNode  m_triangles[*itr].getNodeID(iNode) : " << *itr << " " << iNode << " " << m_triangles[*itr].getNodeID(iNode) << std::endl;
//
//			if(	m_triangles[*itr].getNodeID(iNode) == nodeID2 ){
//				return true;
//			}
//		}
//	}
//
//	return false;
//
//}

// Get ID of the triangle which has the two specified nodes.
// Return -1 if the two specified nodes does not share any triangle.
int TriangleList::shareSameTriangle( const int nodeID1, const int nodeID2 ) const{

	std::vector<int> triangles = m_nodeList.getTriangleIDs( nodeID1 );
	for( std::vector<int>::iterator itr = triangles.begin(); itr != triangles.end(); ++itr ){
		for( int iNode = 0; iNode < 3; ++iNode ){
			if(	m_triangles[*itr].getNodeID(iNode) == nodeID2 ){
				return *itr;
			}
		}
	}

	return -1;

}

// Get ID of the triangle which has the two specified nodes except domain boundary.
// Return -1 if the two specified nodes does not share any triangle or edge connecting the two nodes corresponds to a domain boundary.
int TriangleList::shareSameTriangleExceptDomainBoundary( const int nodeID1, const int nodeID2 ) const{

	int triangleIDs[2];
	int icount(0);
	std::vector<int> triangles = m_nodeList.getTriangleIDs( nodeID1 );
	for( std::vector<int>::iterator itr = triangles.begin(); itr != triangles.end(); ++itr ){
		for( int iNode = 0; iNode < 3; ++iNode ){
			if(	m_triangles[*itr].getNodeID(iNode) == nodeID2 ){
				//return *itr;
				triangleIDs[icount++] = *itr;
				break;
			}
		}
		if( icount >= 2 ){
			break;
		}
	}
	
	if( icount < 2 ){
		return -1;
	}
	else if( m_triangles[triangleIDs[0]].getDomainType() != m_triangles[triangleIDs[1]].getDomainType() ){
		//OutputFiles::m_logFile << " Warning : Edge connecting the two nodes ( " << nodeID1 << ", " << nodeID2 << " ) corresponds to a domain boundary." << std::endl;
		return -1;
	}

	return triangleIDs[0];

}

// Subdivide polygon to keep boundary
void TriangleList::subdividePolygon( std::vector<int>& triangles, const int nodeIDStart , const int nodeIDEnd ){

	if( static_cast<int>( triangles.size() ) == 2 ){// Swap triangles
		flipEdgeOfTwoTriangles( triangles, nodeIDStart, nodeIDEnd );
	}
	else{
		std::vector<int> nodeGroup1;
		std::vector<int> nodeGroup2;
		std::vector< std::pair<int,int> > surroundingElemEdges;
		searchNodesOfPolygon( triangles, nodeIDStart, nodeIDEnd, nodeGroup1, nodeGroup2, surroundingElemEdges );

		// Delete the data relating to the polygon
		deleteDataRelatingToPolygon( triangles );

		// Fill triangles in polygon
		int iElem(0);
		fillTrianglesInPolygon( triangles, nodeGroup1, nodeIDStart, nodeIDEnd, iElem );
		fillTrianglesInPolygon( triangles, nodeGroup2, nodeIDEnd, nodeIDStart, iElem );

		// Reconstruct adjacency relationship
		reconstructAdjacencyRelationship( triangles, surroundingElemEdges );
	}

}

// Edge flipping
void TriangleList::flipEdgeOfTwoTriangles( const std::vector<int>& triangles, const int nodeIDStart , const int nodeIDEnd ){

		if( static_cast<int>( triangles.size() ) != 2 ){
			OutputFiles::m_logFile << " Error : Number of element " << static_cast<int>( triangles.size() ) << " is not equal to 2 in " << __FUNCTION__ << " ." << std::endl;
			exit(1);
		}

		// It is confirmed that the first triangle includes the first node 
		int triangleID = -1;
		int triangleIDAdj = -1;
		for( int iNode = 0; iNode < 3; ++iNode ){
			if( m_triangles[ triangles[0] ].getNodeID( iNode ) == nodeIDStart ){
				triangleID = triangles[0];
				triangleIDAdj = triangles[1];
				break;
			}
			else if( m_triangles[ triangles[1] ].getNodeID( iNode ) == nodeIDStart ){
				triangleID = triangles[1];
				triangleIDAdj = triangles[0];
				break;
			}
		}
		if( triangleID == -1 || triangleIDAdj == -1 ){
			OutputFiles::m_logFile << " Error : Neither triangle " << triangles[0] << " nor " << triangles[1] << " include the start node " << nodeIDStart << std::endl;
			exit(1);
		}

		int edgeID = 0;
		bool found(false);
		for( ; edgeID < 3; ++edgeID ){
			if( m_triangles[triangleID].getAdjacentTriangleID(edgeID) == triangleIDAdj ){
				found = true;
				break;
			}
		}
		if( !found ){
			OutputFiles::m_logFile << " Error : Triangle " <<  triangleID << " is NOT adjacent to " << triangleIDAdj << " ." << std::endl;
			exit(1);
		}

		const int nodeIDsOrg[3] = {
			m_triangles[triangleID].getNodeID( 0 ),
			m_triangles[triangleID].getNodeID( 1 ),
			m_triangles[triangleID].getNodeID( 2 )
		};

		const int adjTriangleIDsOrg[3] = {
			m_triangles[triangleID].getAdjacentTriangleID( 0 ),
			m_triangles[triangleID].getAdjacentTriangleID( 1 ),
			m_triangles[triangleID].getAdjacentTriangleID( 2 )
		};

		const int edgeIDsOfAdjTriangleOrg[3] = {
			m_triangles[triangleID].getEdgeIDOfAdjacentTriangle( 0 ),
			m_triangles[triangleID].getEdgeIDOfAdjacentTriangle( 1 ),
			m_triangles[triangleID].getEdgeIDOfAdjacentTriangle( 2 )
		};

		const int nodeIDsAdjOrg[3] = {
			m_triangles[triangleIDAdj].getNodeID( 0 ),
			m_triangles[triangleIDAdj].getNodeID( 1 ),
			m_triangles[triangleIDAdj].getNodeID( 2 )
		};

		const int adjTriangleIDsAdjOrg[3] = {
			m_triangles[triangleIDAdj].getAdjacentTriangleID( 0 ),
			m_triangles[triangleIDAdj].getAdjacentTriangleID( 1 ),
			m_triangles[triangleIDAdj].getAdjacentTriangleID( 2 )
		};

		const int edgeIDsOfAdjTriangleAdjOrg[3] = {
			m_triangles[triangleIDAdj].getEdgeIDOfAdjacentTriangle( 0 ),
			m_triangles[triangleIDAdj].getEdgeIDOfAdjacentTriangle( 1 ),
			m_triangles[triangleIDAdj].getEdgeIDOfAdjacentTriangle( 2 )
		};

#ifdef _DEBUG_WRITE
		std::cout << "Original triangle first" << std::endl;
		std::cout << "triangleID : " << triangleID << std::endl;
		std::cout << "NodeID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getNodeID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "AdjacentTriangleID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getAdjacentTriangleID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "EdgeIDOfAdjacentTriangle : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getEdgeIDOfAdjacentTriangle(iNode) << " ";
		}
		std::cout << std::endl;

		std::cout << "Original triangle second" << std::endl;
		std::cout << "triangleID : " << triangleIDAdj << std::endl;
		std::cout << "NodeID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleIDAdj].getNodeID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "AdjacentTriangleID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleIDAdj].getAdjacentTriangleID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "EdgeIDOfAdjacentTriangle : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleIDAdj].getEdgeIDOfAdjacentTriangle(iNode) << " ";
		}
		std::cout << std::endl;
#endif

		const int edgeIDAdj = edgeIDsOfAdjTriangleOrg[edgeID];

		//---------------------------------------------------------------------
		m_triangles[triangleID].setNodeID( 0, nodeIDStart );
		m_triangles[triangleID].setNodeID( 1, nodeIDsOrg[edgeID % 3]);
		m_triangles[triangleID].setNodeID( 2, nodeIDEnd );

		const int outerAdjTriangleID1 = adjTriangleIDsOrg[ ( edgeID + 2 ) % 3 ];
		const int outerAdjTriangleID2 = adjTriangleIDsAdjOrg[ ( edgeIDAdj + 1 ) % 3 ];
		m_triangles[triangleID].setAdjacentTriangleID( 0, outerAdjTriangleID1 );
		m_triangles[triangleID].setAdjacentTriangleID( 1, outerAdjTriangleID2 );
		m_triangles[triangleID].setAdjacentTriangleID( 2, triangleIDAdj );

		const int edgeIDOfOuterAdjacentTriangle1 = edgeIDsOfAdjTriangleOrg[ ( edgeID + 2 ) % 3 ];
		const int edgeIDOfOuterAdjacentTriangle2 = edgeIDsOfAdjTriangleAdjOrg[ ( edgeIDAdj + 1 ) % 3 ];
		m_triangles[triangleID].setEdgeIDOfAdjacentTriangle( 0, edgeIDOfOuterAdjacentTriangle1 );
		m_triangles[triangleID].setEdgeIDOfAdjacentTriangle( 1, edgeIDOfOuterAdjacentTriangle2 );
		m_triangles[triangleID].setEdgeIDOfAdjacentTriangle( 2, 0 );

		if( outerAdjTriangleID1 >= 0 ){
			m_triangles[outerAdjTriangleID1].setAdjacentTriangleID( edgeIDOfOuterAdjacentTriangle1, triangleID );
			m_triangles[outerAdjTriangleID1].setEdgeIDOfAdjacentTriangle( edgeIDOfOuterAdjacentTriangle1, 0 );
		}

		if( outerAdjTriangleID2 >= 0 ){
			m_triangles[outerAdjTriangleID2].setAdjacentTriangleID( edgeIDOfOuterAdjacentTriangle2, triangleID );
			m_triangles[outerAdjTriangleID2].setEdgeIDOfAdjacentTriangle( edgeIDOfOuterAdjacentTriangle2, 1 );
		}
		//---------------------------------------------------------------------

		//---------------------------------------------------------------------
		m_triangles[triangleIDAdj].setNodeID( 0, nodeIDStart );
		m_triangles[triangleIDAdj].setNodeID( 1, nodeIDEnd );
		m_triangles[triangleIDAdj].setNodeID( 2, nodeIDsAdjOrg[edgeIDAdj % 3] );

		const int outerAdjTriangleIDAdj1 = adjTriangleIDsAdjOrg[ ( edgeIDAdj + 2 ) % 3 ];
		const int outerAdjTriangleIDAdj2 = adjTriangleIDsOrg[ ( edgeID + 1 ) % 3 ];
		m_triangles[triangleIDAdj].setAdjacentTriangleID( 0, triangleID );
		m_triangles[triangleIDAdj].setAdjacentTriangleID( 1, outerAdjTriangleIDAdj1 );
		m_triangles[triangleIDAdj].setAdjacentTriangleID( 2, outerAdjTriangleIDAdj2 );

		const int edgeIDOfOuterAdjacentTriangleAdj1 = edgeIDsOfAdjTriangleAdjOrg[ ( edgeIDAdj + 2 ) % 3 ];
		const int edgeIDOfOuterAdjacentTriangleAdj2 = edgeIDsOfAdjTriangleOrg[ ( edgeID + 1 ) % 3 ];
		m_triangles[triangleIDAdj].setEdgeIDOfAdjacentTriangle( 0, 2 );
		m_triangles[triangleIDAdj].setEdgeIDOfAdjacentTriangle( 1, edgeIDOfOuterAdjacentTriangleAdj1 );
		m_triangles[triangleIDAdj].setEdgeIDOfAdjacentTriangle( 2, edgeIDOfOuterAdjacentTriangleAdj2 );

		if( outerAdjTriangleIDAdj1 >= 0 ){
			m_triangles[outerAdjTriangleIDAdj1].setAdjacentTriangleID( edgeIDOfOuterAdjacentTriangleAdj1, triangleIDAdj );
			m_triangles[outerAdjTriangleIDAdj1].setEdgeIDOfAdjacentTriangle( edgeIDOfOuterAdjacentTriangleAdj1, 1 );
		}

		if( outerAdjTriangleIDAdj2 >= 0 ){
			m_triangles[outerAdjTriangleIDAdj2].setAdjacentTriangleID( edgeIDOfOuterAdjacentTriangleAdj2, triangleIDAdj );
			m_triangles[outerAdjTriangleIDAdj2].setEdgeIDOfAdjacentTriangle( edgeIDOfOuterAdjacentTriangleAdj2, 2 );
		}
		//---------------------------------------------------------------------

#ifdef _DEBUG_WRITE
		std::cout << "Modfied triangle first" << std::endl;
		std::cout << "triangleID : " << triangleID << std::endl;
		std::cout << "NodeID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getNodeID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "AdjacentTriangleID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getAdjacentTriangleID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "EdgeIDOfAdjacentTriangle : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getEdgeIDOfAdjacentTriangle(iNode) << " ";
		}
		std::cout << std::endl;

		std::cout << "Modfied triangle second" << std::endl;
		std::cout << "triangleID : " << triangleIDAdj << std::endl;
		std::cout << "NodeID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleIDAdj].getNodeID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "AdjacentTriangleID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleIDAdj].getAdjacentTriangleID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "EdgeIDOfAdjacentTriangle : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleIDAdj].getEdgeIDOfAdjacentTriangle(iNode) << " ";
		}
		std::cout << std::endl;
#endif

		//---------------------------------------------------------------------
		m_nodeList.getPointerToNode( nodeIDStart )->addTriangleIDs( triangleIDAdj );
		m_nodeList.getPointerToNode( nodeIDEnd )->addTriangleIDs( triangleID );
		m_nodeList.getPointerToNode( nodeIDsOrg[edgeID % 3] )->eraseTriangleIDs( triangleIDAdj );
		m_nodeList.getPointerToNode( nodeIDsAdjOrg[edgeIDAdj % 3] )->eraseTriangleIDs( triangleID );
		//---------------------------------------------------------------------
			
#ifdef _DEBUG_WRITE
		for( int i = 0; i < 3; ++i ){
			std::cout << "nodeIDsOrg[" << i << "] : " << nodeIDsOrg[i] << " ";
			std::vector<int> tempVec = m_nodeList.getPointerToNode( nodeIDsOrg[i] )->getTriangleIDs();
			for( std::vector<int>::iterator itr = tempVec.begin(); itr != tempVec.end(); ++itr ){
				std::cout << *itr << " ";
			}
			std::cout << std::endl;
		}
		std::cout << "nodeIDsAdjOrg[( edgeIDAdj + 2 ) % 3] : " << nodeIDsAdjOrg[( edgeIDAdj + 2 ) % 3] << " ";
		std::vector<int> tempVec = m_nodeList.getPointerToNode( nodeIDsAdjOrg[( edgeIDAdj + 2 ) % 3] )->getTriangleIDs();
		for( std::vector<int>::iterator itr = tempVec.begin(); itr != tempVec.end(); ++itr ){
			std::cout << *itr << " ";
		}
		std::cout << std::endl;
#endif

}

// Delete the data relating to the polygon
void TriangleList::deleteDataRelatingToPolygon( const std::vector<int>& triangles ){

	for( std::vector<int>::const_iterator itr = triangles.begin(); itr != triangles.end(); ++itr ){
		const int nodeIDsOrg[3] = {
			m_triangles[*itr].getNodeID( 0 ),
			m_triangles[*itr].getNodeID( 1 ),
			m_triangles[*itr].getNodeID( 2 )
		};

		const int adjTriangleIDsOrg[3] = {
			m_triangles[*itr].getAdjacentTriangleID( 0 ),
			m_triangles[*itr].getAdjacentTriangleID( 1 ),
			m_triangles[*itr].getAdjacentTriangleID( 2 )
		};

		const int edgeIDsOfAdjTriangleOrg[3] = {
			m_triangles[*itr].getEdgeIDOfAdjacentTriangle( 0 ),
			m_triangles[*itr].getEdgeIDOfAdjacentTriangle( 1 ),
			m_triangles[*itr].getEdgeIDOfAdjacentTriangle( 2 )
		};

		for( int i = 0; i < 3; ++i ){
			m_triangles[*itr].setNodeID( i, -1 );
			m_triangles[*itr].setAdjacentTriangleID( i, -1 );
			m_triangles[*itr].setEdgeIDOfAdjacentTriangle( i, -1 );
			if( adjTriangleIDsOrg[i] >= 0 ){
				m_triangles[ adjTriangleIDsOrg[i] ].setAdjacentTriangleID( edgeIDsOfAdjTriangleOrg[i], -1 );
				m_triangles[ adjTriangleIDsOrg[i] ].setEdgeIDOfAdjacentTriangle( edgeIDsOfAdjTriangleOrg[i], -1 );
			}
			m_nodeList.getPointerToNode( nodeIDsOrg[i] )->eraseTriangleIDs( *itr );
		}
	}
}

// Search nodes of polygon
void TriangleList::searchNodesOfPolygon( const std::vector<int>& triangles, const int nodeIDStart , const int nodeIDEnd, 
	std::vector<int>& nodeGroup1, std::vector<int>& nodeGroup2, std::vector< std::pair<int,int> >& surroundingElemEdges ){


	// For the node group 1
	int nodeIDCur = nodeIDStart;
	int triangleID(-1);
	nodeGroup1.push_back( nodeIDStart );
	for( std::vector<int>::const_iterator itr = triangles.begin(); itr != triangles.end(); ++itr ){
		if( m_triangles[*itr].getLocalNodeID( nodeIDStart ) >= 0 ){
			triangleID = *itr;
			break;
		}
	}

	if( triangleID < 0 ){
		OutputFiles::m_logFile << " Error : Any triangle include the start node " << nodeIDStart << std::endl;
		exit(1);
	}

	const int countMax = static_cast<int>( triangles.size() );

	int icount(0);
	while(true){

		if( nodeIDCur == nodeIDEnd ){
			break;
		}

		if( ++icount > countMax ){
			OutputFiles::m_logFile << " Error : Reach maximum count number ( " << countMax << " ) " << std::endl;
			exit(1);
		}

		int nodeIDLocal = m_triangles[triangleID].getLocalNodeID( nodeIDCur );
		if( nodeIDLocal < 0 || nodeIDLocal >= 3 ){
			OutputFiles::m_logFile << " Error : Local node ID " << nodeIDLocal << " is out of the range." << std::endl;
			exit(1);
		}

		const int triangleIDAdj = m_triangles[triangleID].getAdjacentTriangleID( nodeIDLocal );
		const int nodeIDNext = m_triangles[triangleID].getNodeID( (nodeIDLocal+1) % 3 ); 

		if( triangleIDAdj < -1 ){
			// Domain boundary
			nodeGroup1.push_back( nodeIDNext );
			triangleID = m_triangles[triangleID].getAdjacentTriangleID( (nodeIDLocal+1)%3 );
			nodeIDCur = nodeIDNext;
		}else if( std::find( triangles.begin(), triangles.end(), triangleIDAdj ) == triangles.end() ){
			// Adjacent element is not included in the polygon
			nodeGroup1.push_back( nodeIDNext );
			surroundingElemEdges.push_back( std::make_pair( triangleIDAdj, m_triangles[triangleID].getEdgeIDOfAdjacentTriangle(nodeIDLocal) ) );
			triangleID = m_triangles[triangleID].getAdjacentTriangleID( (nodeIDLocal+1)%3 );
			nodeIDCur = nodeIDNext;
		}else{
			triangleID = triangleIDAdj;
		}

	}

	// For the node group 2
	nodeIDCur = nodeIDEnd;
	triangleID = -1;
	nodeGroup2.push_back( nodeIDEnd );
	for( std::vector<int>::const_iterator itr = triangles.begin(); itr != triangles.end(); ++itr ){
		if( m_triangles[*itr].getLocalNodeID( nodeIDEnd ) >= 0 ){
			triangleID = *itr;
			break;
		}
	}

	if( triangleID < 0 || nodeIDCur < 0 ){
		OutputFiles::m_logFile << " Error : Any triangle include the end node " << nodeIDEnd << std::endl;
		exit(1);
	}

	icount = 0;
	while(true){

		if( nodeIDCur == nodeIDStart ){
			break;
		}

		if( ++icount> countMax ){
			OutputFiles::m_logFile << " Error : Reach maximum count number ( " << countMax << " ) " << std::endl;
			exit(1);
		}

		int nodeIDLocal = m_triangles[triangleID].getLocalNodeID( nodeIDCur );
		if( nodeIDLocal < 0 || nodeIDLocal >= 3 ){
			OutputFiles::m_logFile << " Error : Local node ID " << nodeIDLocal << " is out of the range." << std::endl;
			exit(1);
		}

		const int triangleIDAdj = m_triangles[triangleID].getAdjacentTriangleID( nodeIDLocal );
		const int nodeIDNext = m_triangles[triangleID].getNodeID( (nodeIDLocal+1) % 3 ); 

		if( triangleIDAdj < -1 ){
			// Domain boundary
			nodeGroup2.push_back( nodeIDNext );
			triangleID = m_triangles[triangleID].getAdjacentTriangleID( (nodeIDLocal+1)%3 );
			nodeIDCur = nodeIDNext;
		}else if( std::find( triangles.begin(), triangles.end(), triangleIDAdj ) == triangles.end() ){
			// Adjacent element is not included in the polygon
			nodeGroup2.push_back( nodeIDNext );
			surroundingElemEdges.push_back( std::make_pair( triangleIDAdj, m_triangles[triangleID].getEdgeIDOfAdjacentTriangle(nodeIDLocal) ) );
			triangleID = m_triangles[triangleID].getAdjacentTriangleID( (nodeIDLocal+1)%3 );
			nodeIDCur = nodeIDNext;
		}else{
			triangleID = triangleIDAdj;
		}

	}

#ifdef _DEBUG_WRITE
	std::cout << "nodeGroup1" << std::endl;
	for( std::vector<int>::iterator itr = nodeGroup1.begin(); itr != nodeGroup1.end(); ++itr ){
		std::cout << *itr << std::endl;
	}

	std::cout << "nodeGroup2" << std::endl;
	for( std::vector<int>::iterator itr = nodeGroup2.begin(); itr != nodeGroup2.end(); ++itr ){
		std::cout << *itr << std::endl;
	}

	std::cout << "surroundingElemEdges" << std::endl;
	for( std::vector< std::pair<int,int> >::iterator itr = surroundingElemEdges.begin(); itr != surroundingElemEdges.end(); ++itr ){
		std::cout << itr->first << " " << itr->second << std::endl;
	}
#endif

}

// Fill triangles in polygon
void TriangleList::fillTrianglesInPolygon( const std::vector<int>& triangles, std::vector<int>& nodeGroup, const int nodeIDStart , const int nodeIDEnd, int& iElem ){

	if( static_cast<int>( nodeGroup.size() ) < 3 ){
		OutputFiles::m_logFile << " Error : Totao number of node ( " << nodeGroup.size()<< " ) of the polygon is less than three." << std::endl;
		exit(1);
	}

	int iNode(0);
	const int countMax = 10000;
	int icount = 0;
	while(static_cast<int>( nodeGroup.size() ) >= 3){

		if( ++icount >= countMax ){
			OutputFiles::m_logFile << " Error : Reach maximum count number ( " << countMax << " ) " << std::endl;
			exit(1);
		}

		const int nodeID0 = nodeGroup[iNode];
		const int nodeID1 = nodeGroup[(iNode+1) % static_cast<int>( nodeGroup.size() )];
		const int nodeID2 = nodeGroup[(iNode+2) % static_cast<int>( nodeGroup.size() )];

		if( Util::locateRightHandSide( m_nodeList.getCoordXYOfPoints(nodeID0),  m_nodeList.getCoordXYOfPoints(nodeID1), m_nodeList.getCoordXYOfPoints(nodeID2) ) ){
			// The third point locate right of the line connecting the first and the second node

			const int triangleIDModified = triangles[iElem];
			m_triangles[triangleIDModified].setNodeID( 0, nodeID0 );
			m_triangles[triangleIDModified].setNodeID( 1, nodeID1 );
			m_triangles[triangleIDModified].setNodeID( 2, nodeID2 );

			m_nodeList.getPointerToNode( nodeID0 )->addTriangleIDs(triangleIDModified);
			m_nodeList.getPointerToNode( nodeID1 )->addTriangleIDs(triangleIDModified);
			m_nodeList.getPointerToNode( nodeID2 )->addTriangleIDs(triangleIDModified);

			nodeGroup.erase( std::find( nodeGroup.begin(), nodeGroup.end(), nodeID1 ) ); 
			++iElem;
		}

		iNode = (iNode + 1) % static_cast<int>( nodeGroup.size() );
	}

}

// Reconstruct adjacency relationship
void TriangleList::reconstructAdjacencyRelationship( const std::vector<int>& triangles, const std::vector< std::pair<int,int> >& surroundingElemEdges ){

	for( std::vector<int>::const_iterator itr0 = triangles.begin(); itr0 != triangles.end(); ++itr0 ){
		for( int iEdge0 = 0; iEdge0 < 3; ++iEdge0 ){
			if( m_triangles[*itr0].getAdjacentTriangleID(iEdge0) >= 0 ){
				continue;// Adjacent element has already been found
			}
			const int node0[2] = { m_triangles[*itr0].getNodeID(iEdge0), m_triangles[*itr0].getNodeID((iEdge0+1)%3) };
			for( std::vector<int>::const_iterator itr1 = triangles.begin(); itr1 != triangles.end(); ++itr1 ){
				if( *itr0 == *itr1 ){
					continue;
				}
				for( int iEdge1 = 0; iEdge1 < 3; ++iEdge1 ){
					const int node1[2] = { m_triangles[*itr1].getNodeID(iEdge1), m_triangles[*itr1].getNodeID((iEdge1+1)%3) };
					if( ( node0[0] == node1[0] && node0[1] == node1[1] ) ||( node0[0] == node1[1] && node0[1] == node1[0] ) ){
						m_triangles[*itr0].setAdjacentTriangleID( iEdge0, *itr1 );
						m_triangles[*itr0].setEdgeIDOfAdjacentTriangle( iEdge0, iEdge1 );
						m_triangles[*itr1].setAdjacentTriangleID( iEdge1, *itr0 );
						m_triangles[*itr1].setEdgeIDOfAdjacentTriangle( iEdge1, iEdge0 );
					}
				}
			}
			for( std::vector< std::pair<int,int> >::const_iterator itrSurr = surroundingElemEdges.begin(); itrSurr != surroundingElemEdges.end(); ++itrSurr ){
				const int elemIDSurr = itrSurr->first;
				const int iEdgeSurr = itrSurr->second;
				const int node1[2] = { m_triangles[elemIDSurr].getNodeID(iEdgeSurr), m_triangles[elemIDSurr].getNodeID((iEdgeSurr+1)%3) };
				if( ( node0[0] == node1[0] && node0[1] == node1[1] ) ||( node0[0] == node1[1] && node0[1] == node1[0] ) ){
					m_triangles[*itr0].setAdjacentTriangleID( iEdge0, elemIDSurr );
					m_triangles[*itr0].setEdgeIDOfAdjacentTriangle( iEdge0, iEdgeSurr );
					m_triangles[elemIDSurr].setAdjacentTriangleID( iEdgeSurr, *itr0 );
					m_triangles[elemIDSurr].setEdgeIDOfAdjacentTriangle( iEdgeSurr, iEdge0 );
				}
			}
		}
	}

#ifdef _DEBUG_WRITE
	std::cout << "Triangles in polygon" << std::endl;
	for( std::vector<int>::const_iterator itr = triangles.begin(); itr != triangles.end(); ++itr ){
		const int triangleID = *itr;
		std::cout << "triangleID : " << triangleID << std::endl;
		std::cout << "NodeID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getNodeID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "AdjacentTriangleID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getAdjacentTriangleID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "EdgeIDOfAdjacentTriangle : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getEdgeIDOfAdjacentTriangle(iNode) << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "Triangles surrounding polygon" << std::endl;
	for( std::vector< std::pair<int,int> >::const_iterator itrSurr = surroundingElemEdges.begin(); itrSurr != surroundingElemEdges.end(); ++itrSurr ){
		const int triangleID = itrSurr->first;
		std::cout << "triangleID : " << triangleID << " iEdge : " << itrSurr->second << std::endl;
		std::cout << "NodeID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getNodeID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "AdjacentTriangleID : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getAdjacentTriangleID(iNode) << " ";
		}
		std::cout << std::endl;
		std::cout << "EdgeIDOfAdjacentTriangle : ";
		for( int iNode = 0; iNode < 3; ++iNode ){
			std::cout << m_triangles[triangleID].getEdgeIDOfAdjacentTriangle(iNode) << " ";
		}
		std::cout << std::endl;
	}
#endif

}

// Get the nodes surroudding the analysis domain for setting the extended region
void TriangleList::getNodesSurroundingAnalysisDomain( const int* additionalCornerNodes, std::vector<CommonParameters::XYZ>& nodeList, std::vector<int>& nodes ) const{

	const BoundaryCurveList* const ptrBoundaryCurveList = &m_boundaryCurveList;

	const AnalysisDomain* const ptrAnalysisDomain = AnalysisDomain::getInstance();

	const double xMin = ptrAnalysisDomain->getMinCoordX();
	const double xMax = ptrAnalysisDomain->getMaxCoordX(); 
	const double yMin = ptrAnalysisDomain->getMinCoordY(); 
	const double yMax = ptrAnalysisDomain->getMaxCoordY(); 

	std::vector< std::pair<double,int> > distanceAndNodeID;

	const double EPS = 1.0e-6;

	const int numOuterBoundary = ptrBoundaryCurveList->getTotalNumberOuterBoundaries();
	for( int iBounOuter = 0; iBounOuter < numOuterBoundary; ++iBounOuter ){

		if( ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getGeologicalType() != CommonParameters::LAND ){
			continue;
		}

		//bool isIntersectPre = false;
		const int numNodes = ptrBoundaryCurveList->getNumNodeOuterBoundary( iBounOuter );
		for( int iNode = 0; iNode < numNodes; ++iNode ){

			const int nodeIDBoundCurveList = ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getNodeID(iNode);	
			const int nodeID = convertNodeIDBoundCurve2Triangle( nodeIDBoundCurveList );
			const CommonParameters::XY coordXY = ptrBoundaryCurveList->getPointCoordOuterBoundary(iBounOuter, iNode);

			if( ptrAnalysisDomain->doesIntersectWithBoundary( coordXY ) ){
				if( ptrAnalysisDomain->doesIntersectWithPlusXEdgeOfAnalysisDomain(coordXY) ){
					const double offset = xMax - xMin;
					distanceAndNodeID.push_back( std::make_pair( offset + coordXY.Y - yMin, nodeID ) );
				}
				else if( ptrAnalysisDomain->doesIntersectWithMinusXEdgeOfAnalysisDomain(coordXY) ){
					const double offset = 2.0 * ( xMax - xMin ) + yMax - yMin;
					distanceAndNodeID.push_back( std::make_pair( offset + yMax - coordXY.Y, nodeID ) );
				}
				else if( ptrAnalysisDomain->doesIntersectWithPlusYEdgeOfAnalysisDomain(coordXY) ){
					const double offset = xMax - xMin + yMax - yMin;
					distanceAndNodeID.push_back( std::make_pair( offset + xMax - coordXY.X, nodeID ) );
				}
				else if( ptrAnalysisDomain->doesIntersectWithMinusYEdgeOfAnalysisDomain(coordXY) ){
					distanceAndNodeID.push_back( std::make_pair( coordXY.X - xMin, nodeID ) );
				}

				//if( ( Control::getInstance() )->isExtendedRegionAdded() && 
				//	( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getGeologicalType() == CommonParameters::LAND ){
				if( ( Control::getInstance() )->isExtendedRegionAdded() ){
					const double zCoord = nodeList[nodeID].Z;
					if( fabs(zCoord) > EPS ){
						OutputFiles::m_logFile << " Warning : The z coordinate ( " << zCoord << " ) of node " << nodeID << " is changed to be 0 since the node belongs to the outer boundary of the analysis domain." << std::endl;
					}
					nodeList[nodeID].Z = 0.0;
				}
			}

		}

	}

	if( additionalCornerNodes[TriangleList::XPLUS_YMINUS_SURFACE] >= 0 ){
		distanceAndNodeID.push_back( std::make_pair( xMax - xMin, additionalCornerNodes[TriangleList::XPLUS_YMINUS_SURFACE] ) );
	}
	if( additionalCornerNodes[TriangleList::XPLUS_YPLUS_SURFACE] >= 0 ){
		distanceAndNodeID.push_back( std::make_pair( xMax - xMin + yMax - yMin, additionalCornerNodes[TriangleList::XPLUS_YPLUS_SURFACE] ) );
	}
	if( additionalCornerNodes[TriangleList::XMINUS_YPLUS_SURFACE] >= 0 ){
		distanceAndNodeID.push_back( std::make_pair( 2.0 * ( xMax - xMin ) + yMax - yMin, additionalCornerNodes[TriangleList::XMINUS_YPLUS_SURFACE] ) );
	}
	if( additionalCornerNodes[TriangleList::XMINUS_YMINUS_SURFACE]>= 0 ){
		distanceAndNodeID.push_back( std::make_pair( 0.0, additionalCornerNodes[TriangleList::XMINUS_YMINUS_SURFACE] ) );
	}

	std::sort(distanceAndNodeID.begin(),distanceAndNodeID.end());

	int iNodeIDPre(-99999);
	for( std::vector< std::pair<double,int> >::iterator itr = distanceAndNodeID.begin(); itr != distanceAndNodeID.end(); ++itr ){
		if( itr->second < 0 || itr->second == iNodeIDPre ){
			continue;
		}
		nodes.push_back(itr->second);
		iNodeIDPre = itr->second;

		//std::cout << itr->first << " " << itr->second << std::endl;
	}

}

// Edit height of nodes deviating from the ones of neighbors
void TriangleList::editDeviatingHeight(){

	const double threthould = Control::getInstance()->getThresholdDeviationOfHegith();
	if( Control::getInstance()->getThresholdDeviationOfHegith() < 0.0 ){
		return;
	}

	OutputFiles::m_logFile << "# Edit height of nodes deviating from the ones of neighbors" << std::endl;

	const int numNodes = m_nodeList.getTotalNumberOfNode();
	double* deviation = new double[numNodes];

	calcDeviationOfHeight(deviation);

	for( int iNode = 0; iNode < numNodes; ++iNode ){
		if( deviation[iNode] > threthould ){
			switch( m_nodeList.getPointerToNode(iNode)->getLocation() ){
				case Node::SEA:// Go through
				case Node::LAND:
					OutputFiles::m_logFile << " Warning : Deviation of node ( " << iNode << " ) is " << deviation[iNode] << " and exceeds the threthold value." << std::endl;
					OutputFiles::m_logFile << "           Its height is calculated again from the ones of neighbor nodes." << std::endl;
					m_nodeList.getPointerToNode(iNode)->setCoordZ( calcHeightFromTheOnesOfNeighbors(iNode) );
					break;
				case Node::COAST_LINE:// Go through
				case Node::LAKE_LINE:// Go through
				case Node::LAKE:
					break;
				case Node::UNKNOWN:
					OutputFiles::m_logFile << " Error : Location of the node " << iNode << " is unknowon." << std::endl;
					exit(1);
					break;
				default:
					OutputFiles::m_logFile << " Error : Unknown location." << std::endl;
					exit(1);
					break;
			}
		}
	}

	delete [] deviation;

}

// Calculate height from the ones of neighbors
double TriangleList::calcHeightFromTheOnesOfNeighbors( const int nodeID ) const{

	const std::vector<int> triangleIDs = m_nodeList.getTriangleIDs(nodeID);
	 
	std::vector<int> neigoborNodes;	 
	for( std::vector<int>::const_iterator itr = triangleIDs.begin(); itr != triangleIDs.end(); ++itr ){
		for( int iNode = 0; iNode < 3; ++iNode ){
			const int nodeIDAdj = m_triangles[*itr].getNodeID(iNode);
			if( nodeID != nodeIDAdj ){
				neigoborNodes.push_back(nodeIDAdj);
			}
		}
	}

	std::sort(neigoborNodes.begin(), neigoborNodes.end());
	neigoborNodes.erase( std::unique( neigoborNodes.begin(), neigoborNodes.end() ), neigoborNodes.end() );

	const CommonParameters::XYZ coord0 = m_nodeList.getCoordXYZOfPoints(nodeID);
	const double distanceUsedToAvoidTooSmallDenominator = ( Control::getInstance() )->getDistanceUsedToAvoidTooSmallDenominator();

	double zCoord(0.0);
	double weightSum(0.0);
	for( std::vector<int>::const_iterator itr = neigoborNodes.begin(); itr != neigoborNodes.end(); ++itr ){
		const CommonParameters::XYZ coord = m_nodeList.getCoordXYZOfPoints(*itr);
		const double weight = 1.0 / ( distanceUsedToAvoidTooSmallDenominator + hypot( coord.X - coord0.X, coord.Y - coord0.Y ) );
		zCoord += weight * coord.Z;
		weightSum += weight;
	}	

	return zCoord / weightSum;
}

// Calculate deviation of height
void TriangleList::calcDeviationOfHeight( double* deviation ) const{

	const int numNodes = m_nodeList.getTotalNumberOfNode();

	int* numNeighNodes = new int[numNodes];
	for( int iNode = 0; iNode < numNodes; ++iNode ){
		numNeighNodes[iNode] = 0;
		deviation[iNode] = 0.0;
	}

	for( std::vector< Triangle >::const_iterator itr = m_triangles.begin(); itr != m_triangles.end(); ++itr ){
		int nodeID[3];
		CommonParameters::XYZ coord[3];
		for( int iNode = 0; iNode < 3; ++iNode ){
			nodeID[iNode] = itr->getNodeID(iNode);
			coord[iNode] = m_nodeList.getCoordXYZOfPoints(nodeID[iNode]);			
		}

		deviation[nodeID[0]] += (coord[0].Z - coord[1].Z) / hypot( coord[0].X - coord[1].X, coord[0].Y - coord[1].Y );
		deviation[nodeID[0]] += (coord[0].Z - coord[2].Z) / hypot( coord[0].X - coord[2].X, coord[0].Y - coord[2].Y );
		numNeighNodes[nodeID[0]] += 2;

		deviation[nodeID[1]] += (coord[1].Z - coord[2].Z) / hypot( coord[1].X - coord[2].X, coord[1].Y - coord[2].Y );
		deviation[nodeID[1]] += (coord[1].Z - coord[0].Z) / hypot( coord[1].X - coord[0].X, coord[1].Y - coord[0].Y );
		numNeighNodes[nodeID[1]] += 2;

		deviation[nodeID[2]] += (coord[2].Z - coord[0].Z) / hypot( coord[2].X - coord[0].X, coord[2].Y - coord[0].Y );
		deviation[nodeID[2]] += (coord[2].Z - coord[1].Z) / hypot( coord[2].X - coord[1].X, coord[2].Y - coord[1].Y );
		numNeighNodes[nodeID[2]] += 2;
	}

	for( int iNode = 0; iNode < numNodes; ++iNode ){
		deviation[iNode] = fabs( deviation[iNode] / static_cast<double>(numNeighNodes[iNode] ) );
	}

	delete [] numNeighNodes;

}

// Make boundary curve for the surface meshes of the sea and lakes
void TriangleList::makeBoundaryCurveForSeaAndLakeMeshes( std::vector<CommonParameters::XYZ>& nodeList, std::map<int,int>*& outerBounNode2MeshNode, std::map<int,int>*& innerBounNode2MeshNode, std::map<int,int>*& subInnerBounNode2MeshNode ){

	const BoundaryCurveList* const ptrBoundaryCurveList = &m_boundaryCurveList;

	std::map<int,int> alreadyInserted;

	// Outer boundary
	const int numOuterBoundary = ptrBoundaryCurveList->getTotalNumberOuterBoundaries();
	if( outerBounNode2MeshNode != NULL ){
		delete [] outerBounNode2MeshNode;
		outerBounNode2MeshNode = NULL;
	}
	if( numOuterBoundary > 0 ){
		outerBounNode2MeshNode = new std::map<int,int>[numOuterBoundary];
		for( int iBoun = 0; iBoun < numOuterBoundary; ++iBoun ){
			const BoundaryCurveOuter* const ptrBound = ptrBoundaryCurveList->getPointerToOuterBoundary(iBoun);
			if( ptrBound->getGeologicalType() != CommonParameters::SEA && ptrBound->getGeologicalType() != CommonParameters::LAKE ){
				continue;// Treat only the sea or lakes
			}
			const int numNodesBoun = ptrBound->getNumOfPoints();
			for( int iNodeBoun = 0; iNodeBoun < numNodesBoun; ++iNodeBoun ){
				const int nodeIDBoun = ptrBound->getNodeID(iNodeBoun);
				const int nodeIDTriangle = convertNodeIDBoundCurve2Triangle(nodeIDBoun);

				if( m_nodeList.getPointerToNode(nodeIDTriangle)->getLocation() == Node::COAST_LINE || m_nodeList.getPointerToNode(nodeIDTriangle)->getLocation() == Node::LAKE_LINE ){
					outerBounNode2MeshNode[iBoun].insert(std::make_pair(nodeIDBoun, nodeIDTriangle));
				}
				else{
					const CommonParameters::XY coordXY = ptrBoundaryCurveList->getPointCoordOuterBoundary(iBoun, iNodeBoun);
					const CommonParameters::XYZ coordXYZ = { coordXY.X, coordXY.Y, 0.0 };
					nodeList.push_back( coordXYZ );
					std::map<int,int>::const_iterator itrAlreadyInserted = alreadyInserted.find(nodeIDBoun);
					if( itrAlreadyInserted != alreadyInserted.end() ){
						outerBounNode2MeshNode[iBoun].insert(std::make_pair(nodeIDBoun, itrAlreadyInserted->second));
					}
					else{
						std::pair<int,int> inserted(nodeIDBoun, static_cast<int>(nodeList.size())-1);
						outerBounNode2MeshNode[iBoun].insert(inserted);
						alreadyInserted.insert(inserted);
					}
				}
			}
		}
	}

//#ifdef _DEBUG_WRITE
	//for( int iBoun = 0; iBoun < numOuterBoundary; ++iBoun ){
	//	for( std::map<int,int>::const_iterator itr = outerBounNode2MeshNode[iBoun].begin(); itr != outerBounNode2MeshNode[iBoun].end(); ++itr ){
	//		std::cout << "iBoun first second : " << iBoun << " " << itr->first << " " << itr->second << std::endl;
	//	}
	//}
//#endif

	// Inner boundary
	const int numInnerBoundary = ptrBoundaryCurveList->getTotalNumberInnerBoundaries();
	if( innerBounNode2MeshNode != NULL ){
		delete [] innerBounNode2MeshNode;
		innerBounNode2MeshNode = NULL;
	}
	if( numInnerBoundary > 0 ){
		innerBounNode2MeshNode = new std::map<int,int>[numInnerBoundary];
		for( int iBoun = 0; iBoun < numInnerBoundary; ++iBoun ){
			const BoundaryCurveInner* const ptrBound = ptrBoundaryCurveList->getPointerToInnerBoundary(iBoun);
			//if( ptrBound->getGeologicalType() != CommonParameters::SEA && ptrBound->getGeologicalType() != CommonParameters::LAKE ){
			//	continue;// Treat only the sea or lakes
			//}
			const int numNodesBoun = ptrBound->getNumOfPoints();
			for( int iNodeBoun = 0; iNodeBoun < numNodesBoun; ++iNodeBoun ){
				const int nodeIDBoun = ptrBound->getNodeID(iNodeBoun);
				const int nodeIDTriangle = convertNodeIDBoundCurve2Triangle(nodeIDBoun);

				if( m_nodeList.getPointerToNode(nodeIDTriangle)->getLocation() == Node::COAST_LINE || m_nodeList.getPointerToNode(nodeIDTriangle)->getLocation() == Node::LAKE_LINE ){
					innerBounNode2MeshNode[iBoun].insert(std::make_pair(nodeIDBoun, nodeIDTriangle));
				}
				else{
					std::map<int,int>::const_iterator itrAlreadyInserted = alreadyInserted.find(nodeIDBoun);
					if( itrAlreadyInserted != alreadyInserted.end() ){
						innerBounNode2MeshNode[iBoun].insert(std::make_pair(nodeIDBoun, itrAlreadyInserted->second));
					}
					else{
						const CommonParameters::XY coordXY = ptrBoundaryCurveList->getPointCoordInnerBoundary(iBoun, iNodeBoun);
						const CommonParameters::XYZ coordXYZ = { coordXY.X, coordXY.Y, 0.0 };
						nodeList.push_back( coordXYZ );
						std::pair<int,int> inserted(nodeIDBoun, static_cast<int>(nodeList.size())-1);
						innerBounNode2MeshNode[iBoun].insert(inserted);
						alreadyInserted.insert(inserted);
					}
				}
			}
		}
	}

	// Sub-inner boundary
	const int numSubInnerBoundary = ptrBoundaryCurveList->getTotalNumberSubInnerBoundaries();
	if( subInnerBounNode2MeshNode != NULL ){
		delete [] subInnerBounNode2MeshNode;
		subInnerBounNode2MeshNode = NULL;
	}
	if( numSubInnerBoundary > 0 ){
		subInnerBounNode2MeshNode = new std::map<int,int>[numSubInnerBoundary];
		for( int iBoun = 0; iBoun < numSubInnerBoundary; ++iBoun ){
			const BoundaryCurveSubInner* const ptrBound = ptrBoundaryCurveList->getPointerToSubInnerBoundary(iBoun);
			//if( ptrBound->getGeologicalType() != CommonParameters::SEA && ptrBound->getGeologicalType() != CommonParameters::LAKE ){
			//	continue;// Treat only the sea or lakes
			//}
			const int numNodesBoun = ptrBound->getNumOfPoints();
			for( int iNodeBoun = 0; iNodeBoun < numNodesBoun; ++iNodeBoun ){
				const int nodeIDBoun = ptrBound->getNodeID(iNodeBoun);
				const int nodeIDTriangle = convertNodeIDBoundCurve2Triangle(nodeIDBoun);

				if( m_nodeList.getPointerToNode(nodeIDTriangle)->getLocation() == Node::COAST_LINE || m_nodeList.getPointerToNode(nodeIDTriangle)->getLocation() == Node::LAKE_LINE ){
					subInnerBounNode2MeshNode[iBoun].insert(std::make_pair(nodeIDBoun, nodeIDTriangle));
				}
				else{
					std::map<int,int>::const_iterator itrAlreadyInserted = alreadyInserted.find(nodeIDBoun);
					if( itrAlreadyInserted != alreadyInserted.end() ){
						subInnerBounNode2MeshNode[iBoun].insert(std::make_pair(nodeIDBoun, itrAlreadyInserted->second));
					}
					else{
						const CommonParameters::XY coordXY = ptrBoundaryCurveList->getPointCoordSubInnerBoundary(iBoun, iNodeBoun);
						const CommonParameters::XYZ coordXYZ = { coordXY.X, coordXY.Y, 0.0 };
						nodeList.push_back( coordXYZ );
						std::pair<int,int> inserted(nodeIDBoun, static_cast<int>(nodeList.size())-1);
						subInnerBounNode2MeshNode[iBoun].insert(inserted);
						alreadyInserted.insert(inserted);
					}
				}
			}
		}
	}
}

// Make facet from the surrounding nodes
void TriangleList::makeFacetFromSurroundingNodes( const std::vector<int>& nodes, const bool reverseOrder, const CommonParameters::Boundary& planeType,
	const CommonParameters::DomainType& domainType, std::vector<CommonParameters::XYZ>& nodeList, std::vector<TriangleList::Facet>& facetList,
	const int iLayer ){

	if( static_cast<int>( nodes.size() ) < 3 ){
		OutputFiles::m_logFile << " Warning : Number of nodes of the faces is less than 3 !!" << std::endl;
		return;
	}

	std::vector<int> nodesWork = nodes;
	if( reverseOrder ){
		std::reverse( nodesWork.begin(), nodesWork.end() );
	}

	const int numNodes = static_cast<int>( nodesWork.size() );
	int* nodeIDsTemp = new int[numNodes];
	CommonParameters::XY* coord2D = new CommonParameters::XY[numNodes];
	std::map<int,int> nodeIDTemp2Org;

	for( int iNode = 0; iNode < numNodes; ++iNode ){
		nodeIDsTemp[iNode] = iNode;
		nodeIDTemp2Org.insert( std::make_pair(iNode, nodesWork[iNode]) );

		const CommonParameters::XYZ coord3D = nodeList[nodesWork[iNode]];
		switch(planeType){
			case CommonParameters::SURFACE:	// Go through
			case CommonParameters::TOP:	// Go through
			case CommonParameters::LAYER:	// Go through
			case CommonParameters::BOT:
				coord2D[iNode].X = coord3D.X;
				coord2D[iNode].Y = coord3D.Y;
				break;
			case CommonParameters::YZ_MINUS:// Go through
			case CommonParameters::YZ_PLUS:
				coord2D[iNode].X = coord3D.Y;
				coord2D[iNode].Y = coord3D.Z;
				break;
			case CommonParameters::ZX_MINUS:// Go through
			case CommonParameters::ZX_PLUS:
				coord2D[iNode].X = coord3D.Z;
				coord2D[iNode].Y = coord3D.X;
				break;
			default:
				OutputFiles::m_logFile << " Error : Unknown plane type : " << m_locationOfPlane << std::endl;
				exit(1);
				break;
		}
	}

	TriangleList triangleListForWork;
	if( planeType == CommonParameters::LAYER ){
		triangleListForWork.setPlaneLocation( planeType, iLayer );
	}else{
		triangleListForWork.setPlaneLocation( planeType );
	}
	triangleListForWork.setDataOfOneOuterBoundaryCurve( numNodes, nodeIDsTemp, coord2D, domainType );
	triangleListForWork.createSurfaceTrianglesOnAdditionalPlane();

	const double coordZ = planeType == CommonParameters::LAYER ? (Control::getInstance())->getDepthLayerInterfaces(iLayer) : -1.0;

	int nodeIDCur = static_cast<int>( nodeList.size() );
	const int numNodesTemp = triangleListForWork.m_nodeList.getTotalNumberOfNode();
	for( int iNode = 0; iNode < numNodesTemp; ++iNode ){
		if( iNode < 3 ){// Exclude super triangles
			continue;
		}
		const int nodeID = iNode - 3;
		if( nodeIDTemp2Org.find( nodeID ) != nodeIDTemp2Org.end() ){// Already inserted
			continue;
		}
		if( planeType == CommonParameters::LAYER ){
			nodeList.push_back( Util::calcCoordOf3DModel( triangleListForWork.m_nodeList.getCoordXYOfPoints( iNode ), planeType, coordZ ) );
		}
		else{
			nodeList.push_back( Util::calcCoordOf3DModel( triangleListForWork.m_nodeList.getCoordXYOfPoints( iNode ), planeType ) );
		}
		nodeIDTemp2Org.insert( std::make_pair(nodeID, nodeIDCur++) );
	}

	// Surface of the earth
	for( std::vector< Triangle >::const_iterator itr = triangleListForWork.m_triangles.begin(); itr != triangleListForWork.m_triangles.end(); ++itr ){

		if( itr->getDomainType() == CommonParameters::OUTSIDE_OF_DOMAIN || itr->getDomainType() == CommonParameters::UNKNOWN ){// Exclude super triangles
			continue;
		}

		TriangleList::PolygonNodes nodes;
		for( int iNode = 0; iNode < 3; ++iNode ){
			const int nodeIDTemp = itr->getNodeID(iNode) - 3;
			std::map<int,int>::const_iterator itr = nodeIDTemp2Org.find(nodeIDTemp);
			if( itr == nodeIDTemp2Org.end() ){// Not found
				OutputFiles::m_logFile << " Error : Node ID " << nodeIDTemp << " cannot be found !! " << std::endl;				
				exit(1);
			}
			nodes.push_back( itr->second );
		}
		TriangleList::Facet facetTmp;
		facetTmp.polygons.push_back( nodes );
		facetList.push_back( facetTmp );

	}

	delete [] nodeIDsTemp;
	delete [] coord2D;

}

// Search the type of the side boundary where specified nodes locate on
CommonParameters::Boundary TriangleList::searchTypeOfSideBoundary( const std::vector<CommonParameters::XYZ>& nodeList, const std::vector<int>& nodes ) const{

	const AnalysisDomain* const ptrAnalysisDomain = AnalysisDomain::getInstance();

	// Check whether all nodes locate on the minus side of the YZ plane
	bool flag(true);
	for( std::vector<int>::const_iterator itr = nodes.begin(); itr != nodes.end(); ++itr ){
		const CommonParameters::XY coord = { nodeList[*itr].X, nodeList[*itr].Y };
		if( !ptrAnalysisDomain->doesIntersectWithMinusXEdgeOfAnalysisDomain(coord) ){
			flag = false;
			break;
		}
	}
	if(flag){
		return CommonParameters::YZ_MINUS;
	}

	// Check whether all nodes locate on the plus side of the YZ plane
	flag = true;
	for( std::vector<int>::const_iterator itr = nodes.begin(); itr != nodes.end(); ++itr ){
		const CommonParameters::XY coord = { nodeList[*itr].X, nodeList[*itr].Y };
		if( !ptrAnalysisDomain->doesIntersectWithPlusXEdgeOfAnalysisDomain(coord) ){
			flag = false;
			break;
		}
	}
	if(flag){
		return CommonParameters::YZ_PLUS;
	}

	// Check whether all nodes locate on the minus side of the ZX plane
	flag = true;
	for( std::vector<int>::const_iterator itr = nodes.begin(); itr != nodes.end(); ++itr ){
		const CommonParameters::XY coord = { nodeList[*itr].X, nodeList[*itr].Y };
		if( !ptrAnalysisDomain->doesIntersectWithMinusYEdgeOfAnalysisDomain(coord) ){
			flag = false;
			break;
		}
	}
	if(flag){
		return CommonParameters::ZX_MINUS;
	}

	// Check whether all nodes locate on the plus side of the ZX plane
	flag = true;
	for( std::vector<int>::const_iterator itr = nodes.begin(); itr != nodes.end(); ++itr ){
		const CommonParameters::XY coord = { nodeList[*itr].X, nodeList[*itr].Y };
		if( !ptrAnalysisDomain->doesIntersectWithPlusYEdgeOfAnalysisDomain(coord) ){
			flag = false;
			break;
		}
	}
	if(flag){
		return CommonParameters::ZX_PLUS;
	}

	OutputFiles::m_logFile << " Error : Any side boundary don't include all of the specified nodes !! " << std::endl;				
	exit(1);

	return CommonParameters::UNDEFINED_BOUNDARY;
}

// Calculate corner coordinates of the land surface
void TriangleList::calcCornerCoordsOfLandSurface( CommonParameters::XYZ coords[4] ){

	const BoundaryCurveList* const ptrBoundaryCurveList = &m_boundaryCurveList;

	const AnalysisDomain* const ptrAnalysisDomain = AnalysisDomain::getInstance();

	const int numOuterBoundary = ptrBoundaryCurveList->getTotalNumberOuterBoundaries();
	for( int iBounOuter = 0; iBounOuter < numOuterBoundary; ++iBounOuter ){

		const int numNodes = ptrBoundaryCurveList->getNumNodeOuterBoundary( iBounOuter );

		for( int iNode = 0; iNode < numNodes; ++iNode ){
			const CommonParameters::XY coordXY = ptrBoundaryCurveList->getPointCoordOuterBoundary(iBounOuter, iNode);
			const bool isIntersect = ptrAnalysisDomain->doesIntersectWithBoundary( coordXY ); 

			if( isIntersect ){
				const int nodeIDBoundCurveList = ( ptrBoundaryCurveList->getPointerToOuterBoundary( iBounOuter ) )->getNodeID(iNode);
				const int nodeID = convertNodeIDBoundCurve2Triangle( nodeIDBoundCurveList );

				if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, AnalysisDomain::XPLUS_YMINUS ) ){
					coords[AnalysisDomain::XPLUS_YMINUS] = m_nodeList.getPointerToNode(nodeID)->getCoordXYZ();
				}
				else if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, AnalysisDomain::XPLUS_YPLUS ) ){
					coords[AnalysisDomain::XPLUS_YPLUS] = m_nodeList.getPointerToNode(nodeID)->getCoordXYZ();
				}
				else if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, AnalysisDomain::XMINUS_YPLUS ) ){
					coords[AnalysisDomain::XMINUS_YPLUS] = m_nodeList.getPointerToNode(nodeID)->getCoordXYZ();
				}
				else if( ptrAnalysisDomain->doesLocateOnGivenCorner( coordXY, AnalysisDomain::XMINUS_YMINUS ) ){
					coords[AnalysisDomain::XMINUS_YMINUS] = m_nodeList.getPointerToNode(nodeID)->getCoordXYZ();
				}
			}

		}
	}

}

//// Add boundary curvers from facet list
//void TriangleList::addBoundaryCurvesFromFacetList( const std::vector<CommonParameters::XYZ>& nodeList, const std::vector<TriangleList::Facet>& facetList ){
//
//	for( std::vector<TriangleList::Facet>::const_iterator itrFacetList = facetList.begin(); itrFacetList != facetList.end(); ++itrFacetList ){
//		BoundaryCurveOuter boundCurve;
//	
//		boundCurve.readBoudaryCurve()
//	}
//
//}
//
//// Calculate coordinate of 3D model from the coordinate of the plane on which 2D mesh is created
//CommonParameters::XYZ TriangleList::calcCoordOf3DModel( const CommonParameters::XY& coord2D ) const{
//
//	CommonParameters::XYZ coord3D = { coord2D.X, coord2D.Y, 0.0 };
//
//	const AnalysisDomain* const ptrAnalysisDomain = AnalysisDomain::getInstance();
//
//	switch(m_locationOfPlane){
//		case CommonParameters::SURFACE:
//			break;
//		case CommonParameters::TOP:
//			coord3D.Z = ptrAnalysisDomain->getMinCoordZ();
//			break;
//		case CommonParameters::BOT:
//			coord3D.Z = ptrAnalysisDomain->getMaxCoordZ();
//			break;
//		case CommonParameters::YZ_MINUS:
//			coord3D.X = ptrAnalysisDomain->getMinCoordX();
//			coord3D.Y = coord2D.X;
//			coord3D.Z = coord2D.Y;
//			break;
//		case CommonParameters::YZ_PLUS:
//			coord3D.X = ptrAnalysisDomain->getMaxCoordX();
//			coord3D.Y = coord2D.X;
//			coord3D.Z = coord2D.Y;
//			break;
//		case CommonParameters::ZX_MINUS:
//			coord3D.Z = coord2D.X;
//			coord3D.X = coord2D.Y;
//			coord3D.Y = ptrAnalysisDomain->getMinCoordY();
//			break;
//		case CommonParameters::ZX_PLUS:
//			coord3D.Z = coord2D.X;
//			coord3D.X = coord2D.Y;
//			coord3D.Y = ptrAnalysisDomain->getMaxCoordY();
//			break;
//		default:
//			OutputFiles::m_logFile << " Error : Unknown plane type : " << m_locationOfPlane << std::endl;
//			exit(1);
//			break;
//	}
//
//	return coord3D;
//
//};
//
//// Coordinate transform to the one on XY plane from the one of the plane on which 2D mesh is created
//CommonParameters::XY TriangleList::coordinateTransform( const CommonParameters::XY& coord2D ) const{
//
//	const CommonParameters::XYZ coord3D = calcCoordOf3DModel(coord2D);
//	const CommonParameters::XY coordTranformed = { coord3D.X, coord3D.Y };
//	return coordTranformed;
//
//}