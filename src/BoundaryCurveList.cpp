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
#include "BoundaryCurveList.h"
#include "CoastLine.h"
#include "CoastLineList.h"
#ifdef _MOD_FOR_NMT
#include "ObservingSiteList.h"
#endif
#include "AnalysisDomain.h"
#include "OutputFiles.h"
#include "Node.h"
#include "Control.h"
#include "LakeList.h"
#include "Util.h"
#include "math.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
const double BoundaryCurveList::m_eps = 1.0e-12;

// Default constructer
BoundaryCurveList::BoundaryCurveList():
	m_outer2InnerBoundaries(NULL),
	m_inner2SubInnerBoundaries(NULL)
{

}

// Destructer
BoundaryCurveList::~BoundaryCurveList(){

	if( m_outer2InnerBoundaries != NULL ){
		delete [] m_outer2InnerBoundaries;
		m_outer2InnerBoundaries = NULL;
	}

	if( m_inner2SubInnerBoundaries != NULL ){
		delete [] m_inner2SubInnerBoundaries;
		m_inner2SubInnerBoundaries = NULL;
	}

}

	//------ For debug >>>>>
#ifdef _DEBUG_WRITE

//void BoundaryCurveList::debugWriteBoundaryCurveListVTK( const std::string& fineName ){
//
//	enum BoundaryType{
//		OUTER_BOUNDARY = 0,
//		INNER_BOUNDARY,
//	};
//
//	// Open output vtk file -----
//	std::ofstream ofsVTK( fineName );
//	if( !ofsVTK ) {
//		std::cerr << "Cannot open file temp.vtk" << std::endl;
//		exit(1);
//	}
//
//	ofsVTK << "# vtk DataFile Version 2.0" << std::endl;
//	ofsVTK << "ForDebug" << std::endl;
//	ofsVTK << "ASCII" << std::endl;
//	ofsVTK << "DATASET UNSTRUCTURED_GRID" << std::endl;
//	
//	ofsVTK.precision(9);
//	ofsVTK << std::fixed;
//
//	int numSegmentsAll = 0;
//
//	// output data to vtk file -----
//	//std::vector<CommonParameters::XY> stack;
//	std::vector<Node> stack;
//
//	for( std::vector<BoundaryCurve>::iterator itr = m_outerBoundaries.begin(); itr != m_outerBoundaries.end(); ++itr ){
//		const int numPointOfSegment = itr->getNumOfPoints();
//		if( numPointOfSegment > 0 ){
//			for( int i = 0; i < numPointOfSegment; ++i ){
//				//stack.push_back( itr->getCoordXYOfPoints( i ) );
//				stack.push_back( itr->getNode2D( i ) );
//			}
//			++numSegmentsAll;
//		}
//	}
//	for( std::vector<BoundaryCurve>::iterator itr = m_innerBoundaries.begin(); itr != m_innerBoundaries.end(); ++itr ){
//		const int numPointOfSegment = itr->getNumOfPoints();
//		if( numPointOfSegment > 0 ){
//			for( int i = 0; i < numPointOfSegment; ++i ){
//				//stack.push_back( itr->getCoordXYOfPoints( i ) );
//				stack.push_back( itr->getNode2D( i ) );
//			}
//			++numSegmentsAll;
//		}
//	}
//
//	const int numPointsAll = static_cast<int>( stack.size() );
//
//	ofsVTK << "POINTS " << numPointsAll << " float" << std::endl;
//	//for( std::vector<CommonParameters::XY>::iterator itr = stack.begin(); itr != stack.end(); ++itr ){	
//	for( std::vector<Node>::iterator itr = stack.begin(); itr != stack.end(); ++itr ){
//		ofsVTK << itr->getCoordX() << " " << itr->getCoordY() << " " << 0.0 << std::endl;
//	}
//
//	ofsVTK << "CELLS " << numSegmentsAll << " " << numSegmentsAll * 2 + numPointsAll << std::endl;
//
//	int icount = 0;
//	for( std::vector<BoundaryCurve>::iterator itr = m_outerBoundaries.begin(); itr != m_outerBoundaries.end(); ++itr ){
//		const int numPointOfSegment = itr->getNumOfPoints();
//		if( numPointOfSegment > 0 ){
//			ofsVTK << numPointOfSegment + 1 << " ";
//			const int icountFirst = icount;
//			for( int i = 0; i < numPointOfSegment; ++i ){
//				ofsVTK << icount++ << " ";
//			}
//			ofsVTK << icountFirst << std::endl;
//		}
//	}
//	for( std::vector<BoundaryCurve>::iterator itr = m_innerBoundaries.begin(); itr != m_innerBoundaries.end(); ++itr ){
//		const int numPointOfSegment = itr->getNumOfPoints();
//		if( numPointOfSegment > 0 ){
//			ofsVTK << numPointOfSegment + 1 << " ";
//			const int icountFirst = icount;
//			for( int i = 0; i < numPointOfSegment; ++i ){
//				ofsVTK << icount++ << " ";
//			}
//			ofsVTK << icountFirst << std::endl;
//		}
//	}
//
//	ofsVTK << "CELL_TYPES " << numSegmentsAll << std::endl;
//	for( int i = 0; i < numSegmentsAll; ++i ){
//		ofsVTK << "4" << std::endl;
//	}
//
//	ofsVTK << "CELL_DATA " << numSegmentsAll << std::endl;
//	ofsVTK << "SCALARS ID int" <<  std::endl;
//	ofsVTK << "LOOKUP_TABLE default" <<  std::endl;
//	icount = 0;
//	for( std::vector<BoundaryCurve>::iterator itr = m_outerBoundaries.begin(); itr != m_outerBoundaries.end(); ++itr ){
//		if( itr->getNumOfPoints() > 0 ){
//			ofsVTK << icount++ << std::endl;
//		}
//	}
//	for( std::vector<BoundaryCurve>::iterator itr = m_innerBoundaries.begin(); itr != m_innerBoundaries.end(); ++itr ){
//		if( itr->getNumOfPoints() > 0 ){
//			ofsVTK << icount++ << std::endl;
//		}
//	}
//
//	ofsVTK << "SCALARS BoundaryType int" <<  std::endl;
//	ofsVTK << "LOOKUP_TABLE default" <<  std::endl;
//	for( std::vector<BoundaryCurve>::iterator itr = m_outerBoundaries.begin(); itr != m_outerBoundaries.end(); ++itr ){
//		if( itr->getNumOfPoints() > 0 ){
//			ofsVTK << OUTER_BOUNDARY << std::endl;
//		}
//	}
//	for( std::vector<BoundaryCurve>::iterator itr = m_innerBoundaries.begin(); itr != m_innerBoundaries.end(); ++itr ){
//		if( itr->getNumOfPoints() > 0 ){
//			ofsVTK << INNER_BOUNDARY << std::endl;
//		}
//	}
//
//	ofsVTK << "SCALARS Type int" <<  std::endl;
//	ofsVTK << "LOOKUP_TABLE default" <<  std::endl;
//	for( std::vector<BoundaryCurve>::iterator itr = m_outerBoundaries.begin(); itr != m_outerBoundaries.end(); ++itr ){
//		if( itr->getNumOfPoints() > 0 ){
//			ofsVTK << itr->getType() << std::endl;
//		}
//	}
//	for( std::vector<BoundaryCurve>::iterator itr = m_innerBoundaries.begin(); itr != m_innerBoundaries.end(); ++itr ){
//		if( itr->getNumOfPoints() > 0 ){
//			ofsVTK << itr->getType() << std::endl;
//		}
//	}
//
//	ofsVTK << "POINT_DATA " << numPointsAll << std::endl;
//	ofsVTK << "SCALARS Fix int" <<  std::endl;
//	ofsVTK << "LOOKUP_TABLE default" <<  std::endl;
//	//for( std::vector<CommonParameters::XY>::iterator itr = stack.begin(); itr != stack.end(); ++itr ){	
//	for( std::vector<Node>::iterator itr = stack.begin(); itr != stack.end(); ++itr ){	
//		ofsVTK << itr->isFixed() << std::endl;
//	}
//
//	ofsVTK.close();
//
//}

#endif
	//------ For debug <<<<<

void BoundaryCurveList::writeBoundaryCurveListVTK( const std::string& fileName ){

	enum BoundaryType{
		OUTER_BOUNDARY = 0,
		INNER_BOUNDARY,
		SUB_INNER_BOUNDARY,
	};

	// Open output vtk file -----
	std::ofstream ofsVTK( fileName.c_str() );
	if( !ofsVTK ) {
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << std::endl;
		exit(1);
	}

	ofsVTK << "# vtk DataFile Version 2.0" << std::endl;
	ofsVTK << "BoundaryList2D" << std::endl;
	ofsVTK << "ASCII" << std::endl;
	ofsVTK << "DATASET UNSTRUCTURED_GRID" << std::endl;
	
	ofsVTK.precision(9);
	ofsVTK << std::fixed;

	// output data to vtk file -----
	std::vector<Node> stack;

	for( std::vector<BoundaryCurveOuter>::iterator itr = m_outerBoundaries.begin(); itr != m_outerBoundaries.end(); ++itr ){
		const int numPointOfSegment = itr->getNumOfPoints();
		if( numPointOfSegment > 0 ){
			for( int i = 0; i < numPointOfSegment; ++i ){
				//stack.push_back( itr->getCoordXYOfPoints( i ) );
				//stack.push_back( *(itr->getPointerToNode(i)) );
				stack.push_back( *m_nodeList.getPointerToNode(itr->getNodeID(i)) );
			}
		}
	}
	for( std::vector<BoundaryCurveInner>::iterator itr = m_innerBoundaries.begin(); itr != m_innerBoundaries.end(); ++itr ){
		const int numPointOfSegment = itr->getNumOfPoints();
		if( numPointOfSegment > 0 ){
			for( int i = 0; i < numPointOfSegment; ++i ){
				//stack.push_back( itr->getCoordXYOfPoints( i ) );
				//stack.push_back( *(itr->getPointerToNode(i)) );
				stack.push_back( *m_nodeList.getPointerToNode(itr->getNodeID(i)) );
			}
		}
	}
	for( std::vector<BoundaryCurveSubInner>::iterator itr = m_subInnerBoundaries.begin(); itr != m_subInnerBoundaries.end(); ++itr ){
		const int numPointOfSegment = itr->getNumOfPoints();
		if( numPointOfSegment > 0 ){
			for( int i = 0; i < numPointOfSegment; ++i ){
				//stack.push_back( itr->getCoordXYOfPoints( i ) );
				//stack.push_back( *(itr->getPointerToNode(i)) );
				stack.push_back( *m_nodeList.getPointerToNode(itr->getNodeID(i)) );
			}
		}
	}

	const int numPointsAll = static_cast<int>( stack.size() );

	ofsVTK << "POINTS " << numPointsAll << " float" << std::endl;
	for( std::vector<Node>::iterator itr = stack.begin(); itr != stack.end(); ++itr ){
		ofsVTK << itr->getCoordX() << " " << itr->getCoordY() << " " << 0.0 << std::endl;
	}

	ofsVTK << "CELLS " << numPointsAll << " " << numPointsAll * 3 << std::endl;

	int icount = 0;
	for( std::vector<BoundaryCurveOuter>::iterator itr = m_outerBoundaries.begin(); itr != m_outerBoundaries.end(); ++itr ){
		const int numPointOfSegment = itr->getNumOfPoints();
		if( numPointOfSegment <= 0 ){
			continue;
		}
		const int IDFirst = icount; 
		for( int i = 0; i < numPointOfSegment - 1; ++i ){
			const int ID = icount++;
			ofsVTK << 2 << " " << ID << " " << ID + 1 << std::endl;
		}
		ofsVTK << 2 << " " << icount++ << " " << IDFirst << std::endl;
	}
	for( std::vector<BoundaryCurveInner>::iterator itr = m_innerBoundaries.begin(); itr != m_innerBoundaries.end(); ++itr ){
		const int numPointOfSegment = itr->getNumOfPoints();
		if( numPointOfSegment <= 0 ){
			continue;
		}
		const int IDFirst = icount; 
		for( int i = 0; i < numPointOfSegment - 1; ++i ){
			const int ID = icount++;
			ofsVTK << 2 << " " << ID << " " << ID + 1 << std::endl;
		}
		ofsVTK << 2 << " " << icount++ << " " << IDFirst << std::endl;
	}
	for( std::vector<BoundaryCurveSubInner>::iterator itr = m_subInnerBoundaries.begin(); itr != m_subInnerBoundaries.end(); ++itr ){
		const int numPointOfSegment = itr->getNumOfPoints();
		if( numPointOfSegment <= 0 ){
			continue;
		}
		const int IDFirst = icount; 
		for( int i = 0; i < numPointOfSegment - 1; ++i ){
			const int ID = icount++;
			ofsVTK << 2 << " " << ID << " " << ID + 1 << std::endl;
		}
		ofsVTK << 2 << " " << icount++ << " " << IDFirst << std::endl;
	}

	ofsVTK << "CELL_TYPES " << numPointsAll << std::endl;
	for( int i = 0; i < numPointsAll; ++i ){
		ofsVTK << "3" << std::endl;
	}

	ofsVTK << "CELL_DATA " << numPointsAll << std::endl;
	ofsVTK << "SCALARS ID int" <<  std::endl;
	ofsVTK << "LOOKUP_TABLE default" <<  std::endl;
	for( int i = 0; i < numPointsAll; ++i ){
		ofsVTK << i << std::endl;
	}

	ofsVTK << "SCALARS BoundaryType int" <<  std::endl;
	ofsVTK << "LOOKUP_TABLE default" <<  std::endl;
	for( std::vector<BoundaryCurveOuter>::iterator itr = m_outerBoundaries.begin(); itr != m_outerBoundaries.end(); ++itr ){
		const int numPointOfSegment = itr->getNumOfPoints();
		if( numPointOfSegment <= 0 ){
			continue;
		}
		for( int i = 0; i < numPointOfSegment; ++i ){
			ofsVTK << OUTER_BOUNDARY << std::endl;
		}
	}
	for( std::vector<BoundaryCurveInner>::iterator itr = m_innerBoundaries.begin(); itr != m_innerBoundaries.end(); ++itr ){
		const int numPointOfSegment = itr->getNumOfPoints();
		if( numPointOfSegment <= 0 ){
			continue;
		}
		for( int i = 0; i < numPointOfSegment; ++i ){
			ofsVTK << INNER_BOUNDARY << std::endl;
		}
	}
	for( std::vector<BoundaryCurveSubInner>::iterator itr = m_subInnerBoundaries.begin(); itr != m_subInnerBoundaries.end(); ++itr ){
		const int numPointOfSegment = itr->getNumOfPoints();
		if( numPointOfSegment <= 0 ){
			continue;
		}
		for( int i = 0; i < numPointOfSegment; ++i ){
			ofsVTK << SUB_INNER_BOUNDARY << std::endl;
		}
	}

	ofsVTK << "SCALARS Type int" <<  std::endl;
	ofsVTK << "LOOKUP_TABLE default" <<  std::endl;
	for( std::vector<BoundaryCurveOuter>::iterator itr = m_outerBoundaries.begin(); itr != m_outerBoundaries.end(); ++itr ){
		const int numPointOfSegment = itr->getNumOfPoints();
		const int itype = itr->getGeologicalType();
		for( int i = 0; i < numPointOfSegment; ++i ){
			ofsVTK << itype << std::endl;
		}
	}
	for( std::vector<BoundaryCurveInner>::iterator itr = m_innerBoundaries.begin(); itr != m_innerBoundaries.end(); ++itr ){
		const int numPointOfSegment = itr->getNumOfPoints();
		const int itype = itr->getGeologicalType();
		for( int i = 0; i < numPointOfSegment; ++i ){
			ofsVTK << itype << std::endl;
		}
	}
	for( std::vector<BoundaryCurveSubInner>::iterator itr = m_subInnerBoundaries.begin(); itr != m_subInnerBoundaries.end(); ++itr ){
		const int numPointOfSegment = itr->getNumOfPoints();
		const int itype = itr->getGeologicalType();
		for( int i = 0; i < numPointOfSegment; ++i ){
			ofsVTK << itype << std::endl;
		}
	}

	ofsVTK << "SCALARS Length float" <<  std::endl;
	ofsVTK << "LOOKUP_TABLE default" <<  std::endl;
	for( std::vector<BoundaryCurveOuter>::iterator itr = m_outerBoundaries.begin(); itr != m_outerBoundaries.end(); ++itr ){
		const int numPointOfSegment = itr->getNumOfPoints();
		for( int i = 0; i < numPointOfSegment - 1; ++i ){
			//const double length = hypot( (itr->getCoordXYOfPoints(i+1)).X - (itr->getCoordXYOfPoints(i)).X , (itr->getCoordXYOfPoints(i+1)).Y - (itr->getCoordXYOfPoints(i)).Y );
			const CommonParameters::XY coord0 = m_nodeList.getCoordXYOfPoints( itr->getNodeID(i) );
			const CommonParameters::XY coord1 = m_nodeList.getCoordXYOfPoints( itr->getNodeID(i+1) );
			const double length = hypot( coord1.X - coord0.X , coord1.Y - coord0.Y );
			ofsVTK << length << std::endl;
		}
		const CommonParameters::XY coord0 = m_nodeList.getCoordXYOfPoints( itr->getNodeID(numPointOfSegment-1) );
		const CommonParameters::XY coord1 = m_nodeList.getCoordXYOfPoints( itr->getNodeID(0) );
		const double length = hypot( coord1.X - coord0.X , coord1.Y - coord0.Y );
		ofsVTK << length << std::endl;
	}
	for( std::vector<BoundaryCurveInner>::iterator itr = m_innerBoundaries.begin(); itr != m_innerBoundaries.end(); ++itr ){
		const int numPointOfSegment = itr->getNumOfPoints();
		for( int i = 0; i < numPointOfSegment - 1; ++i ){
			const CommonParameters::XY coord0 = m_nodeList.getCoordXYOfPoints( itr->getNodeID(i) );
			const CommonParameters::XY coord1 = m_nodeList.getCoordXYOfPoints( itr->getNodeID(i+1) );
			const double length = hypot( coord1.X - coord0.X , coord1.Y - coord0.Y );
			ofsVTK << length << std::endl;
		}
		const CommonParameters::XY coord0 = m_nodeList.getCoordXYOfPoints( itr->getNodeID(numPointOfSegment-1) );
		const CommonParameters::XY coord1 = m_nodeList.getCoordXYOfPoints( itr->getNodeID(0) );
		const double length = hypot( coord1.X - coord0.X , coord1.Y - coord0.Y );
		ofsVTK << length << std::endl;
	}
	for( std::vector<BoundaryCurveSubInner>::iterator itr = m_subInnerBoundaries.begin(); itr != m_subInnerBoundaries.end(); ++itr ){
		const int numPointOfSegment = itr->getNumOfPoints();
		for( int i = 0; i < numPointOfSegment - 1; ++i ){
			const CommonParameters::XY coord0 = m_nodeList.getCoordXYOfPoints( itr->getNodeID(i) );
			const CommonParameters::XY coord1 = m_nodeList.getCoordXYOfPoints( itr->getNodeID(i+1) );
			const double length = hypot( coord1.X - coord0.X , coord1.Y - coord0.Y );
			ofsVTK << length << std::endl;
		}
		const CommonParameters::XY coord0 = m_nodeList.getCoordXYOfPoints( itr->getNodeID(numPointOfSegment-1) );
		const CommonParameters::XY coord1 = m_nodeList.getCoordXYOfPoints( itr->getNodeID(0) );
		const double length = hypot( coord1.X - coord0.X , coord1.Y - coord0.Y );
		ofsVTK << length << std::endl;
	}

	ofsVTK.close();

}

//// Make boudary curves
//void BoundaryCurveList::makeBoundaryCurvesCoastLine(){
//
//	const CoastLineList* const ptrCoastLineList = CoastLineList::getInstance();
//
//	const AnalysisDomain* const ptrAnalysisDomain = AnalysisDomain::getInstance();
//	const double xMin = ptrAnalysisDomain->getMinCoordX(); 
//	const double xMax = ptrAnalysisDomain->getMaxCoordX(); 
//	const double yMin = ptrAnalysisDomain->getMinCoordY(); 
//	const double yMax = ptrAnalysisDomain->getMaxCoordY(); 
//
//	const CommonParameters::XY coordUpperLeft  = { xMax, yMin };
//	const CommonParameters::XY coordUpperRight = { xMax, yMax };
//	const CommonParameters::XY coordLowerRight = { xMin, yMax };
//	const CommonParameters::XY coordLowerLeft  = { xMin, yMin };
//
//	const int numCoastLine = ptrCoastLineList->getNumCoastLines();
//
//	// No coast line
//	if( numCoastLine == 0 ){
//		std::vector<Node> nodes;
//		nodes.push_back( Node( coordUpperLeft, true ) );
//		nodes.push_back( Node( coordUpperRight, true ) );
//		nodes.push_back( Node( coordLowerRight, true ) );
//		nodes.push_back( Node( coordLowerLeft, true ) );
//
//		m_outerBoundaries.push_back( BoundaryCurve( nodes, BoundaryCurve::LAND ) );
//
//		return;
//	}
//
//	//std::vector< std::vector<CommonParameters::XY> > nodeVecArray;
//	std::vector< std::vector<Node> > nodeVecArray;
//
//	//std::vector<CommonParameters::XY> nodeVec;
//	std::vector<Node> nodeVec;
//
//	for( int iCoastLine = 0; iCoastLine < numCoastLine; ++iCoastLine ){
//
//		const CoastLine* const ptrCoastLine = ptrCoastLineList->getPointerToCoastLine(iCoastLine);
//
//		if( !ptrCoastLine->isClosed() ){// Skip if coast line is not closed
//			continue;
//		}
//
//		bool hasIntersects = false;
//		bool locateInsidePreviously = true;
//		int segmentIDFirst = -1;
//
//		const int numPoint = ptrCoastLine->getNumOfPoints();
//		for( int i = 0; i < numPoint; ++i ){
//			CommonParameters::XY coord = ptrCoastLine->getCoordXYOfPoints(i);
//			if( xMin <= coord.X && coord.X <= xMax && yMin <= coord.Y && coord.Y <= yMax ){
//				// This point locate inside of the analysis domain
//
//				if( !locateInsidePreviously ){
//					nodeVec.push_back( Node( ptrAnalysisDomain->getCoordOfIntersectionPoint( ptrCoastLine->getCoordXYOfPoints(i-1), coord ) , true ) );
//					hasIntersects = true;
//				}
//
//				nodeVec.push_back( Node( coord, false ) );
//
//				if( i == 0 ){// First point
//					segmentIDFirst = static_cast<int>( nodeVecArray.size() );
//				}else if( i == numPoint - 1 && hasIntersects ){// Last point
//					if( segmentIDFirst >= 0 ){
//						for( std::vector<Node>::iterator itr = nodeVecArray[segmentIDFirst].begin(); itr != nodeVecArray[segmentIDFirst].end(); ++itr ){
//							nodeVec.push_back( *itr );
//						}
//						nodeVecArray[segmentIDFirst].swap( nodeVec );
//					}else{
//						OutputFiles::m_logFile << "Error : First segment ID of the coast " << iCoastLine << " is negative !!" << std::endl;
//						exit(1);
//					}
//				}
//
//				locateInsidePreviously = true;
//
//			}else{
//				// This point locate out of the analysis domain
//
//				if( locateInsidePreviously ){
//					nodeVec.push_back( Node( ptrAnalysisDomain->getCoordOfIntersectionPoint( ptrCoastLine->getCoordXYOfPoints(i-1), coord ) , true ) );
//					nodeVecArray.push_back( nodeVec );
//					nodeVec.clear();
//					hasIntersects = true;
//				}
//
//				locateInsidePreviously = false;
//			}
//		}
//
//		if( !hasIntersects ){
//#ifdef _DEBUG_WRITE
//			std::cout << "No intersects" << std::endl;
//#endif
//
//			reverse( nodeVec.begin(), nodeVec.end() );
//			m_innerBoundaries.push_back( BoundaryCurve( nodeVec, BoundaryCurve::LAND ) );
//		}
//
//	}
//
//	if( nodeVecArray.empty() ){
//
//		std::vector<Node> nodes;
//		nodes.push_back( Node( coordUpperLeft, true ) );
//		nodes.push_back( Node( coordUpperRight, true ) );
//		nodes.push_back( Node( coordLowerRight, true ) );
//		nodes.push_back( Node( coordLowerLeft, true ) );
//
//		m_outerBoundaries.push_back( BoundaryCurve( nodes, BoundaryCurve::SEA) );
//
//		return;
//
//	}
//
//	typedef std::vector< std::pair< std::pair<double, double>, int > > ArrayType;
//	ArrayType arrayStart;
//	ArrayType arrayEnd;
//
//	int icount(0);
//	for( std::vector< std::vector<Node> >::iterator itr = nodeVecArray.begin(); itr != nodeVecArray.end(); ++itr, ++icount ){
//
//#ifdef _DEBUG_WRITE
//		std::cout << "itr->front() : " << itr->front().getCoordX() << " " << itr->front().getCoordY() << std::endl;
//		std::cout << "itr->back() : " << itr->back().getCoordX() << " " << itr->back().getCoordY() << std::endl;
//#endif
//		const double disStart = ptrAnalysisDomain->calcDistanceOnBoundaryFromUpperLeft( (itr->front()).getCoordXY() );
//		const double disEnd = ptrAnalysisDomain->calcDistanceOnBoundaryFromUpperLeft( (itr->back()).getCoordXY() );
//		arrayStart.push_back( std::make_pair( std::make_pair( disStart , disEnd ) , icount ) );
//		arrayEnd.push_back( std::make_pair( std::make_pair( disEnd , disStart ) , icount ) );
//	}
//
//	std::sort( arrayStart.begin(), arrayStart.end() );// Ascending order
//	std::sort( arrayEnd.begin(), arrayEnd.end() );// // Ascending order
//
//#ifdef _DEBUG_WRITE
//	for( ArrayType::iterator itr = arrayStart.begin(); itr != arrayStart.end(); ++itr ){
//		std::cout << "arrayStart first second : " << (itr->first).first << " " << (itr->first).second << " " << itr->second << std::endl;
//	}
//	for( ArrayType::iterator itr = arrayEnd.begin(); itr != arrayEnd.end(); ++itr ){
//		std::cout << "arrayEnd first second : " << (itr->first).first << " " << (itr->first).second << " " << itr->second << std::endl;
//	}
//#endif
//
//#ifdef _DEBUG_WRITE
//	std::cout << "--- For the land boundary" << std::endl;
//#endif
//
//	//---------------------------------
//	//----- For the land boundary -----
//	//---------------------------------
//	icount = 0;
//	while( !arrayStart.empty() ){
//
//#ifdef _DEBUG_WRITE
//		std::cout << " icount = " << icount << std::endl;
//#endif
//
//		if( icount > 1.0e6 ){
//			OutputFiles::m_logFile << "Error : icount exceed 1.0e6 in the loop !!" << std::endl;
//			exit(1);
//		}
//
//		//std::vector< ArrayType::iterator > stack;
//		std::vector<int> stack;
//		//std::vector<CommonParameters::XY> coords;
//		std::vector<Node> nodes;
//
//		ArrayType::iterator itrFirst = arrayStart.begin();
//		ArrayType::iterator itrLast  = arrayStart.end();
//
//		for( ArrayType::iterator itr = itrFirst; itr != itrLast; ){
//
//			const int vecID = itr->second;
//
//#ifdef _DEBUG_WRITE
//			std::cout << "vecID = " << vecID << std::endl;
//#endif
//
//			nodes.insert( nodes.end(), nodeVecArray[vecID].begin(), nodeVecArray[vecID].end() );
//			stack.push_back( vecID );
//
//#ifdef _DEBUG_WRITE
//			std::cout << "To be Erased " << stack.back() << std::endl;
//#endif
//
//			int edgeIDEnd = ptrAnalysisDomain->getEdgeID( (nodeVecArray[vecID].back()).getCoordXY() );
//			const double disEnd = (itr->first).second;
//
//#ifdef _DEBUG_WRITE
//			std::cout << "edgeIDEnd disEnd : " << edgeIDEnd << " " << disEnd << std::endl;
//#endif
//
//			bool goNext(false);
//
//			for( ArrayType::iterator itrNext = itrFirst; itrNext != itrLast; ++itrNext ){
//
//				const double disStartNext = (itrNext->first).first; 
//
//#ifdef _DEBUG_WRITE
//				std::cout << "disEnd disStartNext : " << disEnd << " " << disStartNext << std::endl;
//#endif
//
//				if( disEnd < disStartNext ){
//
//					if( std::find( stack.begin(), stack.end(), itrNext->second ) != stack.end() ){// This segment has already been found
//#ifdef _DEBUG_WRITE
//						std::cout << "Already found" << std::endl;
//#endif
//						break;
//					}
//
//					const int edgeIDStartNext = ptrAnalysisDomain->getEdgeID( (nodeVecArray[itrNext->second].front()).getCoordXY() );
//
//					while( edgeIDEnd != edgeIDStartNext ){
//						nodes.push_back( Node( ptrAnalysisDomain->getEndCoord( edgeIDEnd ) , true ) ); 
//						edgeIDEnd = ptrAnalysisDomain->getNextEdgeIDCWR( edgeIDEnd );
//#ifdef _DEBUG_WRITE
//						std::cout << "nodes.back() : " << (nodes.back()).getCoordX() << " " << (nodes.back()).getCoordY() << std::endl;
//						std::cout << "edgeIDEnd : " << edgeIDEnd << std::endl;
//#endif
//					}
//
//					goNext = true;
//					itr = itrNext;
//					break;
//				}
//
//			}
//
//			if( goNext ){
//				continue;
//			}
//
//#ifdef _DEBUG_WRITE
//			std::cout << "Last segment" << std::endl;
//#endif
//
//			const int edgeIDStartFirst = ptrAnalysisDomain->getEdgeID( (nodeVecArray[itrFirst->second].front()).getCoordXY() );
//
//#ifdef _DEBUG_WRITE
//			std::cout << "edgeIDStartFirst : " << edgeIDStartFirst << std::endl;
//#endif
//
//			const double disStartFirst = (itrFirst->first).first;
//
//			if( edgeIDStartFirst == edgeIDEnd && disStartFirst < disEnd ){
//				nodes.push_back( Node( ptrAnalysisDomain->getEndCoord( edgeIDEnd ) , true ) ); 
//				edgeIDEnd = ptrAnalysisDomain->getNextEdgeIDCWR( edgeIDEnd );
//			}
//
//			while( edgeIDEnd != edgeIDStartFirst ){
//				nodes.push_back( Node( ptrAnalysisDomain->getEndCoord( edgeIDEnd ) , true ) ); 
//				edgeIDEnd = ptrAnalysisDomain->getNextEdgeIDCWR( edgeIDEnd );
//
//#ifdef _DEBUG_WRITE
//				std::cout << "nodes.back() : " << (nodes.back()).getCoordX() << " " << (nodes.back()).getCoordY() << std::endl;
//				std::cout << "edgeIDEnd : " << edgeIDEnd << std::endl;
//#endif
//
//			}
//			break;			
//
//		}
//
//		m_outerBoundaries.push_back( BoundaryCurve( nodes, BoundaryCurve::LAND ) );
//
//#ifdef _DEBUG_WRITE
//		std::cout << "Remove already found segment from list" << std::endl;
//#endif
//
//		// Remove already found segment from list
//#ifdef _DEBUG_WRITE
//		std::cout << "stack.size() " << stack.size() << std::endl;
//		for( std::vector<int>::iterator itr = stack.begin(); itr != stack.end(); ++itr ){
//			std::cout << *itr << std::endl;
//		}
//#endif
//
//		for( std::vector<int>::iterator itrStack = stack.begin(); itrStack != stack.end(); ++itrStack ){
//
//			bool found(false);
//
//			for( ArrayType::iterator itrArray = arrayStart.begin(); itrArray != arrayStart.end(); ++itrArray ){
//				if( std::find( stack.begin(), stack.end(), itrArray->second ) != stack.end() ){
//#ifdef _DEBUG_WRITE
//					std::cout << "Erase " << itrArray->second << " from arrayStart " << std::endl;
//#endif
//					arrayStart.erase( itrArray );
//					found = true;
//					break;
//				}
//			}
//
//			if( !found ){
//				OutputFiles::m_logFile << "Error : Element whose second value " << *itrStack << " can't be found !!" << std::endl;
//				exit(1);
//			}
//
//		}
//
//		++icount;
//
//	}
//
//#ifdef _DEBUG_WRITE
//	std::cout << "--- For the sea boundary" << std::endl;
//#endif
//
//	//--------------------------------
//	//----- For the sea boundary -----
//	//--------------------------------
//	icount = 0;
//	while( !arrayEnd.empty() ){
//
//#ifdef _DEBUG_WRITE
//		std::cout << " icount = " << icount << std::endl;
//#endif
//
//		if( icount > 1.0e6 ){
//			OutputFiles::m_logFile << "Error : icount exceed 1.0e6 in the loop !!" << std::endl;
//			exit(1);
//		}
//
//		std::vector<int> stack;
//		//std::vector<CommonParameters::XY> coords;
//		std::vector<Node> nodes;
//
//		ArrayType::iterator itrFirst = arrayEnd.begin();
//		ArrayType::iterator itrLast  = arrayEnd.end();
//
//		for( ArrayType::iterator itr = itrFirst; itr != itrLast; ){
//
//			const int vecID = itr->second;
//
//#ifdef _DEBUG_WRITE
//			std::cout << "vecID = " << vecID << std::endl;
//#endif
//
//			for( std::vector<Node>::reverse_iterator itrArray = nodeVecArray[vecID].rbegin(); itrArray != nodeVecArray[vecID].rend(); ++itrArray ){
//				nodes.push_back( *itrArray );
//			}
//			stack.push_back( vecID );
//
//#ifdef _DEBUG_WRITE
//			std::cout << "To be Erased " << stack.back() << std::endl;
//#endif
//
//			int edgeIDStart = ptrAnalysisDomain->getEdgeID( (nodeVecArray[vecID].front()).getCoordXY() );
//			const double disStart = (itr->first).second;
//
//#ifdef _DEBUG_WRITE
//			std::cout << "edgeIDStart disStart : " << edgeIDStart << " " << disStart << std::endl;
//#endif
//
//			bool goNext(false);
//
//			for( ArrayType::iterator itrNext = itrFirst; itrNext != itrLast; ++itrNext ){
//
//				const double disEndNext = (itrNext->first).first; 
//
//#ifdef _DEBUG_WRITE
//				std::cout << "disStart disEndNext : " << disStart << " " << disEndNext << std::endl;
//#endif
//
//				if( disStart < disEndNext ){
//					
//					if( std::find( stack.begin(), stack.end(), itrNext->second ) != stack.end() ){// This segment has already been found
//#ifdef _DEBUG_WRITE
//						std::cout << "Already found" << std::endl;
//#endif
//						break;
//					}
//
//					const int edgeIDEndNext = ptrAnalysisDomain->getEdgeID( (nodeVecArray[itrNext->second].back()).getCoordXY() );
//
//					while( edgeIDStart != edgeIDEndNext ){
//						nodes.push_back( Node( ptrAnalysisDomain->getEndCoord( edgeIDStart ) , true ) ); 
//						edgeIDStart = ptrAnalysisDomain->getNextEdgeIDCWR( edgeIDStart );
//#ifdef _DEBUG_WRITE
//						std::cout << "nodes.back() : " << (nodes.back()).getCoordX() << " " << (nodes.back()).getCoordY() << std::endl;
//						std::cout << "edgeIDStart : " << edgeIDStart << std::endl;
//#endif
//					}
//
//					goNext = true;
//					itr = itrNext;
//					break;
//				}
//
//			}
//
//			if( goNext ){
//				continue;
//			}
//
//#ifdef _DEBUG_WRITE
//			std::cout << "Last segment" << std::endl;
//#endif
//
//			const int edgeIDEndFirst = ptrAnalysisDomain->getEdgeID( (nodeVecArray[itrFirst->second].back()).getCoordXY() );
//
//#ifdef _DEBUG_WRITE
//			std::cout << "edgeIDEndFirst : " << edgeIDEndFirst << std::endl;
//#endif
//
//			const double disEndFirst = (itrFirst->first).first;
//			if( edgeIDEndFirst == edgeIDStart && disEndFirst < disStart ){
//				nodes.push_back( Node( ptrAnalysisDomain->getEndCoord( edgeIDStart ) , true ) ); 
//				edgeIDStart = ptrAnalysisDomain->getNextEdgeIDCWR( edgeIDStart );
//			}
//
//			while( edgeIDStart != edgeIDEndFirst ){
//				nodes.push_back( Node( ptrAnalysisDomain->getEndCoord( edgeIDStart ) , true ) ); 
//				edgeIDStart = ptrAnalysisDomain->getNextEdgeIDCWR( edgeIDStart );
//
//#ifdef _DEBUG_WRITE
//						std::cout << "nodes.back() : " << (nodes.back()).getCoordX() << " " << (nodes.back()).getCoordY() << std::endl;
//				std::cout << "edgeIDStart : " << edgeIDStart << std::endl;
//#endif
//
//			}
//			break;			
//
//		}
//
//		m_outerBoundaries.push_back( BoundaryCurve( nodes, BoundaryCurve::SEA ) );
//
//#ifdef _DEBUG_WRITE
//		std::cout << "Remove already found segment from list" << std::endl;
//#endif
//
//		// Remove already found segment from list
//#ifdef _DEBUG_WRITE
//		std::cout << "stack.size() " << stack.size() << std::endl;
//		for( std::vector<int>::iterator itr = stack.begin(); itr != stack.end(); ++itr ){
//			std::cout << *itr << std::endl;
//		}
//#endif
//
//		for( std::vector<int>::iterator itrStack = stack.begin(); itrStack != stack.end(); ++itrStack ){
//
//			bool found(false);
//
//			for( ArrayType::iterator itrArray = arrayEnd.begin(); itrArray != arrayEnd.end(); ++itrArray ){
//				if( std::find( stack.begin(), stack.end(), itrArray->second ) != stack.end() ){
//#ifdef _DEBUG_WRITE
//					std::cout << "Erase " << itrArray->second << " from arrayEnd " << std::endl;
//#endif
//					arrayEnd.erase( itrArray );
//					found = true;
//					break;
//				}
//			}
//
//			if( !found ){
//				OutputFiles::m_logFile << "Error : Element whose second value " << *itrStack << " can't be found !!" << std::endl;
//				exit(1);
//			}
//
//		}
//
//		++icount;
//
//	}
//
//}

// Make boudary curves
void BoundaryCurveList::makeBoundaryCurvesCoastLine(){
	
	OutputFiles::m_logFile << "# Make boundary curves" << std::endl;

	m_nodeList.clearList();

	// Include coast line
	const CoastLineList* const ptrCoastLineList = CoastLineList::getInstance();
	const AnalysisDomain* const ptrAnalysisDomain = AnalysisDomain::getInstance();

	const double xMin = ptrAnalysisDomain->getMinCoordX();
	const double xMax = ptrAnalysisDomain->getMaxCoordX(); 
	const double yMin = ptrAnalysisDomain->getMinCoordY(); 
	const double yMax = ptrAnalysisDomain->getMaxCoordY(); 

	const CommonParameters::XY coordUpperLeft  = { xMax, yMin };
	const CommonParameters::XY coordUpperRight = { xMax, yMax };
	const CommonParameters::XY coordLowerRight = { xMin, yMax };
	const CommonParameters::XY coordLowerLeft  = { xMin, yMin };

	const int numCoastLine = ptrCoastLineList->getNumCoastLines();
	// No coast line
	if( numCoastLine == 0 ){

		std::vector<int> nodeIDs;
		nodeIDs.push_back( m_nodeList.addNewNode( Node( coordUpperLeft,  true ) ) );
		addNodesBetweenTwoPoins( coordUpperLeft, coordUpperRight, nodeIDs );
		nodeIDs.push_back( m_nodeList.addNewNode( Node( coordUpperRight, true ) ) );
		addNodesBetweenTwoPoins( coordUpperRight, coordLowerRight, nodeIDs );
		nodeIDs.push_back( m_nodeList.addNewNode( Node( coordLowerRight, true ) ) );
		addNodesBetweenTwoPoins( coordLowerRight, coordLowerLeft, nodeIDs );
		nodeIDs.push_back( m_nodeList.addNewNode( Node( coordLowerLeft,  true ) ) );
		addNodesBetweenTwoPoins( coordLowerLeft, coordUpperLeft, nodeIDs );

		m_outerBoundaries.push_back( BoundaryCurveOuter( nodeIDs, CommonParameters::SEA ) );

		return;
	}

	//std::vector< std::vector<Node> > nodeVecArray;
	std::vector< std::vector<int> > nodeIDVecArray;

	for( int iCoastLine = 0; iCoastLine < numCoastLine; ++iCoastLine ){

		std::vector<int> nodeIDVec;

		const CoastLine* const ptrCoastLine = ptrCoastLineList->getPointerToCoastLine(iCoastLine);

		if( !ptrCoastLine->isClosed() ){// Skip if coast line is not closed
			continue;
		}

		const int numPoint = ptrCoastLine->getNumOfPoints();

		if( numPoint < 3 ){
			continue;
		}

		const Node nodeLast = ptrCoastLine->getNode(numPoint-1);
		bool locateInsidePreviously = ptrAnalysisDomain->doesLocateWithinAnalysisDomain( nodeLast.getCoordXY() ) ? true : false;
		bool hasIntersects = false;

		int segmentIDFirst = -1;

		for( int iPoint = 0; iPoint < numPoint; ++iPoint ){

			Node node = ptrCoastLine->getNode(iPoint);

			if( ptrAnalysisDomain->doesIntersectWithBoundary( node.getCoordXY() ) ){
				// This point locate on the analysis domain

				hasIntersects = true;

				//nodeVec.push_back( node );
				nodeIDVec.push_back( m_nodeList.addNewNode(node) );

				if( iPoint == 0 ){// First point
					segmentIDFirst = static_cast<int>( nodeIDVecArray.size() );
				}
				else if( iPoint == numPoint - 1 && hasIntersects ){// Last point
					if( segmentIDFirst >= 0 ){
						for( std::vector<int>::iterator itr = nodeIDVecArray[segmentIDFirst].begin(); itr != nodeIDVecArray[segmentIDFirst].end(); ++itr ){
							nodeIDVec.push_back( *itr );
						}
						nodeIDVecArray[segmentIDFirst].swap( nodeIDVec );
					}
				}

				locateInsidePreviously = true;

			}
			else if( ptrAnalysisDomain->doesLocateWithinAnalysisDomain( node.getCoordXY() ) ){
				// This point locate inside of the analysis domain
				
				if( !locateInsidePreviously ){
					nodeIDVec.push_back( m_nodeList.addNewNode( Node( Util::getCoordOfIntersectionPoint( ptrCoastLine->getCoordXYOfPoints((iPoint+numPoint-1)%numPoint), node.getCoordXY() ) , true ) ) );
					hasIntersects = true;
				}

				nodeIDVec.push_back( m_nodeList.addNewNode( node ) );

				if( iPoint == 0 ){// First point
					segmentIDFirst = static_cast<int>( nodeIDVecArray.size() );
				}
				else if( iPoint == numPoint - 1 && hasIntersects ){// Last point
					if( segmentIDFirst >= 0 ){
						for( std::vector<int>::iterator itr = nodeIDVecArray[segmentIDFirst].begin(); itr != nodeIDVecArray[segmentIDFirst].end(); ++itr ){
							nodeIDVec.push_back( *itr );
						}
						nodeIDVecArray[segmentIDFirst].swap( nodeIDVec );
					}
				}
				locateInsidePreviously = true;

			}
			else{
				// This point locate out of the analysis domain

				if( locateInsidePreviously ){
					if( iPoint == 0 ){// First point
						segmentIDFirst = static_cast<int>( nodeIDVecArray.size() );
					}
					nodeIDVec.push_back( m_nodeList.addNewNode( Node( Util::getCoordOfIntersectionPoint( ptrCoastLine->getCoordXYOfPoints((iPoint+numPoint-1)%numPoint), node.getCoordXY() ) , true ) ) );
					nodeIDVecArray.push_back( nodeIDVec );
					nodeIDVec.clear();
					hasIntersects = true;
				}

				locateInsidePreviously = false;

			}
		}

		if( !hasIntersects && static_cast<int>( nodeIDVec.size() ) >= 3 ){
			if(isClockWise( nodeIDVec )){
				reverse( nodeIDVec.begin(), nodeIDVec.end() );
				m_innerBoundaries.push_back( BoundaryCurveInner( nodeIDVec, CommonParameters::LAND ) );
			}
			else{
				m_innerBoundaries.push_back( BoundaryCurveInner( nodeIDVec, CommonParameters::LAKE ) );
			}
		}

	}

	if( nodeIDVecArray.empty() ){
		std::vector<int> nodeIDs;

		nodeIDs.push_back( m_nodeList.addNewNode( Node( coordUpperLeft,  true ) ) );
		addNodesBetweenTwoPoins( coordUpperLeft, coordUpperRight, nodeIDs );
		nodeIDs.push_back( m_nodeList.addNewNode( Node( coordUpperRight, true ) ) );
		addNodesBetweenTwoPoins( coordUpperRight, coordLowerRight, nodeIDs );
		nodeIDs.push_back( m_nodeList.addNewNode( Node( coordLowerRight, true ) ) );
		addNodesBetweenTwoPoins( coordLowerRight, coordLowerLeft, nodeIDs );
		nodeIDs.push_back( m_nodeList.addNewNode( Node( coordLowerLeft,  true ) ) );
		addNodesBetweenTwoPoins( coordLowerLeft, coordUpperLeft, nodeIDs );

		if( m_innerBoundaries.empty() ){
			m_outerBoundaries.push_back( BoundaryCurveOuter( nodeIDs, CommonParameters::LAND) );
		}else{
			m_outerBoundaries.push_back( BoundaryCurveOuter( nodeIDs, CommonParameters::SEA) );
		}

		return;
	}

	typedef std::vector< std::pair< std::pair<double, double>, int > > ArrayType;
	ArrayType arrayStart;
	ArrayType arrayEnd;

	int icount(0);
	for( std::vector< std::vector<int> >::iterator itr = nodeIDVecArray.begin(); itr != nodeIDVecArray.end(); ++itr, ++icount ){
		const int startID = itr->front();
		const int endID = itr->back();

		const double disStart = ptrAnalysisDomain->calcDistanceOnBoundaryFromUpperLeft( m_nodeList.getCoordXYOfPoints(startID) );
		
		const double disEnd = ptrAnalysisDomain->calcDistanceOnBoundaryFromUpperLeft( m_nodeList.getCoordXYOfPoints(endID) );
		arrayStart.push_back( std::make_pair( std::make_pair( disStart , disEnd ) , icount ) );
		arrayEnd.push_back( std::make_pair( std::make_pair( disEnd , disStart ) , icount ) );
	}

	std::sort( arrayStart.begin(), arrayStart.end() );// Ascending order
	std::sort( arrayEnd.begin(), arrayEnd.end() );// // Ascending order

#ifdef _DEBUG_WRITE
	for( ArrayType::iterator itr = arrayStart.begin(); itr != arrayStart.end(); ++itr ){
		std::cout << "arrayStart first second : " << (itr->first).first << " " << (itr->first).second << " " << itr->second << std::endl;
	}
	for( ArrayType::iterator itr = arrayEnd.begin(); itr != arrayEnd.end(); ++itr ){
		std::cout << "arrayEnd first second : " << (itr->first).first << " " << (itr->first).second << " " << itr->second << std::endl;
	}
#endif

#ifdef _DEBUG_WRITE
	std::cout << "--- For the land boundary" << std::endl;
#endif
	
	//m_nodeList.writeNode2DListToVTK( "temp.vtk" );
	
	//---------------------------------
	//----- For the land boundary -----
	//---------------------------------
	icount = 0;
	while( !arrayStart.empty() ){

		if( icount > 1.0e6 ){
			OutputFiles::m_logFile << "Error : icount exceed 1.0e6 in the loop !!" << std::endl;
			exit(1);
		}

		std::vector<int> stack;
		std::vector<int> nodeIDs;

		ArrayType::iterator itrFirst = arrayStart.begin();
		ArrayType::iterator itrLast  = arrayStart.end();

		for( ArrayType::iterator itr = itrFirst; itr != itrLast; ){

			const int vecID = itr->second;

			nodeIDs.insert( nodeIDs.end(), nodeIDVecArray[vecID].begin(), nodeIDVecArray[vecID].end() );
			stack.push_back( vecID );

			const CommonParameters::XY coordEnd = m_nodeList.getCoordXYOfPoints( nodeIDVecArray[vecID].back() );
			int edgeIDEnd = ptrAnalysisDomain->getEdgeID( coordEnd );
			const double disEnd = (itr->first).second;

			bool goNext(false);

			for( ArrayType::iterator itrNext = itrFirst; itrNext != itrLast; ++itrNext ){

				const double disStartNext = (itrNext->first).first; 

#ifdef _DEBUG_WRITE
				std::cout << "disEnd : " << disEnd << std::endl;
				std::cout << "disStartNext : " << disStartNext << std::endl;
#endif

				if( disEnd < disStartNext ){

					if( std::find( stack.begin(), stack.end(), itrNext->second ) != stack.end() ){// This segment has already been found
						break;
					}

					const CommonParameters::XY coordStartNext = m_nodeList.getCoordXYOfPoints( nodeIDVecArray[itrNext->second].front() );
					const int edgeIDStartNext = ptrAnalysisDomain->getEdgeID( coordStartNext );

					CommonParameters::XY coord = coordEnd;
					while( edgeIDEnd != edgeIDStartNext ){
						const CommonParameters::XY coordNext = ptrAnalysisDomain->getEndCoord( edgeIDEnd );
						addNodesBetweenTwoPoins( coord, coordNext, nodeIDs );
						nodeIDs.push_back( m_nodeList.addNewNode( Node( coordNext, true ) ) ); 
						edgeIDEnd = ptrAnalysisDomain->getNextEdgeIDCWR( edgeIDEnd );
						coord = coordNext;
					}
					addNodesBetweenTwoPoins( coord, coordStartNext, nodeIDs );

					goNext = true;
					itr = itrNext;
					break;
				}

			}

			if( goNext ){
				continue;
			}

			const CommonParameters::XY coordFirst = m_nodeList.getCoordXYOfPoints( nodeIDVecArray[itrFirst->second].front() );
			const int edgeIDStartFirst = ptrAnalysisDomain->getEdgeID( coordFirst );

			const double disStartFirst = (itrFirst->first).first;

			CommonParameters::XY coord = coordEnd;

			if( edgeIDStartFirst == edgeIDEnd && disStartFirst < disEnd ){
				const CommonParameters::XY coordNext = ptrAnalysisDomain->getEndCoord( edgeIDEnd );
				addNodesBetweenTwoPoins( coord, coordNext, nodeIDs );
				nodeIDs.push_back( m_nodeList.addNewNode( Node( coordNext, true ) ) ); 
				edgeIDEnd = ptrAnalysisDomain->getNextEdgeIDCWR( edgeIDEnd );
				coord = coordNext;
			}

			while( edgeIDEnd != edgeIDStartFirst ){
				const CommonParameters::XY coordNext = ptrAnalysisDomain->getEndCoord( edgeIDEnd );
				addNodesBetweenTwoPoins( coord, coordNext, nodeIDs );
				nodeIDs.push_back( m_nodeList.addNewNode( Node( coordNext, true ) ) ); 
				edgeIDEnd = ptrAnalysisDomain->getNextEdgeIDCWR( edgeIDEnd );
				coord = coordNext;
			}

			addNodesBetweenTwoPoins( coord, coordFirst, nodeIDs );

			break;			

		}

		m_outerBoundaries.push_back( BoundaryCurveOuter( nodeIDs, CommonParameters::LAND ) );

		// Remove already found segment from list
		for( std::vector<int>::iterator itrStack = stack.begin(); itrStack != stack.end(); ++itrStack ){

			bool found(false);

			for( ArrayType::iterator itrArray = arrayStart.begin(); itrArray != arrayStart.end(); ++itrArray ){
				if( std::find( stack.begin(), stack.end(), itrArray->second ) != stack.end() ){
					arrayStart.erase( itrArray );
					found = true;
					break;
				}
			}

			if( !found ){
				OutputFiles::m_logFile << "Error : Element whose second value " << *itrStack << " can't be found !!" << std::endl;
				exit(1);
			}

		}

		++icount;

	}

#ifdef _DEBUG_WRITE
	std::cout << "--- For the sea boundary" << std::endl;
#endif

	//--------------------------------
	//----- For the sea boundary -----
	//--------------------------------
	icount = 0;
	while( !arrayEnd.empty() ){

		if( icount > 1.0e6 ){
			OutputFiles::m_logFile << "Error : icount exceed 1.0e6 in the loop !!" << std::endl;
			exit(1);
		}

		std::vector<int> stack;
		std::vector<int> nodeIDs;

		ArrayType::iterator itrFirst = arrayEnd.begin();
		ArrayType::iterator itrLast  = arrayEnd.end();

		for( ArrayType::iterator itr = itrFirst; itr != itrLast; ){

			const int vecID = itr->second;

			for( std::vector<int>::reverse_iterator itrArray = nodeIDVecArray[vecID].rbegin(); itrArray != nodeIDVecArray[vecID].rend(); ++itrArray ){
				nodeIDs.push_back( *itrArray );
			}
			stack.push_back( vecID );
			
			const CommonParameters::XY coordStart = m_nodeList.getCoordXYOfPoints( nodeIDVecArray[vecID].front() );
			int edgeIDStart = ptrAnalysisDomain->getEdgeID( coordStart );
			const double disStart = (itr->first).second;

			bool goNext(false);

			for( ArrayType::iterator itrNext = itrFirst; itrNext != itrLast; ++itrNext ){

				const double disEndNext = (itrNext->first).first; 

				if( disStart < disEndNext ){

					if( std::find( stack.begin(), stack.end(), itrNext->second ) != stack.end() ){// This segment has already been found
						break;
					}
					
					const CommonParameters::XY coordEndNext = m_nodeList.getCoordXYOfPoints( nodeIDVecArray[itrNext->second].back() );
					const int edgeIDEndNext = ptrAnalysisDomain->getEdgeID( coordEndNext );

					CommonParameters::XY coord = coordStart;
					while( edgeIDStart != edgeIDEndNext ){
						const CommonParameters::XY coordNext = ptrAnalysisDomain->getEndCoord( edgeIDStart );
						addNodesBetweenTwoPoins( coord, coordNext, nodeIDs );
						nodeIDs.push_back( m_nodeList.addNewNode( Node( coordNext, true ) ) ); 
						edgeIDStart = ptrAnalysisDomain->getNextEdgeIDCWR( edgeIDStart );
						coord = coordNext;
					}

					addNodesBetweenTwoPoins( coord, coordEndNext, nodeIDs );

					goNext = true;
					itr = itrNext;
					break;
				}

			}

			if( goNext ){
				continue;
			}
			
			const CommonParameters::XY coordFirst = m_nodeList.getCoordXYOfPoints( nodeIDVecArray[itrFirst->second].back() );
			const int edgeIDEndFirst = ptrAnalysisDomain->getEdgeID( coordFirst );

			const double disEndFirst = (itrFirst->first).first;

			CommonParameters::XY coord = coordStart;

			if( edgeIDEndFirst == edgeIDStart && disEndFirst < disStart ){
				const CommonParameters::XY coordNext = ptrAnalysisDomain->getEndCoord( edgeIDStart );
				addNodesBetweenTwoPoins( coord, coordNext, nodeIDs );
				nodeIDs.push_back( m_nodeList.addNewNode( Node( coordNext, true ) ) ); 
				edgeIDStart = ptrAnalysisDomain->getNextEdgeIDCWR( edgeIDStart );
				coord = coordNext;
			}

			while( edgeIDStart != edgeIDEndFirst ){
				const CommonParameters::XY coordNext = ptrAnalysisDomain->getEndCoord( edgeIDStart );
				addNodesBetweenTwoPoins( coord, coordNext, nodeIDs );
				nodeIDs.push_back( m_nodeList.addNewNode( Node( coordNext, true ) ) ); 
				edgeIDStart = ptrAnalysisDomain->getNextEdgeIDCWR( edgeIDStart );
				coord = coordNext;
			}

			addNodesBetweenTwoPoins( coord, coordFirst, nodeIDs );

			break;			

		}

		m_outerBoundaries.push_back( BoundaryCurveOuter( nodeIDs, CommonParameters::SEA ) );

		// Remove already found segment from list
		for( std::vector<int>::iterator itrStack = stack.begin(); itrStack != stack.end(); ++itrStack ){

			bool found(false);

			for( ArrayType::iterator itrArray = arrayEnd.begin(); itrArray != arrayEnd.end(); ++itrArray ){
				if( std::find( stack.begin(), stack.end(), itrArray->second ) != stack.end() ){
					arrayEnd.erase( itrArray );
					found = true;
					break;
				}
			}

			if( !found ){
				OutputFiles::m_logFile << "Error : Element whose second value " << *itrStack << " can't be found !!" << std::endl;
				exit(1);
			}

		}

		++icount;

	}


}

#ifdef _MOD_FOR_NMT
// Add NMT dipoles to boudary curve list
void BoundaryCurveList::addNMTDipolesToBoundaryCurveList(){

	//int nodeIDFirst = m_nodeList.getTotalNumberOfNode();

	const  ObservingSiteList* ptrObservingSiteList = ObservingSiteList::getInstance();
	const int numObsLine = ptrObservingSiteList->getNumObsLine();

	for( int iobs = 0; iobs < numObsLine; ++iobs ){
		std::vector<int> nodeIDs;
		const CommonParameters::XY startCoord = ptrObservingSiteList->getCoordOfStartPointObsLine(iobs);
		const CommonParameters::XY endCoord = ptrObservingSiteList->getCoordOfEndPointObsLine(iobs);

		nodeIDs.push_back( m_nodeList.addNewNode( Node( startCoord,  true ) ) );
		addNodesBetweenTwoPoins( startCoord, endCoord, nodeIDs );
		nodeIDs.push_back( m_nodeList.addNewNode( Node( endCoord, true ) ) );

		const int numNode = static_cast<int>( nodeIDs.size() );
		for( int iNode = 1; iNode < numNode; ++iNode ){
			const int nodeID0 = nodeIDs[iNode - 1];
			const int nodeID1 = nodeIDs[iNode];
			addNodesBetweenTwoPoins( m_nodeList.getCoordXYOfPoints(nodeID0), m_nodeList.getCoordXYOfPoints(nodeID1), nodeIDs );
		}

		const int numNodeMod = static_cast<int>( nodeIDs.size() );
		std::map<double, int> distanceToNodeID;
		for( int iNode = 1; iNode < numNodeMod; ++iNode ){
			const int nodeID = nodeIDs[iNode];
			const double distance = Util::calcEdgeLength( startCoord, m_nodeList.getCoordXYOfPoints(nodeID) );
			distanceToNodeID.insert( std::make_pair( distance, nodeID ) );
		}
		int icount = 1;
		for( std::map<double, int>::const_iterator itr = distanceToNodeID.begin(); itr != distanceToNodeID.end(); ++itr, ++icount ){
			nodeIDs[icount] = itr->second;
		}

		m_innerBoundaries.push_back( BoundaryCurveInner( nodeIDs, CommonParameters::NMT_DIPOLE) );
	}
	
}
#endif


// Relate inner boundaries to outer boundaries containing them
// [note] This function does not support sub-inner boundary
void BoundaryCurveList::relateInnerBoundToOuterBound(){

	OutputFiles::m_logFile << "# Relate inner boundaries to outer boundaries" << std::endl;

	const int numOuterBoundary = static_cast<int>( m_outerBoundaries.size() );

	m_outer2InnerBoundaries = new std::vector<int>[numOuterBoundary];

	int icountIn(0);
	for( std::vector<BoundaryCurveInner>::iterator itrIn = m_innerBoundaries.begin(); itrIn != m_innerBoundaries.end(); ++itrIn, ++icountIn ){// Loop of inner boundary
		bool found(false);
		int icountOut(0);
		for( std::vector<BoundaryCurveOuter>::iterator itrOut = m_outerBoundaries.begin(); itrOut != m_outerBoundaries.end(); ++itrOut, ++icountOut ){// Loop of outer boundary
			if( itrOut->BoundaryCurve::include( *itrIn, &m_nodeList ) ){
				if( itrIn->getGeologicalType() == itrOut->getGeologicalType() ){
					OutputFiles::m_logFile << "Error : Geological types of the inner boundary " << icountIn << " is the same as the one of the outer boundary " << icountOut << " !!" << std::endl;
					exit(1);
				}
				// The outer boundary include the inner boundary
				m_outer2InnerBoundaries[icountOut].push_back(icountIn);
				found = true;
				break;
			}
		}
		if( !found ){
			OutputFiles::m_logFile << "Error : Inner boundary " << icountIn << " is not included by any outer boundary !!" << std::endl;
			exit(1);
		}
	}


#ifdef _DEBUG_WRITE
	//for( std::multimap<int,int>::iterator itr = m_outer2InnerBoundaries.begin(); itr != m_outer2InnerBoundaries.end(); ++itr ){
	//	std::cout << itr->first << " " << itr->second << std::endl;
	//}
	for( int iBoun = 0; iBoun < numOuterBoundary; ++iBoun ){
		for( std::vector<int>::iterator	itr = m_outer2InnerBoundaries[iBoun].begin();
			itr != m_outer2InnerBoundaries[iBoun].end(); ++itr ){
			std::cout << iBoun << " " << *itr << std::endl;
		}
	}

	std::cout << "End" << std::endl;
#endif

};

// Relate sb-inner boundaries to inner boundaries containing them
// [note] Implementation has not been completed
void BoundaryCurveList::relateSubInnerBoundToInnerBound(){

	const int numInnerBoundary = static_cast<int>( m_innerBoundaries.size() );

	m_inner2SubInnerBoundaries = new std::vector<int>[numInnerBoundary];


}

//----- Don't delete for future use >>>>>
//// Add boundary curves of anomalies to list
//// [note] Change future to the function adding boundary curves of lakes
//void BoundaryCurveList::addBoudaryCurveOfAnomalies(){
//
////#ifdef _DEBUG_WRITE
////	std::cout << "m_outer2InnerBoundaries" << std::endl;
////	for( std::multimap<int,int>::iterator itr = m_outer2InnerBoundaries.begin(); itr != m_outer2InnerBoundaries.end(); ++itr ){
////		std::cout << itr->first << " " << itr->second << std::endl;
////	}
////	std::cout << "m_inner2SubInnerBoundaries" << std::endl;
////	for( std::multimap<int,int>::iterator itr = m_inner2SubInnerBoundaries.begin(); itr != m_inner2SubInnerBoundaries.end(); ++itr ){
////		std::cout << itr->first << " " << itr->second << std::endl;
////	}
////#endif
//
//	const int numInnerBoundary = static_cast<int>( m_innerBoundaries.size() );
//	m_inner2SubInnerBoundaries = new std::vector<int>[numInnerBoundary];
//
//	const AnomalyList* const ptrAnomalyList = AnomalyList::getInstance();
//	const int numAnomaly = ptrAnomalyList->getNumOfAnomalies();
//
//	for( int iano = 0; iano < numAnomaly; ++iano ){
//
//		const Anomaly* const ptrAnomaly = ptrAnomalyList->getPtrAnomaly( iano );
//		std::vector<CommonParameters::XY> coordVec;
//		ptrAnomaly->calcCoordOfCrossLine( coordVec );
//
//		CommonParameters::XY centerCoord = ptrAnomaly->getCenterCoordXY();
//
//		bool found(false);
//
//		int icountIn(0);
//		for( std::vector<BoundaryCurveInner>::iterator itrIn = m_innerBoundaries.begin(); itrIn != m_innerBoundaries.end(); ++itrIn ){// Loop of inner boundary
//
//			if( itrIn->BoundaryCurve::include( centerCoord, &m_nodeList ) ){
//				// The outer boundary include the anomaly
//
//				std::vector<int> nodeIDs;
//				for( std::vector<CommonParameters::XY>::iterator itr = coordVec.begin(); itr != coordVec.end(); ++itr ){
//					nodeIDs.push_back( m_nodeList.addNewNode( Node( *itr,  true ) ) );
//				}
//
//				m_subInnerBoundaries.push_back( BoundaryCurveSubInner( nodeIDs, BoundaryCurve::LAKE ) );
//
//				m_inner2SubInnerBoundaries[icountIn].push_back( static_cast<int>(m_subInnerBoundaries.size()) - 1 );
//
//				found = true;
//				break;
//			}
//
//			++icountIn;
//		}
//
//		if( found ){
//			continue;
//		}
//
//		int icountOut(0);
//		for( std::vector<BoundaryCurveOuter>::iterator itrOut = m_outerBoundaries.begin(); itrOut != m_outerBoundaries.end(); ++itrOut ){// Loop of outer boundary
//
//			if( itrOut->BoundaryCurve::include( centerCoord, &m_nodeList ) ){
//				// The outer boundary include the anomaly
//
//				std::vector<int> nodeIDs;
//				for( std::vector<CommonParameters::XY>::iterator itr = coordVec.begin(); itr != coordVec.end(); ++itr ){
//					nodeIDs.push_back( m_nodeList.addNewNode( Node( *itr,  true ) ) );
//				}
//
//				m_innerBoundaries.push_back( BoundaryCurveInner( nodeIDs, BoundaryCurve::LAKE ) );
//
//				//m_outer2innerBoundaries[icountOut].push_back( icountOut );
//				m_outer2InnerBoundaries[icountOut].push_back( static_cast<int>(m_innerBoundaries.size()) - 1 );
//
//				found = true;
//				break;
//			}
//
//			++icountOut;
//		}
//
//		if( !found ){
//			OutputFiles::m_logFile << "Error : Anomaly " << iano << " is not included by the boudary of the coast lines !!" << std::endl;
//			exit(1);
//		}
//
//	}
//
////#ifdef _DEBUG_WRITE
////	std::cout << "m_outer2InnerBoundaries" << std::endl;
////	for( std::multimap<int,int>::iterator itr = m_outer2InnerBoundaries.begin(); itr != m_outer2InnerBoundaries.end(); ++itr ){
////		std::cout << itr->first << " " << itr->second << std::endl;
////	}
////	std::cout << "m_inner2SubInnerBoundaries" << std::endl;
////	for( std::multimap<int,int>::iterator itr = m_inner2SubInnerBoundaries.begin(); itr != m_inner2SubInnerBoundaries.end(); ++itr ){
////		std::cout << itr->first << " " << itr->second << std::endl;
////	}
////#endif
//
//	
//}
//----- Don't delete for future use <<<<<

//// Relate node to boundary curve
//void BoundaryCurveList::relateNodeToBoundaryCurve(){
//
//	//NodeList* ptrNode2DList = NodeList::getInstance();
//
//	int icount(0);
//	for( std::vector<BoundaryCurveOuter>::const_iterator itr = m_outerBoundaries.begin(); itr != m_outerBoundaries.end(); ++itr ){
//		const int numPointOfSegment = itr->getNumOfPoints();
//		for( int i = 0; i < numPointOfSegment; ++i ){
//			//(itr->getPointerToNode(i))->addBoundaryCurveID( BoundaryCurve::OUTER, icount );
//			m_nodeList.getPointerToNode( itr->getNodeID(i) )->addBoundaryCurveID( BoundaryCurve::OUTER, icount );
//		}
//		++icount;
//	}
//
//	icount = 0;
//	for( std::vector<BoundaryCurveInner>::const_iterator itr = m_innerBoundaries.begin(); itr != m_innerBoundaries.end(); ++itr ){
//		const int numPointOfSegment = itr->getNumOfPoints();
//		for( int i = 0; i < numPointOfSegment; ++i ){
//			//(itr->getPointerToNode(i))->addBoundaryCurveID( BoundaryCurve::INNER, icount );
//			m_nodeList.getPointerToNode( itr->getNodeID(i) )->addBoundaryCurveID( BoundaryCurve::INNER, icount );
//		}
//		++icount;
//	}
//
//	icount = 0;
//	for( std::vector<BoundaryCurveSubInner>::const_iterator itr = m_subInnerBoundaries.begin(); itr != m_subInnerBoundaries.end(); ++itr ){
//		const int numPointOfSegment = itr->getNumOfPoints();
//		for( int i = 0; i < numPointOfSegment; ++i ){
//			//(itr->getPointerToNode(i))->addBoundaryCurveID( BoundaryCurve::SUB_INNER, icount );
//			m_nodeList.getPointerToNode( itr->getNodeID(i) )->addBoundaryCurveID( BoundaryCurve::SUB_INNER, icount );
//		}
//		++icount;
//	}
//
//}

// Write boundary curve list
void BoundaryCurveList::writeBoundaryCurveList( const std::string& fileName ) const{

	std::ofstream ofs( fileName.c_str() ); 

	if( !ofs.is_open() ) {
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << " !!" << std::endl;
		exit(1);
	}

	ofs << "*OuterBoundary" << std::endl;
	const int numOuterBoundary = static_cast<int>( m_outerBoundaries.size() );
	ofs << std::setw(10) << numOuterBoundary << std::endl;
	int icount(0);
	for( std::vector<BoundaryCurveOuter>::const_iterator itr = m_outerBoundaries.begin(); itr != m_outerBoundaries.end(); ++itr ){
		//ofs << "*ID " << std::setw(6) << icount++ << std::endl;
		itr->writeBoudaryCurve( ofs );
	}

	ofs << "*InnerBoundary" << std::endl;
	const int numInnerBoundary = static_cast<int>( m_innerBoundaries.size() );
	ofs << std::setw(10) << numInnerBoundary << std::endl;
	icount = 0;
	for( std::vector<BoundaryCurveInner>::const_iterator itr = m_innerBoundaries.begin(); itr != m_innerBoundaries.end(); ++itr ){
		//ofs << "*ID " << std::setw(6) << icount++ << std::endl;
		itr->writeBoudaryCurve( ofs );
	}

	ofs << "*Sub-innerBoundary" << std::endl;
	const int numSubInnerBoundary = static_cast<int>( m_subInnerBoundaries.size() );
	ofs << std::setw(10) << numSubInnerBoundary << std::endl;
	icount = 0;
	for( std::vector<BoundaryCurveSubInner>::const_iterator itr = m_subInnerBoundaries.begin(); itr != m_subInnerBoundaries.end(); ++itr ){
		//ofs << "*ID " << std::setw(6) << icount++ << std::endl;
		itr->writeBoudaryCurve( ofs );
	}
	
	ofs.close();

}

// Write containment relationship of boundary curves
void BoundaryCurveList::writeBoundaryCurveRelatios( const std::string& fileName ) const{

	std::ofstream ofs( fileName.c_str() ); 

	if( !ofs.is_open() ) {
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << " !!" << std::endl;
		exit(1);
	}

	//ofs << std::setw(10) << static_cast<int>( m_outer2InnerBoundaries.size() ) << std::endl;
	//for( std::multimap<int,int>::const_iterator itr = m_outer2InnerBoundaries.begin(); itr != m_outer2InnerBoundaries.end(); ++itr ){
	//	ofs << std::setw(10) << itr->first << std::setw(10) << itr->second << std::endl;
	//}

	//ofs << std::setw(10) << static_cast<int>( m_inner2SubInnerBoundaries.size() ) << std::endl;
	//for( std::multimap<int,int>::const_iterator itr = m_inner2SubInnerBoundaries.begin(); itr != m_inner2SubInnerBoundaries.end(); ++itr ){
	//	ofs << std::setw(10) << itr->first << std::setw(10) << itr->second << std::endl;
	//}

	const int numBounOuter = static_cast<int>( m_outerBoundaries.size() );
	ofs << std::setw(10) << numBounOuter << std::endl;
	for( int iBoun = 0; iBoun < numBounOuter; ++iBoun ){
		ofs << std::setw(10) << m_outer2InnerBoundaries[iBoun].size();
		if( m_outer2InnerBoundaries[iBoun].empty() ){
			ofs << std::endl;
			continue;
		}
		for( std::vector<int>::const_iterator itr = m_outer2InnerBoundaries[iBoun].begin(); itr != m_outer2InnerBoundaries[iBoun].end(); ++itr ){
			ofs << std::setw(10) << *itr;
		}
		ofs << std::endl;
	}

	const int numBounInner = static_cast<int>( m_innerBoundaries.size() );

	ofs << std::setw(10) << numBounInner << std::endl;
	for( int iBoun = 0; iBoun < numBounInner; ++iBoun ){
		ofs << std::setw(10) << m_inner2SubInnerBoundaries[iBoun].size();
		if( m_inner2SubInnerBoundaries[iBoun].empty() ){
			ofs << std::endl;
			continue;
		}
		for( std::vector<int>::const_iterator itr = m_inner2SubInnerBoundaries[iBoun].begin(); itr != m_inner2SubInnerBoundaries[iBoun].end(); ++itr ){
			ofs << std::setw(10) << *itr;
		}
		ofs << std::endl;
	}

	ofs.close();

}

// Read boundary curve list
void BoundaryCurveList::readBoundaryCurveList( const std::string& fileName ){

	std::ifstream ifs( fileName.c_str() ); 

	if( !ifs.is_open() ){
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << " !!" << std::endl;
		exit(1);
	}

	OutputFiles::m_logFile << "# Read boundary curve list from " << fileName.c_str() << std::endl;

	std::string sbuf;

	// Outer boundary
	ifs >> sbuf;

#ifdef _DEBUG_WRITE
	std::cout << " sbuf = " << sbuf << std::endl;
#endif

	int numOuterBoundary(0);
	ifs >> numOuterBoundary;

	for( int iBOut = 0; iBOut < numOuterBoundary; ++iBOut ){
		m_outerBoundaries.push_back( BoundaryCurveOuter() );
		(m_outerBoundaries.back()).readBoudaryCurve( ifs );
	}

	// Inner boundary
	ifs >> sbuf;

#ifdef _DEBUG_WRITE
	std::cout << " sbuf = " << sbuf << std::endl;
#endif

	int numInnerBoundary(0);
	ifs >> numInnerBoundary;

	for( int iBIn = 0; iBIn < numInnerBoundary; ++iBIn ){
		m_innerBoundaries.push_back( BoundaryCurveInner() );
		(m_innerBoundaries.back()).readBoudaryCurve( ifs );
	}

	// Sub-inner boundary
	ifs >> sbuf;

#ifdef _DEBUG_WRITE
	std::cout << " sbuf = " << sbuf << std::endl;
#endif

	int numSubInnerBoundary(0);
	ifs >> numSubInnerBoundary;

	for( int iBSubIn = 0; iBSubIn < numSubInnerBoundary; ++iBSubIn ){
		m_subInnerBoundaries.push_back( BoundaryCurveSubInner() );
		(m_subInnerBoundaries.back()).readBoudaryCurve( ifs );
	}

	ifs.close();

}

// Read containment relationship of boundary curves 
void BoundaryCurveList::readBoundaryCurveRelatios( const std::string& fileName ){

	std::ifstream ifs( fileName.c_str() ); 

	if( !ifs.is_open() ){
		OutputFiles::m_logFile << " Error : Cannot open file " << fileName.c_str() << " !!" << std::endl;
		exit(1);
	}

	OutputFiles::m_logFile << "# Read containment relationship of boundary curves from " << fileName.c_str() << std::endl;

	//int numOuterToInner(0);
	//ifs >> numOuterToInner;

	//for( int imap = 0; imap < numOuterToInner; ++imap ){
	//	std::pair<int,int> bondPair;
	//	ifs >> bondPair.first >> bondPair.second; 
	//	m_outer2InnerBoundaries.insert( bondPair );
	//}

	//int numInnerToSub(0);
	//ifs >> numInnerToSub;

	//for( int imap = 0; imap < numInnerToSub; ++imap ){
	//	std::pair<int,int> bondPair;
	//	ifs >> bondPair.first >> bondPair.second; 
	//	m_inner2SubInnerBoundaries.insert( bondPair );
	//}

	int numBounOuter(0);
	ifs >> numBounOuter;
	m_outer2InnerBoundaries = new std::vector<int>[numBounOuter];
	for( int iBoun = 0; iBoun < numBounOuter; ++iBoun ){
		int nInner(0);
		ifs >> nInner;
		for( int i = 0; i < nInner; ++i ){
			int iInner(0);
			ifs >> iInner;
			m_outer2InnerBoundaries[iBoun].push_back(iInner);
		}
	}

	int numBounInner(0);
	ifs >> numBounInner;
	m_inner2SubInnerBoundaries = new std::vector<int>[numBounInner];
	for( int iBoun = 0; iBoun < numBounInner; ++iBoun ){
		int nSubInner(0);
		ifs >> nSubInner;
		for( int i = 0; i < nSubInner; ++i ){
			int iSubInner(0);
			ifs >> iSubInner;
			m_inner2SubInnerBoundaries[iBoun].push_back(iSubInner);
		}
	}

	ifs.close();
	
}

// Set one outer boundary curve
void BoundaryCurveList::setOneOuterBoundaryCurve( const int numNodes, const int* nodeIDs, const int geologicalType ){
	m_outerBoundaries.push_back( BoundaryCurveOuter(geologicalType) );
	(m_outerBoundaries.back()).setOneBoundaryCurve( numNodes, nodeIDs );
	m_outer2InnerBoundaries = new std::vector<int>[1];
}

// Make all node fixed
void BoundaryCurveList::makeAllNodeFixed(){
	
	OutputFiles::m_logFile << "# Make all node fixed" << std::endl;

	m_nodeList.makeAllNodeFixed();

}

// Make node data to VTK file
void BoundaryCurveList::writeNode2DListToVTK( const std::string& fileName ) const{

	m_nodeList.writeNode2DListToVTK( fileName );

}

// Write node 2D data to the Intermediate file of step 1
void BoundaryCurveList::writeNode2DListToIntermediateFileStep1( const std::string& fileName ) const{

	m_nodeList.writeNode2DListToIntermediateFileStep1( fileName );

}

// Read node list of 2D mesh from intermediate file after step 1
void BoundaryCurveList::readNode2DListFromIntermediateFile( const std::string& fileName ){

	m_nodeList.readNode2DListFromIntermediateFile( fileName );

}

// Set node list of 2D mesh
void BoundaryCurveList::setNode2DList( const int numNodes, const CommonParameters::XY* coord2D ){
	
	m_nodeList.setNode2DList( numNodes, coord2D );

}

// Get pointer to the list of the node 2D
const NodeList* BoundaryCurveList::getPointerToNodeList() const{

	return &m_nodeList;

}

// Get total number of the outer boundaries
int BoundaryCurveList::getTotalNumberOuterBoundaries() const{
	return static_cast<int>( m_outerBoundaries.size() );
}

// Get total number of the inner boundaries
int BoundaryCurveList::getTotalNumberInnerBoundaries() const{
	return static_cast<int>( m_innerBoundaries.size() );
}

// Get total number of the sub-inner boundaries
int BoundaryCurveList::getTotalNumberSubInnerBoundaries() const{
	return static_cast<int>( m_subInnerBoundaries.size() );
}

// Get pointer to the outer boundaries
const BoundaryCurveOuter* BoundaryCurveList::getPointerToOuterBoundary( const int id ) const{

	if( id < 0 || id >= getTotalNumberOuterBoundaries() ){
		OutputFiles::m_logFile << "Error : ID of boudary is out of range !! ID = " << id << std::endl;
		exit(1);
	}

	return &m_outerBoundaries[id];
}

// Get pointer to the inner boundaries
const BoundaryCurveInner* BoundaryCurveList::getPointerToInnerBoundary( const int id ) const{

	if( id < 0 || id >= getTotalNumberInnerBoundaries() ){
		OutputFiles::m_logFile << "Error : ID of boudary is out of range !! ID = " << id << std::endl;
		exit(1);
	}

	return &m_innerBoundaries[id];
}

// Get pointer to the sub-inner boundaries
const BoundaryCurveSubInner* BoundaryCurveList::getPointerToSubInnerBoundary( const int id ) const{

	if( id < 0 || id >= getTotalNumberSubInnerBoundaries() ){
		OutputFiles::m_logFile << "Error : ID of boudary is out of range !! ID = " << id << std::endl;
		exit(1);
	}

	return &m_subInnerBoundaries[id];
}

// Get total number of nodes belong to the outer boundaries
int BoundaryCurveList::getNumNodeOuterBoundary( const int iBoun ) const{

	return m_outerBoundaries[iBoun].getNumOfPoints();

}

// Get total number of nodes belong to the inner boundaries
int BoundaryCurveList::getNumNodeInnerBoundary( const int iBoun ) const{

	return m_innerBoundaries[iBoun].getNumOfPoints();

}

// Get total number of nodes belong to the sub-inner boundaries
int BoundaryCurveList::getNumNodeSubInnerBoundary( const int iBoun ) const{

	return m_subInnerBoundaries[iBoun].getNumOfPoints();

}

// Get coordinate of points belong to the outer boundaries
CommonParameters::XY BoundaryCurveList::getPointCoordOuterBoundary( const int iBoun, const int iNode ) const{

	if( iBoun < 0 || iBoun >= getTotalNumberOuterBoundaries() ){
		OutputFiles::m_logFile << "Error : ID of boudary is out of range !! ID = " << iBoun << std::endl;
		exit(1);
	}

	return m_outerBoundaries[iBoun].getCoordXYOfPoints( iNode, &m_nodeList );

}

// Get coordinate of points belong to the inner boundaries
CommonParameters::XY BoundaryCurveList::getPointCoordInnerBoundary( const int iBoun, const int iNode ) const{

	if( iBoun < 0 || iBoun >= getTotalNumberInnerBoundaries() ){
		OutputFiles::m_logFile << "Error : ID of boudary is out of range !! ID = " << iBoun << std::endl;
		exit(1);
	}

	return m_innerBoundaries[iBoun].getCoordXYOfPoints( iNode, &m_nodeList );

}

// Get coordinate of points belong to the sub-inner boundaries
CommonParameters::XY BoundaryCurveList::getPointCoordSubInnerBoundary( const int iBoun, const int iNode ) const{

	if( iBoun < 0 || iBoun >= getTotalNumberSubInnerBoundaries() ){
		OutputFiles::m_logFile << "Error : ID of boudary is out of range !! ID = " << iBoun << std::endl;
		exit(1);
	}

	return m_subInnerBoundaries[iBoun].getCoordXYOfPoints( iNode, &m_nodeList );

}

// Get number of inner boundaries within the specified outer boundary
int BoundaryCurveList::getNumInnerBoundaryIncluded( const int iBounOuter ) const{

	if( iBounOuter < 0 || iBounOuter >= getTotalNumberOuterBoundaries() ){
		OutputFiles::m_logFile << "Error : ID of boudary is out of range !! ID = " << iBounOuter << std::endl;
		exit(1);
	}

	return static_cast<int>( m_outer2InnerBoundaries[iBounOuter].size() );

}

// Get number of sub-inner boundaries within the specified inner boundary
int BoundaryCurveList::getNumSubInnerBoundaryIncluded( const int iBounInner ) const{

	if( iBounInner < 0 || iBounInner >= getTotalNumberInnerBoundaries() ){
		OutputFiles::m_logFile << "Error : ID of boudary is out of range !! ID = " << iBounInner << std::endl;
		exit(1);
	}

	return static_cast<int>( m_inner2SubInnerBoundaries[iBounInner].size() );

}

// Get ID of an inner boundary within the specified outer boundary
int BoundaryCurveList::getInnerBoundaryIDIncluded( const int iBounOuter, const int id ) const{

	if( id < 0 || id >= getNumInnerBoundaryIncluded(iBounOuter) ){
		OutputFiles::m_logFile << "Error : ID of boudary is out of range !! ID = " << id << std::endl;
		exit(1);
	}

	return m_outer2InnerBoundaries[iBounOuter][id];

}

// Get ID of a sub-inner boundary within the specified inner boundary
int BoundaryCurveList::getSubInnerBoundaryIDIncluded( const int iBounInner, const int id ) const{

	if( id < 0 || id >= getNumSubInnerBoundaryIncluded(iBounInner) ){
		OutputFiles::m_logFile << "Error : ID of boudary is out of range !! ID = " << id << std::endl;
		exit(1);
	}

	return m_inner2SubInnerBoundaries[iBounInner][id];

}

// Add node between the specified two points
void BoundaryCurveList::addNodesBetweenTwoPoins( const CommonParameters::XY& startCoord, const CommonParameters::XY& endCoord, 
											 const double maxLength, std::vector<int>& nodeIDs ){

	const double distance = hypot( endCoord.X - startCoord.X, endCoord.Y - startCoord.Y );

#ifdef _DEBUG_WRITE
	std::cout << "startCoord : " << startCoord.X << " " << startCoord.Y << std::endl;
	std::cout << "endCoord : " << endCoord.X << " " << endCoord.Y << std::endl;
	std::cout << "distance : " << distance << std::endl;
	std::cout << "maxLength : " << maxLength << std::endl;
#endif

	if( distance <= maxLength ){
		return;
	}

	const int numDiv = static_cast<int>(distance/maxLength) + 1;

	const double incX = ( endCoord.X - startCoord.X ) / static_cast<double>( numDiv );
	const double incY = ( endCoord.Y - startCoord.Y ) / static_cast<double>( numDiv );

	//NodeList* ptrNode2DList = NodeList::getInstance();

	for( int i = 1; i < numDiv; ++i ){
		CommonParameters::XY coord = { startCoord.X + incX * static_cast<double>(i), startCoord.Y + incY * static_cast<double>(i) };
		nodeIDs.push_back( m_nodeList.addNewNode( Node( coord, true ) ) );
#ifdef _DEBUG_WRITE
		std::cout << "nodeIDs.back() : " << nodeIDs.back() << std::endl;
#endif
	}

}

// Add node between the specified two points ( maximum edge is calculated automatically )
void BoundaryCurveList::addNodesBetweenTwoPoins( const CommonParameters::XY& startCoord, const CommonParameters::XY& endCoord, std::vector<int>& nodeIDs ){
	addNodesBetweenTwoPoins( startCoord, endCoord, calculateMaxEdgeLength( startCoord, endCoord ), nodeIDs );
}

// Calculate maximum edge length from specified two points
double BoundaryCurveList::calculateMaxEdgeLength( const CommonParameters::XY& startCoord, const CommonParameters::XY& endCoord ) const{
	const double leng1 = ( CoastLineList::getInstance() )->calcMaximumEdgeLengthForCoastLine( startCoord );
	const double leng2 = ( CoastLineList::getInstance() )->calcMaximumEdgeLengthForCoastLine( endCoord );
	return std::min(leng1, leng2);
}

// Check whether coastline is clockwise
bool BoundaryCurveList::isClockWise( std::vector<int>& nodeIDVec ) const{

	const int numNodes = static_cast<int>( nodeIDVec.size() );
	double areaSum(0.0);

	int nodeIDPre = nodeIDVec.back();
	for( std::vector<int>::iterator itr = nodeIDVec.begin(); itr != nodeIDVec.end(); ++itr ){
		const int nodeIDCur = *itr;
		const CommonParameters::XY coord0 = { 0.0, 0.0 };
		const CommonParameters::XY coord1 = m_nodeList.getCoordXYOfPoints(nodeIDPre);
		const CommonParameters::XY coord2 = m_nodeList.getCoordXYOfPoints(nodeIDCur);
		areaSum += Util::calcOuterProduct( coord0, coord1, coord2 );
		nodeIDPre = nodeIDCur;
	}
	if( areaSum > 0.0 ){
		return true;
	}
	else{
		return false;
	}
}

// Get flag specifing whether the specified nodes are the two end of a segment of a boundary curve
bool BoundaryCurveList::nodePairOfBoundaryCurve( const int nodeID0, const int nodeID1 ){

	for( std::vector<BoundaryCurveOuter>::iterator itr = m_outerBoundaries.begin(); itr != m_outerBoundaries.end(); ++itr ){
		if( itr->nodePairOfBoundaryCurve( nodeID0, nodeID1 ) ){
			return true;
		}
	}

	for( std::vector<BoundaryCurveInner>::iterator itr = m_innerBoundaries.begin(); itr != m_innerBoundaries.end(); ++itr ){
		if( itr->nodePairOfBoundaryCurve( nodeID0, nodeID1 ) ){
			return true;
		}
	}

	for( std::vector<BoundaryCurveSubInner>::iterator itr = m_subInnerBoundaries.begin(); itr != m_subInnerBoundaries.end(); ++itr ){
		if( itr->nodePairOfBoundaryCurve( nodeID0, nodeID1 ) ){
			return true;
		}
	}

	return false;

}

// Add new node
//int BoundaryCurveList::addNewNode( const CommonParameters::XY coord ){
int BoundaryCurveList::addNewNodeBetweenSegment( const int nodeID0, const int nodeID1 ){

	const CommonParameters::XY coord0 = m_nodeList.getCoordXYOfPoints( nodeID0 );
	const CommonParameters::XY coord1 = m_nodeList.getCoordXYOfPoints( nodeID1 );
	const CommonParameters::XY coordInsert = { 0.5*( coord0.X + coord1.X ), 0.5*( coord0.Y + coord1.Y ) };

	const int nodeInserted = m_nodeList.addNewNode( Node( coordInsert, false ) );

	bool found(false);
	for( std::vector<BoundaryCurveOuter>::iterator itr = m_outerBoundaries.begin(); itr != m_outerBoundaries.end(); ++itr ){
		if( itr->nodePairOfBoundaryCurve( nodeID0, nodeID1 ) ){
			itr->addNewNode( nodeID0, nodeID1, nodeInserted );
			found = true;
		}
	}

	for( std::vector<BoundaryCurveInner>::iterator itr = m_innerBoundaries.begin(); itr != m_innerBoundaries.end(); ++itr ){
		if( itr->nodePairOfBoundaryCurve( nodeID0, nodeID1 ) ){
			itr->addNewNode( nodeID0, nodeID1, nodeInserted );
			found = true;
		}
	}

	for( std::vector<BoundaryCurveSubInner>::iterator itr = m_subInnerBoundaries.begin(); itr != m_subInnerBoundaries.end(); ++itr ){
		if( itr->nodePairOfBoundaryCurve( nodeID0, nodeID1 ) ){
			itr->addNewNode( nodeID0, nodeID1, nodeInserted );
			found = true;
		}
	}

	if( !found ){
		m_nodeList.removeNode( nodeInserted );
		return -1;
	}

	return nodeInserted;

}

// Remove specified node
void BoundaryCurveList::removeNode( const int nodeID ){

	m_nodeList.removeNode( nodeID );

	for( std::vector<BoundaryCurveOuter>::iterator itr = m_outerBoundaries.begin(); itr != m_outerBoundaries.end(); ++itr ){
		itr->removeNode( nodeID );
	}

	for( std::vector<BoundaryCurveInner>::iterator itr = m_innerBoundaries.begin(); itr != m_innerBoundaries.end(); ++itr ){
		itr->removeNode( nodeID );
	}

	for( std::vector<BoundaryCurveSubInner>::iterator itr = m_subInnerBoundaries.begin(); itr != m_subInnerBoundaries.end(); ++itr ){
		itr->removeNode( nodeID );
	}

}

// Search the point of an outer boundary having specified coordinates and return whether the point belong to the boundary curve
bool BoundaryCurveList::hasPointOfGivenCoordInOuterBoundary( const CommonParameters::XY& coord, const int& iBounOuter, int& pointID ) const{

	if( iBounOuter < 0 || iBounOuter >= getTotalNumberOuterBoundaries() ){
		OutputFiles::m_logFile << "Error : ID of boudary is out of range !! ID = " << iBounOuter << std::endl;
		exit(1);
	}

	return m_outerBoundaries[iBounOuter].hasPointOfGivenCoord( coord, &m_nodeList, pointID );	

}

// Search point of outer boundaries having specified coordinates and return whether the point belong to the boundary curve
void BoundaryCurveList::searchPointOfOuterBoundaryByCoord( const CommonParameters::XY& coord, int& boundCurveID, int& pointID ) const{

	int icount(0);
	for( std::vector<BoundaryCurveOuter>::const_iterator itr = m_outerBoundaries.begin(); itr != m_outerBoundaries.end(); ++itr, ++icount ){
		if( hasPointOfGivenCoordInOuterBoundary( coord, icount, pointID ) ){
			boundCurveID = icount;
			return;
		}
	}

	OutputFiles::m_logFile << "Error : Could't find the point having specified coordinates in outer boundaries at the function " << __FUNCTION__ << " ." << std::endl;
	exit(1);

}

// Search point of outer boundaries having specified coordinates except the specified curve
void BoundaryCurveList::searchPointOfOuterBoundaryByCoordExceptACurve( const CommonParameters::XY& coord, const int exceptCurveID, int& boundCurveID, int& pointID ) const{

	int icount(0);
	for( std::vector<BoundaryCurveOuter>::const_iterator itr = m_outerBoundaries.begin(); itr != m_outerBoundaries.end(); ++itr, ++icount ){
		if( icount == exceptCurveID ){
			continue;
		}
		if( hasPointOfGivenCoordInOuterBoundary( coord, icount, pointID ) ){
			boundCurveID = icount;
			return;
		}
	}

	OutputFiles::m_logFile << "Error : Could't find the point having specified coordinates in outer boundaries at the function " << __FUNCTION__ << " ." << std::endl;
	exit(1);

}

// Exchange the types of boundary curves (Sea => Land, Land => Sea)
void BoundaryCurveList::exchangeTypesOfBoundaryCurves(){

	OutputFiles::m_logFile << "# Exchange the types of boundary curves (Sea => Land, Land => Sea)" << std::endl;

	for( std::vector<BoundaryCurveOuter>::iterator itr = m_outerBoundaries.begin(); itr != m_outerBoundaries.end(); ++itr ){
		itr->exchangeTypesOfBoundaryCurves();
	}

	for( std::vector<BoundaryCurveInner>::iterator itr = m_innerBoundaries.begin(); itr != m_innerBoundaries.end(); ++itr ){
		itr->exchangeTypesOfBoundaryCurves();
	}

	for( std::vector<BoundaryCurveSubInner>::iterator itr = m_subInnerBoundaries.begin(); itr != m_subInnerBoundaries.end(); ++itr ){
		itr->exchangeTypesOfBoundaryCurves();
	}

}

// Setup relation between nodes on boundary and lake data
void BoundaryCurveList::setupRelationBetweenNodesOnBoundaryAndLakeData(){

	LakeList* ptrLakeList = LakeList::getInstance();
	const int numLakeSurfPoints = ptrLakeList->getNumLakeSurfacePoints();

	for( std::vector<BoundaryCurveOuter>::const_iterator itr = m_outerBoundaries.begin(); itr != m_outerBoundaries.end(); ++itr ){
		if( itr->getGeologicalType() == CommonParameters::LAKE ){
			for( int iLake = 0; iLake < numLakeSurfPoints; ++iLake ){
				if( itr->included( ptrLakeList->getPointCoord(iLake), &m_nodeList ) ){
					const int numPoints = itr->getNumOfPoints();
					for( int iPoint = 0; iPoint < numPoints; ++iPoint ){
						m_nodeIDBoundCurveToLakeIndex.insert( std::make_pair( itr->getNodeID(iPoint), iLake ) );
					}
					break;
				}
			}
		}
	}

	for( std::vector<BoundaryCurveInner>::const_iterator itr = m_innerBoundaries.begin(); itr != m_innerBoundaries.end(); ++itr ){
		if( itr->getGeologicalType() == CommonParameters::LAKE ){
			for( int iLake = 0; iLake < numLakeSurfPoints; ++iLake ){
				if( itr->included( ptrLakeList->getPointCoord(iLake), &m_nodeList ) ){
					const int numPoints = itr->getNumOfPoints();
					for( int iPoint = 0; iPoint < numPoints; ++iPoint ){
						m_nodeIDBoundCurveToLakeIndex.insert( std::make_pair( itr->getNodeID(iPoint), iLake ) );
					}
					break;
				}
			}
		}
	}

	for( std::vector<BoundaryCurveSubInner>::const_iterator itr = m_subInnerBoundaries.begin(); itr != m_subInnerBoundaries.end(); ++itr ){
		if( itr->getGeologicalType() == CommonParameters::LAKE ){
			for( int iLake = 0; iLake < numLakeSurfPoints; ++iLake ){
				if( itr->included( ptrLakeList->getPointCoord(iLake), &m_nodeList ) ){
					const int numPoints = itr->getNumOfPoints();
					for( int iPoint = 0; iPoint < numPoints; ++iPoint ){
						m_nodeIDBoundCurveToLakeIndex.insert( std::make_pair( itr->getNodeID(iPoint), iLake ) );
					}
					break;
				}
			}
		}
	}

}

// Get lake height from coordinate
int BoundaryCurveList::getLakeIndexFromCoordinate( const CommonParameters::XY& coord ) const{

	for( std::vector<BoundaryCurveOuter>::const_iterator itr = m_outerBoundaries.begin(); itr != m_outerBoundaries.end(); ++itr ){
		if( itr->getGeologicalType() == CommonParameters::LAKE && itr->included( coord, &m_nodeList )){
			std::map<int,int>::const_iterator itrMap = m_nodeIDBoundCurveToLakeIndex.find(itr->getNodeID(1));
			if( itrMap != m_nodeIDBoundCurveToLakeIndex.end() ){
				return itrMap->second;
			}
			else{
				OutputFiles::m_logFile << "Error : Node " << itr->getNodeID(1) << " of boundary curves is not on lake lines." << std::endl;
				exit(1);
			}
		}
	}

	for( std::vector<BoundaryCurveInner>::const_iterator itr = m_innerBoundaries.begin(); itr != m_innerBoundaries.end(); ++itr ){
		if( itr->getGeologicalType() == CommonParameters::LAKE && itr->included( coord, &m_nodeList )){
			std::map<int,int>::const_iterator itrMap = m_nodeIDBoundCurveToLakeIndex.find(itr->getNodeID(1));
			if( itrMap != m_nodeIDBoundCurveToLakeIndex.end() ){
				return itrMap->second;
			}
			else{
				OutputFiles::m_logFile << "Error : Node " << itr->getNodeID(1) << " of boundary curves is not on lake lines." << std::endl;
				exit(1);
			}
		}
	}

	for( std::vector<BoundaryCurveSubInner>::const_iterator itr = m_subInnerBoundaries.begin(); itr != m_subInnerBoundaries.end(); ++itr ){
		if( itr->getGeologicalType() == CommonParameters::LAKE && itr->included( coord, &m_nodeList )){
			std::map<int,int>::const_iterator itrMap = m_nodeIDBoundCurveToLakeIndex.find(itr->getNodeID(1));
			if( itrMap != m_nodeIDBoundCurveToLakeIndex.end() ){
				return itrMap->second;
			}
			else{
				OutputFiles::m_logFile << "Error : Node " << itr->getNodeID(1) << " of boundary curves is not on lake lines." << std::endl;
				exit(1);
			}
		}
	}

	return -1;
}

// Get lake index from node ID on boundary curves
int BoundaryCurveList::getLakeIndexFromNodeID( const int nodeID ){

	std::map<int,int>::iterator itr = m_nodeIDBoundCurveToLakeIndex.find(nodeID);
	if( itr != m_nodeIDBoundCurveToLakeIndex.end() ){
		return itr->second;
	}
	else{
		OutputFiles::m_logFile << "Error : Node ( " << nodeID << " ) on boundary curves is not located on lake lines." << std::endl;
		exit(1);
		return -1;
	}

}

//
//// Add node between the specified two points
//void BoundaryCurveList::addNodesBetweenTwoPoins( const CommonParameters::XY& startCoord, const CommonParameters::XY& endCoord, 
//											 const double maxLength, std::vector<int>& nodeIDs ){
//
//	const double distance = hypot( endCoord.X - startCoord.X, endCoord.Y - startCoord.Y );
//
//#ifdef _DEBUG_WRITE
//	std::cout << "startCoord : " << startCoord.X << " " << startCoord.Y << std::endl;
//	std::cout << "endCoord : " << endCoord.X << " " << endCoord.Y << std::endl;
//	std::cout << "distance : " << distance << std::endl;
//	std::cout << "maxLength : " << maxLength << std::endl;
//#endif
//
//	if( distance <= maxLength ){
//		return;
//	}
//
//	const int numDiv = static_cast<int>(distance/maxLength) + 1;
//
//	const double incX = ( endCoord.X - startCoord.X ) / static_cast<double>( numDiv );
//	const double incY = ( endCoord.Y - startCoord.Y ) / static_cast<double>( numDiv );
//
//	//NodeList* ptrNode2DList = NodeList::getInstance();
//
//	for( int i = 1; i < numDiv; ++i ){
//		CommonParameters::XY coord = { startCoord.X + incX * static_cast<double>(i), startCoord.Y + incY * static_cast<double>(i) };
//		nodeIDs.push_back( m_nodeList.addNewNode( Node( coord, true ) ) );
//#ifdef _DEBUG_WRITE
//		std::cout << "nodeIDs.back() : " << nodeIDs.back() << std::endl;
//#endif
//	}
//
//}
