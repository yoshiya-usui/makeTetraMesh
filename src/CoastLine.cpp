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
#include "CoastLine.h"
#include "OutputFiles.h"
#include "AnalysisDomain.h"
#include "Util.h"
#include "Control.h"
#include "math.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

// Default constructer
CoastLine::CoastLine():
	m_isClosed(false)
{
}

// Constructer
CoastLine::CoastLine( const bool isClosed ):
	m_isClosed(isClosed)
{
}

// Destructer
CoastLine::~CoastLine(){

}

// Read data of coast lines from input file
void CoastLine::readCoastLineData( std::ifstream& ifs ){

	if( !m_nodes.empty() ){
		m_nodes.clear();
	}

	if( !m_fix.empty() ){
		m_fix.clear();
	}

	const bool invertYCoord = ( Control::getInstance() )->getInvertSignYcoord();

	int icount(0);
	while( !ifs.eof() ) {

		CommonParameters::XY coord = { 0.0 , 0.0 };
		int iflag(0);
		int ifix(0);
		
		ifs >> coord.X >> coord.Y >> iflag >> ifix;

		if( invertYCoord ){
			coord.Y *= -1.0;
		}

		bool bfix(false);
		if( ifix == 1 ){
			bfix = true;
		}

		m_nodes.push_back( Node( coord, bfix ) );
		
		if( iflag == 1 ){
			m_isClosed = true;
			return;
		}else if( iflag == -1 ){
			m_isClosed = false;
			return;
		}

		++icount;

#ifdef _DEBUG_WRITE
		std::cout << "icount m_nodes m_fix m_isClosed : " << icount << " " << (m_nodes.back()).getCoordX() << " " << (m_nodes.back()).getCoordY() << " " << (m_nodes.back()).isFixed() << " " << m_isClosed << std::endl;
#endif

	}

	OutputFiles::m_logFile << "Error : Reach end of coast_line.dat while reading data !!" << std::endl;
	exit(1);
	
}

// Get total number of points
int CoastLine::getNumOfPoints() const{
	return static_cast<int>( m_nodes.size() );
}

// Get coordinate X of points
double CoastLine::getCoordXOfPoints( const int id ) const{

	if( id < 0 || id >= getNumOfPoints() ){
		//OutputFiles::m_logFile << "Error : ID of point is out of range !! ID = " << id << std::endl;
		OutputFiles::m_logFile << "Error : ID of point is out of range at the function " << __FUNCTION__  << " !! ID = " << id << std::endl;
		exit(1);
	}

	return m_nodes[id].getCoordX();

}

// Get coordinate Y of points
double CoastLine::getCoordYOfPoints( const int id ) const{

	if( id < 0 || id >= getNumOfPoints() ){
		//OutputFiles::m_logFile << "Error : ID of point is out of range !! ID = " << id << std::endl;
		OutputFiles::m_logFile << "Error : ID of point is out of range at the function " << __FUNCTION__  << " !! ID = " << id << std::endl;
		exit(1);
	}

	return m_nodes[id].getCoordY();

}

// Get coordinate of point
CommonParameters::XY CoastLine::getCoordXYOfPoints( const int id ) const{

	if( id < 0 || id >= getNumOfPoints() ){
		OutputFiles::m_logFile << "Error : ID of point is out of range at the function " << __FUNCTION__  << " !! ID = " << id << std::endl;
		exit(1);
	}

	return m_nodes[id].getCoordXY();

}

// Get first coordinate of points
CommonParameters::XY CoastLine::getFirstCoordOfPoints() const{

	return (m_nodes.front()).getCoordXY();

}

// Get last coordinate of points
CommonParameters::XY CoastLine::getLastCoordOfPoints() const{

	return (m_nodes.back()).getCoordXY();

}

// Get node data
Node CoastLine::getNode( const int id ) const{

	if( id < 0 || id >= getNumOfPoints() ){
		OutputFiles::m_logFile << "Error : ID of point is out of range at the function " << __FUNCTION__  << " !! ID = " << id << std::endl;
		exit(1);
	}

	return m_nodes[id];

}

// Get flags specifing whether the point must be included or not
bool CoastLine::isFix( const int id ) const{

	if( id < 0 || id >= getNumOfPoints() ){
		OutputFiles::m_logFile << "Error : ID of point is out of range at the function " << __FUNCTION__  << " !! ID = " << id << std::endl;
		exit(1);
	}

	return m_nodes[id].isFixed();

}

// Get flag specifing the coast line is closed or not
bool CoastLine::isClosed() const{
	return m_isClosed;
}

// Remove the nodes consisting to too sharp angle and the nodes whose distances are too short
void CoastLine::coastLineSmoothing( std::vector<Node>& stack ){

	OutputFiles::m_logFile << "# Coast line smoothing " << std::endl;

	// Remove the nodes consisting to too sharp angle and the nodes whose distances are too short
	const double numIterMax = ( Control::getInstance() )->getMaxIterNumCoastLineSmoothing();
	const double thresholdAngle = CommonParameters::DEG2RAD * ( Control::getInstance() )->getThresholdAngleCoastLineSmoothing();
	const double thresholdDistance = ( Control::getInstance() )->getThresholdDistanceCoastLineSmoothing();

	for( int iter = 0; iter < numIterMax; ++iter ){

		bool satisifyThresholds(true);

		// Delete the nodes consisting to too sharp angle
		for( std::vector<Node>::iterator itr1 = stack.begin(); itr1 != stack.end(); ++itr1 ){

			std::vector<Node>::iterator itr2 = itr1;
			if( ++itr2 == stack.end() ){
				break;
			}
			std::vector<Node>::iterator itr3 = itr2;
			if( ++itr3 == stack.end() ){
				break;
			}

			const double angle = fabs( Util::calcAngle( itr1->getCoordXY(), itr2->getCoordXY(), itr3->getCoordXY() ) );

			if( angle < thresholdAngle ){
#ifdef _DEBUG_WRITE
				std::cout << "Too sharp angle !!" << std::endl;
				std::cout << "Angle is " << angle << ". Threshold value is " << thresholdAngle << std::endl;
#endif

				itr1 = stack.erase( itr2 );
				satisifyThresholds = false;
			}

		}

		// Delete the nodes whose distances are too short
		std::vector<Node> stack2;
		for( std::vector<Node>::iterator itr1 = stack.begin(); itr1 != stack.end(); ++itr1 ){

			std::vector<Node>::iterator itr2 = itr1;
			++itr2;
			for( ; itr2 != stack.end(); ++itr2 ){
				if( itr2->calcHorizontalDistance( *itr1 ) < thresholdDistance ){// Two node is too close
#ifdef _DEBUG_WRITE
					std::cout << "Too short !!" << std::endl;
					std::cout << "Distance betwenn two node is " << itr2->calcHorizontalDistance( *itr1 ) << ". Threshold value is " << thresholdDistance << std::endl;
#endif
					itr1 = itr2;
					satisifyThresholds = false;
				}
			}

			stack2.push_back( *itr1 );

		}
		stack.swap( stack2 );

		if( satisifyThresholds ){
			return;
		}

	}

}

// Roughen coast line
void CoastLine::roughenCoastLine(){

	OutputFiles::m_logFile << "# Roughen coast line" << std::endl;

	if( !m_isClosed ){// Skip if coast line is not closed
		OutputFiles::m_logFile << "Warning : This coast line is not closed." << std::endl;
	}
	
	// Not delete for future use >>>>>
	//const AnalysisDomain* const ptrAnalysisDomain = AnalysisDomain::getInstance();

	//const double xMin = ptrAnalysisDomain->getMinCoordX(); 
	//const double xMax = ptrAnalysisDomain->getMaxCoordX(); 
	//const double yMin = ptrAnalysisDomain->getMinCoordY(); 
	//const double yMax = ptrAnalysisDomain->getMaxCoordY(); 

	//// Add intersection point with boundary analysis domain to coast line data
	//CommonParameters::XY coordPre = (m_nodes.back()).getCoordXY();
	//bool locateInsidePreviously = false;
	//if( xMin <= coordPre.X && coordPre.X <= xMax && yMin <= coordPre.Y && coordPre.Y <= yMax ){
	//	locateInsidePreviously = true;
	//}

	//for( std::vector< Node >::iterator itr = m_nodes.begin(); itr != m_nodes.end(); ++itr ){
	//	
	//	CommonParameters::XY coord = itr->getCoordXY();

	//	if( xMin <= coord.X && coord.X <= xMax && yMin <= coord.Y && coord.Y <= yMax ){
	//		// This point locate inside of the analysis domain

	//		if( !locateInsidePreviously ){
	//			m_nodes.insert( itr, Node( Util::getCoordOfIntersectionPoint( coordPre, coord ) , true ) );
	//		}

	//		locateInsidePreviously = true;

	//	}else{
	//		// This point locate out of the analysis domain

	//		if( locateInsidePreviously ){
	//			m_nodes.insert( itr, Node( Util::getCoordOfIntersectionPoint( coordPre, coord ) , true ) );
	//		}

	//		locateInsidePreviously = false;
	//	}

	//	coordPre = coord;

	//}
	// Not delete for future use <<<<<

	const double thresholdEdgeLengthRatio = ( Control::getInstance() )->getThresholdEdgeLengthRatio();

	std::vector<Node>::iterator itr = m_nodes.begin();
	
	std::vector<Node> stack;
	stack.push_back( *itr );
	//std::pair<Node*,Node*> nodePair( &(*itr), &(*itr));
	Node* ptrSegmentStart = &(*itr);
	Node* ptrNodePre = &(*itr);
	double lengthCur( 0.0 );
	bool addNode2DPre = true;
	++itr;

	for( ; itr != m_nodes.end(); ++itr ){

		const double segmentLength = itr->calcHorizontalDistance( *ptrNodePre );
		lengthCur += segmentLength;
		const double maxDistance = std::min( ptrSegmentStart->calcMaximumEdgeLengthForCoastLine(), itr->calcMaximumEdgeLengthForCoastLine() );

#ifdef _DEBUG_WRITE
		std::cout << "segmentLength maxDistance : " << segmentLength << " " << maxDistance << std::endl;
#endif

		if( segmentLength > maxDistance ){

			// Replace last node of stack with previous node if length ratio of the two segments having last node of stack is greater than ths specified value
			if( !addNode2DPre && !(stack.back()).isFixed() && static_cast<int>( stack.size() ) > 2 ){// If previous point is free
				std::vector<Node>::iterator itrBfr = stack.end();
				--itrBfr;
				if( itrBfr->calcHorizontalDistance( stack.back() ) / ptrNodePre->calcHorizontalDistance( stack.back() ) > thresholdEdgeLengthRatio ){
					stack.back() = *ptrNodePre;
					addNode2DPre = true;
				}
			}

			// Divide segment if segment is too long
			const int numDiv = static_cast<int>( segmentLength / maxDistance ) + 1;
			const double incrementX = ( itr->getCoordX() - ptrNodePre->getCoordX() ) / numDiv;
			const double incrementY = ( itr->getCoordY() - ptrNodePre->getCoordY() ) / numDiv;
			for( int i = 1; i < numDiv; ++i ){
				CommonParameters::XY coord = { ptrNodePre->getCoordX() + incrementX * static_cast<double>(i), ptrNodePre->getCoordY() + incrementY * static_cast<double>(i) };
				stack.push_back( Node( coord , false ) );
			}

			stack.push_back( *itr );
			ptrSegmentStart = &(*itr);
			ptrNodePre = &(*itr);
			lengthCur = 0.0;
			addNode2DPre = true;
			continue;
		}

		if( itr->isFixed() ){

			// Delete previous node if length ratio of the two segments having this node is greater than ths specified value
			if( !(stack.back()).isFixed() && static_cast<int>( stack.size() ) > 2 ){// If previous point is free
				std::vector<Node>::iterator itrBfr = stack.end();
				--itrBfr;
				if( itrBfr->calcHorizontalDistance( stack.back() ) / itr->calcHorizontalDistance( stack.back() ) > thresholdEdgeLengthRatio ){
					stack.pop_back();
				}
			}

			stack.push_back( *itr );
			ptrSegmentStart = &(*itr);
			ptrNodePre = &(*itr);
			lengthCur = 0.0;
			addNode2DPre = true;
			continue;
		}
	
#ifdef _DEBUG_WRITE
		std::cout << "lengthCur = "<< lengthCur << std::endl;
#endif

		if( lengthCur > maxDistance ){
			if( !addNode2DPre ){
				stack.push_back( *ptrNodePre );
				ptrSegmentStart = ptrNodePre;
				ptrNodePre = &(*itr);
				lengthCur = segmentLength;
				addNode2DPre = false;
				continue;
			}else{
				OutputFiles::m_logFile << "Error : addNode2DPre must be false at this point !!" << std::endl;
				exit(1);
			}
		}

		ptrNodePre = &(*itr);
		addNode2DPre = false;

	}

	if( !stack.empty() ){

		const double segmentLength = stack.front().calcHorizontalDistance( stack.back() );
		const double maxDistance = std::min( (stack.front()).calcMaximumEdgeLengthForCoastLine(), (stack.back()).calcMaximumEdgeLengthForCoastLine() );

#ifdef _DEBUG_WRITE
		std::cout << "segmentLength maxDistance : " << segmentLength << " " << maxDistance << std::endl;
#endif

		// Divide segment if segment is too long
		if( segmentLength > maxDistance ){
			const int numDiv = static_cast<int>( segmentLength / maxDistance ) + 1;
			const double X0 = (stack.back()).getCoordX();
			const double Y0 = (stack.back()).getCoordY();
			const double incrementX = ( (stack.front()).getCoordX() - X0 ) / numDiv;
			const double incrementY = ( (stack.front()).getCoordY() - Y0 ) / numDiv;

			for( int i = 1; i < numDiv; ++i ){
				CommonParameters::XY coord = { X0 + incrementX * static_cast<double>(i), Y0 + incrementY * static_cast<double>(i) };
#ifdef _DEBUG_WRITE
				std::cout << "coord " << coord.X << " " << coord.Y << std::endl;
#endif
				stack.push_back( Node( coord , false ) );
			}
		}

		// Delete first node if length ratio of the two segments having this node is greater than ths specified value
		if( !(stack.front()).isFixed() && static_cast<int>( stack.size() ) > 2 ){// If first point is free
			if( stack[1].calcHorizontalDistance( stack.front() ) / (stack.back()).calcHorizontalDistance( stack.front() ) > thresholdEdgeLengthRatio ){
				stack.erase( stack.begin() );
			}
		}

	}

	// Remove the nodes consisting to too sharp angle and the nodes whose distances are too short
	coastLineSmoothing( stack );

	m_nodes.swap( stack );

	if( static_cast<int>(m_nodes.size() )< 3 ){
		m_isClosed = false;
	}

}
