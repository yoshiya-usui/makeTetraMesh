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
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include "OutputFiles.h"
#include "CommonParameters.h"

std::ofstream OutputFiles::m_logFile;

// Return the the instance of the class
OutputFiles* OutputFiles::getInstance(){
   	static OutputFiles instance;// The only instance
  	return &instance;
}

// Constructer
OutputFiles::OutputFiles(){

	std::ostringstream logFileName;
	logFileName << CommonParameters::programName << ".log";
	m_logFile.open( logFileName.str().c_str(), std::ios::out );

	if( m_logFile.fail() )
	{
		std::cerr << "File open error !! : " << logFileName.str() << std::endl;
		exit(1);
	}

	OutputFiles::m_logFile << "# " << CommonParameters::programName << " " << CommonParameters::version << std::endl;

}

// Destructer
OutputFiles::~OutputFiles(){

	if( m_logFile.is_open() ){
		m_logFile.close();
	}

}
