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
