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
#ifndef DBLDEF_OUTPUT_FILES
#define DBLDEF_OUTPUT_FILES

#include <fstream>

// Class of output files
class OutputFiles{

public:
	// Return the the instance of the class
	static OutputFiles* getInstance();

	// Log file
	static std::ofstream m_logFile;

private:
	// Constructer
	OutputFiles();

	// Destructer
	~OutputFiles();

	// Copy constructer
	OutputFiles(const OutputFiles& rhs);

	// Assignment operator
	OutputFiles& operator=(const OutputFiles& rhs);
		
};

#endif
