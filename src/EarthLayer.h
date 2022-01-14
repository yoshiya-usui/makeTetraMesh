#ifndef DBLDEF_EARTH_LAYER
#define DBLDEF_EARTH_LAYER

#include <vector>

// Class of the earth layer
class EarthLayer{

public:
	// Return the the instance of the class
    static EarthLayer* getInstance();

	// Read data of anomalies from input file
	void readEarthLayerData();
	
private:
	// Default constructer
	EarthLayer();

	// Destructer
	~EarthLayer();

	// Copy constructer
	EarthLayer(const EarthLayer& rhs);

	// Assignment operator
	EarthLayer& operator=(const EarthLayer& rhs);

	// Total number of layers
	int m_numAnomalyBrick;

	// Thicknesses of layers
	double* m_thickness;

	// Resistivities of layers
	double* m_resistivities;

};

#endif
