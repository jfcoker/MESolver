#pragma once
#include "pch.h"

class site 
{
private:

	// A structure that will be used to hold the position of the site.
	struct vec { double X, Y, Z;  };

	// A structure that will be used to specify which other sites this site interacts with, and the strength of this interaction.
	struct neighbour {
		site* _pSite; 
		double _J;
	};
	
	double energy;
	std::vector<neighbour*> neighbours;

public:
	
	vec pos;
	double occProb = 0;
	
	// Construct a site object
	site(double x, double y, double z, double energy);

	// Give this site a pointer to another site which it interacts with, alongside the associated transfer integral
	void addNeighbour(site* pSite, double J);

	// Check if this site has the passed pointed-to site as a neighbour.
	// If it does, return a pointer to the neighbour object.
	// Otherwise return NULL
	neighbour* hasNeighbour(site* pSite);

	// Calcualte the transfer rate between this site (as initial), and the site passed as pointer (as final).
	// If the passed site is not in the list of interacting neighbours, the rate will be zero.
	double Rate(site* pSite, double fieldZ, double kBT, double reorg);

	~site();

	friend std::ostream& operator<<(std::ostream& os, const site& st);

};

// Arguments specify number of sites in each dimension.
// Currently seperation is hardcoded to 10 Ang.
std::vector<site> CreateSites(size_t X, size_t Y, size_t Z);

std::vector<site> CreateSites(char* XYZfile, char* EDGEfile);



