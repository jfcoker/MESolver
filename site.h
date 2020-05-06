#pragma once
#include "pch.h"

class site;

// A structure that will be used to specify which other sites this site interacts with, and the strength of this interaction.
struct neighbour
{
	//friend double transporter::Rate(site* orig, site* dest);
	site* _pSite = NULL;
	double _J = 0.0;
	//private: // Make sure only the function Rate() has access to the member below.
	double _rate = -1.0; // -1 for not yet set, the rate will be updated the first time it is needed.
};

class site 
{
private:

	// A structure that will be used to hold the position of the site.
	struct vec { double X, Y, Z; };

public:

	// The position of the site.
	vec pos;

	// The occupation probability of this site. -1.0 for not set.
	double occProb = -1.0;

	// The site energy.
	double energy;

	// The list of neighbors this site interacts with.
	std::vector<neighbour*> neighbours;
	
	// Construct a site object
	site(double x, double y, double z, double energy);

	// Give this site a pointer to another site which it interacts with, alongside the associated transfer integral
	void addNeighbour(site* pSite, double J);

	// Check if this site has the passed pointed-to site as a neighbour.
	// If it does, return a pointer to the neighbour object.
	// Otherwise return NULL
	neighbour* hasNeighbour(site* pSite);

	// Destructor
	~site();

	friend std::ostream& operator<<(std::ostream& os, const site& st);

};

std::vector<site> CreateSites(char* XYZfile, char* EDGEfile);
