#pragma once
#include "pch.h"

//class transporter;

class site 
{
private:

	// Transporter can access private members of site.
	friend class transporter;

	// A structure that will be used to hold the position of the site.
	struct vec { double X, Y, Z; };

	// Forward declare a structure that will be used to specify
	// which other sites the parent site interacts with,
	// and the strength of interactions.
	// Full definition in transporter.h
	struct neighbour;

	// The list of neighbors this site interacts with.
	std::vector<neighbour*> neighbours;

	// Check if this site has the passed pointed-to site as a neighbour.
	// If it does, return a pointer to the neighbour object.
	// Otherwise return NULL
	neighbour* hasNeighbour(site* pSite);

public:

	// The position of the site.
	vec pos;

	// The occupation probability of this site. -1.0 for not set.
	double occProb = -1.0;

	// The site energy.
	double energy;
	
	// Construct a site object
	site(double x, double y, double z, double energy);

	// Give this site a pointer to another site which it interacts with, alongside the associated transfer integral
	void addNeighbour(site* pSite, double J);

	// Destructor
	~site();

	friend std::ostream& operator<<(std::ostream& os, const site& st);

};

std::vector<site> CreateSites(char* XYZfile, char* EDGEfile);
