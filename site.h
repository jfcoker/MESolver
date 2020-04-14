#pragma once

#include "pch.h"

class site 
{
private:

	// The
	struct vec { double X, Y, Z;  };

	struct neighbour {
		site* _pSite; 
		double _J;
	};

	vec pos;
	double energy;
	std::vector<neighbour*> neighbours;

public:
	
	// Construct a site object
	site(double x, double y, double z, double energy);

	// Give this site a pointer to another site which it interacts with, alongside the associated transfer integral
	void addNeighbour(site* pSite, double J);

	// Calcualte the transfer rate between this site, and the site passed as pointer.
	// If the passed site is not in the list of interacting neighbours, the rate will be zero.
	double Rate(site* pSite);

	~site();

	friend std::ostream& operator<<(std::ostream& os, const site& st);

};

std::vector<site> CreateSites();



