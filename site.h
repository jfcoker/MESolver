#pragma once

#include "pch.h"

class site 
{
private:

	struct vec { double X, Y, Z;  };

	struct neighbour {
		site* _pSite; 
		double _J;
		//neighbour(site* pSite, double J) : _pSite(pSite), _J(J) {}
	};

	vec pos;
	double energy;
	std::vector<neighbour*> neighbours;

public:
	site(double x, double y, double z, double energy);

	void addNeighbour(site* pSite, double J);

	~site();

	friend std::ostream& operator<<(std::ostream& os, const site& st);

};



