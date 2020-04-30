#pragma once
#include "pch.h"

class site 
{
private:

	// A structure that will be used to hold the position of the site.
	struct vec { double X, Y, Z; };

	// A structure that will be used to specify which other sites this site interacts with, and the strength of this interaction.
	struct neighbour;

	// The list of neighbors this site interacts with.
	std::vector<neighbour*> neighbours;

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

	// Check if this site has the passed pointed-to site as a neighbour.
	// If it does, return a pointer to the neighbour object.
	// Otherwise return NULL
	neighbour* hasNeighbour(site* pSite);

	// Calculate the energetic driving force for the transfer of a charge from this site (as initial) to the site passed as pointer (as final).
	double deltaE(site* pSite, double fieldZ);

	// Calculate the transfer rate between this site (as initial), and the site passed as pointer (as final).
	// If the passed site is not in the list of interacting neighbours, the rate will be zero.
	// Only performs the full calculation the first time this function is called (for this this specific rate).
	double Rate(site* pSite, double fieldZ, double kBT, double reorg);

	// Alternative forms of preconditioning factor
	// (Rather than enum could treat site as an interface, 
	//  with function PrecondFactor as a pure virtual function,
	//  and each form as a different implementation)
	enum class PrecondForm { boltzmann, boltzmannSquared, rateSum};

	// Calculate the preconditioning factor.
	// This is used to transform the rate matrix into a form more suitable for solving numerically.
	double PrecondFactor(double fieldZ, double kBT, double reorg, double E0, PrecondForm form, bool apply);

	// Destructor
	~site();

	friend std::ostream& operator<<(std::ostream& os, const site& st);

private:

	struct neighbour 
	{
		friend double site::Rate(site* pSite, double fieldZ, double kBT, double reorg);
		site* _pSite = NULL;
		double _J = 0.0;
	private: // Make sure only the function Rate() has access to the member below.
		double _rate = -1.0; // -1 for not yet set, the rate will be updated the first time it is needed.
	};

};

std::vector<site> CreateSites(char* XYZfile, char* EDGEfile);



