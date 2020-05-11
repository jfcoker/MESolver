#pragma once
#include "pch.h"
#include "site.h"

class transporter
{
private:

	const double _kBT;
	const double _fieldZ;
	const double _reorg;
	const double _transE;
	const bool _periodic;
	const double _sizeZ;
	const double _rsizeZ;

public:

	// Construct a transporter object
	transporter(double kBT, double fieldZ, double reorg, double transE, bool periodic, double sizeZ);

	// Calculate the energetic driving force for the transfer of a charge from this site (as initial) to the site passed as pointer (as final).
	// If periodic is true then use the minimum image convention to find the shortest path to from this site to the target site.
	double deltaE(site* orig, site* dest);

	// Calculate the transfer rate between two sites.
	// If the passed sites are not interacting, the rate will be zero.
	// Only performs the full calculation the first time this function is called (for this this specific rate).
	double Rate(site* orig, site* dest);

	// Calculate the sum of all transfer rates into the specified site.
	double RateSum(site* dest);

	// Alternative forms of preconditioning factor
	// (Rather than enum could treat site as an interface, 
	//  with function PrecondFactor as a pure virtual function,
	//  and each form as a different implementation)
	enum class PrecondForm { off, boltzmann, boltzmannSquared, rateSum };

	// Calculate the preconditioning factor.
	// This is used to transform the rate matrix into a form more suitable for solving numerically.
	double PrecondFactor(site* pSite, PrecondForm form);

	gsl_matrix* CreateRateMatrix(std::vector<site>& sites, PrecondForm form, bool scale, bool verbose);

	// Use the occupation probability of sites and the rate equation, to find the average velocity of charges in Ang/s
	double velocity_z(std::vector<site>& sites, gsl_matrix* A);

};

// A structure used to specify
// which other sites the parent site interacts with,
// and the strength of interactions.
// Full definition in transporter.h
struct site::neighbour
{
	site* _pSite = NULL;
	double _J = 0.0;
	private: // Make sure only the function Rate() has access to the member below.
		friend double transporter::Rate(site* orig, site* dest);
		double _rate = -1.0; // -1 for not yet set, the rate will be updated the first time it is needed.
};





