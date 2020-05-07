#include "pch.h"
#include "site.h"
#include "transporter.h"
#include "consts.h"
#include "utility.h"


// Construct a transporter object
transporter::transporter(double kBT, double fieldZ, double reorg, double transE, bool periodic, double sizeZ) : 
	_kBT(kBT),
	_fieldZ(fieldZ),
	_reorg(reorg),
    _transE(transE),
	_periodic(periodic),
	_sizeZ(sizeZ),
	_rsizeZ(1.0 / sizeZ)
{}

// Calculate the energetic driving force for the transfer of a charge from this site (as initial) to the site passed as pointer (as final).
// If periodic is true then use the minimum image convention to find the shortest path to from this site to the target site.
double transporter::deltaE(site* orig, site* dest)
{
	double deltaZ = (dest->pos.Z - orig->pos.Z);

	// If periodic boundaries in z, then apply the minimum image convention.
	if (_periodic) deltaZ -= _sizeZ * floor(deltaZ * _rsizeZ + 0.5);

	return (dest->energy - orig->energy) + deltaZ * _fieldZ;

}

// Calculate the transfer rate between two sites.
// If the passed sites are not interacting, the rate will be zero.
// Only performs the full calculation the first time this function is called (for this this specific rate).
double transporter::Rate(site* orig, site* dest)
{
    site::neighbour* pN = orig->hasNeighbour(dest);

    if (!pN)
        // Sites aren't interacting, transfer rate will be zero.
        return 0.0;
    else if (pN->_rate < 0.0) // If rate is -1, then this is the initial call and rate needs to be calculated. Using '<' to avoid equality comparison on floating point values.
    {
        double J2 = std::pow(pN->_J, 2); // |J_if|^2
        pN->_rate = ((2 * pi) / hbar) * J2 * std::pow(4 * pi * _reorg * _kBT, -0.5) * std::exp(-1 * std::pow(deltaE(orig,dest) + _reorg, 2) / (4 * _reorg * _kBT));
        return pN->_rate;
    }
    else // Rate has already been calculated. Just grab it.
        return pN->_rate;

}

// Calculate the sum of all transfer rates into the specified site.
double transporter::RateSum(site* dest)
{
    double sum = 0;
    std::vector<site::neighbour*>::iterator it = dest->neighbours.begin();
    for (int i = 0; it != dest->neighbours.end(); i++, it++)
        sum += Rate((*it)->_pSite, dest);

    return sum;
}

// Calculate the preconditioning factor.
// This is used to transform the rate matrix into a form more suitable for solving numerically.
double transporter::PrecondFactor(site* pSite, PrecondForm form)
{
    // Could use interface / implementation instead of enum / switch
    switch (form)
    {

    case PrecondForm::off:
        return 1.0;

    case PrecondForm::boltzmann:
        return std::exp((_transE - (pSite->energy + pSite->pos.Z * _fieldZ)) / _kBT);

    case PrecondForm::boltzmannSquared:
        return std::pow(std::exp((_transE - (pSite->energy + pSite->pos.Z * _fieldZ)) / _kBT), 2.0);

    case PrecondForm::rateSum:
        return 1.0 / RateSum(pSite);

    default:
        throw std::logic_error("Form of preconditioning factor not implemented.");
        return -1.0;

    }
}

gsl_matrix* transporter::CreateRateMatrix(std::vector<site>& sites, PrecondForm form, bool scale, bool verbose)
{
    size_t M = sites.size();
    gsl_matrix* A = gsl_matrix_alloc(M, M);
    gsl_matrix_set_zero(A);

    int highestO = -999;
    int lowestO = 999;
    int orderOfMag;
    for (int i = 0; i < M; i++)
        for (int f = 0; f < M; f++)
        {
            double el = 0;
            if (i == f)
            {
                double sum = 0.0;
                for (int k = 0; k < M; k++)
                    if (i != k)
                        sum += Rate(&sites[i], &sites[k]) * PrecondFactor(&sites[i], form);

                el = -sum;
            }
            else
            {
                el = Rate(&sites[f], &sites[i]) * PrecondFactor(&sites[f], form);
            }
            gsl_matrix_set(A, i, f, el);
            if (el) // Check el is non-zero otherwise lowestO will equal -inf
            {
                orderOfMag = (int)floor(log10(abs(el)));
                if (orderOfMag > highestO) highestO = orderOfMag;
                if (orderOfMag < lowestO) lowestO = orderOfMag;
            }
        }

    if (verbose)
    {
        printMatrix(A);
        std::cout << "\nValues:\nMax = " << gsl_matrix_max(A) << "\nMin = " << gsl_matrix_min(A) << "\nRange = " << gsl_matrix_max(A) - gsl_matrix_min(A) << "\n";
        std::cout << "\nOrder of magnitude:\nHighest = " << highestO << "\nLowest = " << lowestO << "\nDiff = " << highestO - lowestO << "\n";
    }

    if (scale)
    {
        std::cout << "\nTo reduce precision errors, rescale A by 1e-" << highestO << "\n";
        gsl_matrix_scale(A, pow(10, -highestO));

        if (verbose)
        {
            std::cout << "Reduced A = \n";
            printMatrix(A);
        }
    }

    return A;

}

double transporter::velocity_z(std::vector<site>& sites, gsl_matrix* A)
{
    double sum = 0.0;
    for (int i = 0; i < sites.size(); i++)
        for (int j = 0; j < sites.size(); j++)
            if (i != j)
                sum += (sites[i].pos.Z - sites[j].pos.Z) * gsl_matrix_get(A, j, i) * sites[i].occProb;

    return sum;
}