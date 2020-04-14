// This program should be capable of using SVD to solve the master equation describing the steady state transport of charge carriers in an array of molecules

#include "pch.h"
#include "site.h"
#include "simConsts.h"

int main()
{

    std::cout << "Creating sites...\n";
    std::vector<site> allSites = CreateSites();
    for (int i = 0; i < allSites.size(); i++)
        std::cout << allSites[i] << std::endl;

    std::cout << "\nCalculating transfer rates...\n";
    for (int i = 0; i < allSites.size(); i++)
        for (int f = 0; f < allSites.size(); f++)
        {
            double rate = allSites[i].Rate(&allSites[f]);
            if (rate) // if non-zero
                std::cout << "Rate " << i << "->" << f << " = " << rate << std::endl;
        }


}

std::vector<site> CreateSites()
{
    std::vector<site> sites;

    // For now, hardcode some information about the array
    double E = 0.0;
    double sizeX = 20.0, sizeY = 20.0, sizeZ = 20.0;
    double periodX = 10.0, periodY = 10.0, periodZ = 10.0;

    // axis index
    int i_x = 0, i_y = 0, i_z = 0;

    // Create sites. Currently just using a cubic array. 
    while (i_z * periodZ <= sizeZ)
    {
        while (i_y * periodY <= sizeY)
        {
            while (i_x * periodX <= sizeX)
            {
                sites.push_back(site(i_x * periodX, i_y * periodY, i_z * periodZ, E));
                i_x++;
            }
            i_x = 0;
            i_y++;
        }
        i_y = 0;
        i_z++;
    }

    // Assign interacting neighbours. Currently we are just using the up-to 6 nearest neighbours.
    double J = 0.1;

    i_x = 0, i_y = 0, i_z = 0;

    // highest index in each direction
    int n_x = (int)std::floor(sizeX / periodX);
    int n_y = (int)std::floor(sizeY / periodY);
    int n_z = (int)std::floor(sizeZ / periodZ);

    // Current site index
    int i_p = 0;

    while (i_z * periodZ <= sizeZ)
    {
        while (i_y * periodY <= sizeY)
        {
            while (i_x * periodX <= sizeX)
            {
                // If the current site is not the closest to the upper bound in each axis,
                // then assign it the next site in each direction as a neighbor.
                // This also automatically assigns the current site as a neighbor to each of the 'nexts',
                // therefore we only need to specify 'forward/positive' neighbors here.
                if (i_x != n_x) { sites[i_p].addNeighbour(&sites[i_p + 1], J); }
                if (i_y != n_y) { sites[i_p].addNeighbour(&sites[i_p + n_x + 1], J); }
                if (i_z != n_z) { sites[i_p].addNeighbour(&sites[i_p + (n_x + 1) * (n_y + 1)], J); }

                i_p++;
                i_x++;
            }
            i_x = 0;
            i_y++;
        }
        i_y = 0;
        i_z++;
    }

    return sites;

}

bool CreateA(gsl_matrix* A) 
{
    return false;
}
