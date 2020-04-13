// This program should be capable of using SVD to solve the master equation describing the steady state transport of charge carriers in an array of molecules

#include "pch.h"
#include "site.h"

int main()
{
    const double reorg = 0.01;
    const double kBT = -999;

    const double pi = 3.14;
    const double hbar = -999;

    std::cout << "Creating sites...\n";

    std::vector<site> allSites = CreateSites();

    for (int i = 0; i < allSites.size(); i++)
        std::cout << allSites[i] << std::endl;

}

std::vector<site> CreateSites()
{
    std::vector<site> sites;

    // Hardcode some information about the array
    double E = 0.0;
    double sizeX = 3.0, sizeY = 3.0, sizeZ = 3.0;
    double periodX = 1.0, periodY = 1.0, periodZ = 1.0;

    // index
    int i_x = 0, i_y = 0, i_z = 0;

    // Create sites
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

    // Assign neighbours
    double J = 0.1;

    i_x = 0, i_y = 0, i_z = 0;

    // highest index in each direction
    int n_x = (int)std::floor(sizeX / periodX);
    int n_y = (int)std::floor(sizeY / periodY);
    int n_z = (int)std::floor(sizeZ / periodZ);

    // Particle index
    int i_p = 0;

    while (i_z * periodZ <= sizeZ)
    {
        while (i_y * periodY <= sizeY)
        {
            while (i_x * periodX <= sizeX)
            {

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
