// This program should be capable of using SVD to solve the master equation describing the steady state transport of charge carriers in an array of molecules

#include "pch.h"
#include "site.h"
#include "simConsts.h"

void printMatrix(gsl_matrix* m);

int main()
{

    std::cout << "Creating sites...\n";
    std::vector<site> allSites = CreateSites();
    const int M = allSites.size(); // # sites
    for (int i = 0; i < M; i++)
        std::cout << allSites[i] << std::endl;

    std::cout << "\nCalculating transfer rates...\n";
    for (int i = 0; i < M; i++)
        for (int f = 0; f < M; f++)
        {
            double rate = allSites[i].Rate(&allSites[f]);
            if (rate) // if non-zero
                std::cout << "Rate " << i << "->" << f << " = " << rate << std::endl;
        }

    std::cout << "\nCreating rate matrix A\n";
    gsl_matrix* A = gsl_matrix_alloc(M, M);
    gsl_matrix_set_zero(A);
    for (int i = 0; i < M; i++)
        for (int f = 0; f < M; f++)
        {
            if (i == f)
            {
                double sum = 0.0;
                for (int k = 0; k < M; k++)
                    if (i != k)
                        sum -= allSites[i].Rate(&allSites[k]);
                gsl_matrix_set(A, i, f, sum);
            }

            else
                gsl_matrix_set(A, i, f, allSites[f].Rate(&allSites[i])); // May need to switch i and f?
        }

    printMatrix(A);

    // Free memory
    gsl_matrix_free(A);

    return 0;
}

void printMatrix(gsl_matrix* m)
{

    std::stringstream sstream;
    sstream.setf(std::ios::scientific);
    sstream.precision(1);

    for (int i = 0; i < m->size1; i++) 
    {
        for (int j = 0; j < m->size2; j++)
        {   
            double el = gsl_matrix_get(m, i, j);
            if (el)
                sstream << std::setw(9) << std::left << el;
            else
                sstream << std::setw(9) << std::left << "0";
        }
        std::cout << sstream.str() << "\n";
        sstream.str("");

        
    }

}


std::vector<site> CreateSites()
{
    std::vector<site> sites;

    // For now, hardcode some information about the array
    double E = 0.0;
    double sizeX = 10.0, sizeY = 10.0, sizeZ = 20.0;
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
