// This program should be capable of using SVD to solve the master equation describing the steady state transport of charge carriers in an array of molecules

#include "pch.h"
#include "utility.h"
#include "site.h"
#include "simConsts.h"

struct ArraySize { size_t X, Y, Z;};
ArraySize CmdLineParser(int argc, char* argv[])
{
    ArraySize rtn;
    if (argc == 4)
    {
        rtn.X = std::atoi(argv[1]);
        rtn.Y = std::atoi(argv[2]);
        rtn.Z = std::atoi(argv[3]);
    }

    // If arguments cannot be converted to int, atoi will return 0.
    if (argc != 4 || rtn.X == 0 || rtn.Y == 0 || rtn.Z == 0) 
    {
        // Hardcoded default values
        rtn.X = 3;
        rtn.Y = 3;
        rtn.Z = 15;
    }

    return rtn;
}

int main(int arg, char* argv[])
{
    std::cout << "Simulation parameters"
        << "\nE-Field Strength, F_z (V/Ang) = " << F_Z
        << "\nk_B * Temp (eV) = " << kBT
        << "\nreorg energy (eV) = " << reorg
        << "\n\n";

    ArraySize sz = CmdLineParser(arg, argv);
    std::cout << "Creating sites for "<< sz.X << "x" << sz.Y << "x" << sz.Z << " array...\n";
    std::cout << "Every 5th and 7th site is a trap\n";
    std::vector<site> allSites = CreateSites(sz.X, sz.Y, sz.Z);
    const size_t M = allSites.size(); // # sites
    for (int i = 0; i < M; i++)
        std::cout << allSites[i] << std::endl;

    std::cout << "\nCreating rate matrix A...\n";
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
                gsl_matrix_set(A, i, f, allSites[f].Rate(&allSites[i])); // May need to switch i and f (will n make a difference if rates are non-symmetric)?
        }

    //std::cout << "\nA = \n";
    //printMatrix(A);

    std::cout << "\nSolving ME using SVD...\n";
    gsl_matrix* V = gsl_matrix_alloc(M, M);
    gsl_vector* S = gsl_vector_alloc(M);
    gsl_vector* work = gsl_vector_alloc(M);
    gsl_vector* P = gsl_vector_alloc(M);
    gsl_linalg_SV_decomp(A, V, S, work);

    std::cout << "\nSingular values = \n";
    printVector(S, false);
    std::cout << "\nDisregarding singular values greater than threshold = " << threshold << "\n";
    std::cout << "\nPrinting possible solutions\n\n";

    for (int i = 0; i < S->size; i++)
    {
        int solnum = 1;
        if (gsl_vector_get(S, i) <= threshold)
        {
            std::cout << "Solution " << solnum++ << ": \n";
            gsl_matrix_get_col(P, V, i);
            for (int i = 0; i < P->size; i++) 
                allSites[i].occProb = gsl_vector_get(P, i);
            printOccProbs(allSites);
        }
    }


    //######## Uncomment to be able to interactively select columns from V ##############
    //std::cout << "\nInput column index of V to view occupation probabilities (non-int to exit)\n";
    //std::cout << "\nMax is " << V->size2 - 1 << " : ";

    //size_t col;
    //while (std::cin >> col)
    //{
    //    if (col >= V->size2)
    //    {
    //        std::cout << "\nIndex too large.\n";
    //    }
    //    else
    //    {
    //        std::cout << "\nP =\n";
    //        gsl_matrix_get_col(P, V, col);
    //        for (int i = 0; i < P->size; i++) 
    //            allSites[i].occProb = gsl_vector_get(P, i);
    //        printOccProbs(allSites);

    //    }

    //    std::cout << "\n\nInput column index of V to view occupation probabilities (non-int to exit)\n"; 
    //    std::cout << "\nMax is " << V->size2 - 1 << " : ";

    //}
    //######################################################################################

    // Free memory
    gsl_matrix_free(A);
    gsl_matrix_free(V);
    gsl_vector_free(S);
    gsl_vector_free(work);
    gsl_vector_free(P);

    return 0;
}

std::vector<site> CreateSites(size_t X, size_t Y, size_t Z)
{
    std::vector<site> sites;

    // For now, hardcode some information about the array
    double E = 0.0;
    double Etrap = -0.05;
    double periodX = 10.0, periodY = 10.0, periodZ = 10.0;
    
    double sizeX = (X - 1)*periodX, sizeY = (Y - 1)*periodY, sizeZ = (Z - 1)*periodZ;

    // Current site index
    int i_p = 0;

    // axis indices
    int i_x = 0, i_y = 0, i_z = 0;

    // Create sites. Currently just using a cubic array.
    while (i_z * periodZ <= sizeZ)
    {
        while (i_y * periodY <= sizeY)
        {
            while (i_x * periodX <= sizeX)
            {
                // Make every 5th and 7th site a trap.
                double tmpD;
                bool tmpB;
                if (i_p % 5 == 0 || i_p % 7 == 0) {tmpD = Etrap; tmpB = true;}
                else { tmpD = E; tmpB = false; }
                sites.push_back(site(i_x * periodX, i_y * periodY, i_z * periodZ, tmpD,tmpB));
                i_p++;
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

    i_p = 0; i_x = 0, i_y = 0, i_z = 0;

    // highest index in each direction
    size_t n_x = (int)std::floor(sizeX / periodX);
    size_t n_y = (int)std::floor(sizeY / periodY);
    size_t n_z = (int)std::floor(sizeZ / periodZ);

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
