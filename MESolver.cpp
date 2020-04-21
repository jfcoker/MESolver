// This program should be capable of using SVD to solve the master equation describing the steady state transport of charge carriers in an array of molecules

#include "pch.h"
#include "consts.h"
#include "utility.h"
#include "IO.h"
#include "site.h"

// Simulation parameter labels
const char* label_F_z = "fieldZ"; // Electric field strength in Z direction.
const char* label_T = "temp"; // Temperature
const char* label_reorg = "reorg"; // Reorganisation energy

// Other
const double threshold = 1; // disregard null spaces with associated singular values above this threshold.

int main(int argc, char* argv[])
{
    //Parse command line parameters.
    if (argc < 4) {
        std::cout << "*** ERROR ***: Expect at least three input files: .sim, .xyz, .edge\n";
        exit(-1);
    }
    char sim[128], xyz[128], edge[128], occ[128];
    for (int i = 1; i < argc; i++) {
        if (strstr(argv[i], ".sim"))  strcpy_s(sim, argv[i]);
        if (strstr(argv[i], ".xyz"))  strcpy_s(xyz, argv[i]);
        if (strstr(argv[i], ".edge")) strcpy_s(edge, argv[i]);
        if (strstr(argv[i], ".occ")) strcpy_s(occ, argv[i]);
    }

    std::cout << "Taking input from " << sim << ", " << xyz << ", " << edge << " ...\n";

    std::cout << "\nReading simulation parameters...\n";
    const double F_z = ReadParameter(sim, label_F_z); // V/Ang
    const double temp = ReadParameter(sim, label_T); // K
    const double reorg = ReadParameter(sim, label_reorg); // eV
    const double kBT = kB * temp; // eV
    std::cout << "fieldZ (V/Ang) = " << F_z
        << "\ntemp (K) = " << temp
        << "\nreorg (eV) = " << reorg
        << "\n";

    std::cout << "\nCreating sites...\n";
    std::vector<site> allSites = CreateSites(xyz, edge);
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
                        sum += allSites[i].Rate(&allSites[k], F_z, kBT, reorg);
                gsl_matrix_set(A, i, f, -sum);
            }

            else
                gsl_matrix_set(A, i, f, allSites[f].Rate(&allSites[i], F_z, kBT, reorg)); // May need to switch i and f (will only make a difference if rates are non-symmetric)?
        }
    printMatrix(A);

    std::cout << "\nSolving ME using SVD...\n";
    gsl_matrix* V = gsl_matrix_alloc(M, M);
    gsl_vector* S = gsl_vector_alloc(M);
    gsl_vector* work = gsl_vector_alloc(M);
    gsl_vector* P = gsl_vector_alloc(M);
    //gsl_linalg_SV_decomp(A, V, S, work);
    gsl_linalg_SV_decomp_jacobi(A, V, S);

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

std::vector<site> CreateSites(char* XYZfile, char* EDGEfile)
{
    std::vector<site> sites;

    // Use contents of .xyz file to populate site vector
    // Columns should be formatted as: x (float), y (float), z (float), site type (char, unused), site energy (float)
    std::ifstream in;
    open(XYZfile, in);
    std::string line;
    
    int linenum = 0;
    while (std::getline(in, line) )
    {
        linenum++;
        std::stringstream linestream(line);
        double x, y, z, E;
        char s;
        linestream >> x >> y >> z >> s >> E;

        if (linestream.fail())
        {
            std::cout << "***ERROR***: Unexpected formatting of .xyz file at line " << linenum << ". Expected 5 columns with value types; float, float, float, char, float.\n";
            in.close();
            exit(-1);
        }

        sites.push_back(site(x, y, z, E));
    }

    in.close();

    // Use contents of .edge file to set interacting neighbours
    // Columns should be formatted as: site 1 (int), site 2 (int), J (float)
    open(EDGEfile, in);

    linenum = 0;
    while (std::getline(in, line))
    {
        linenum++;
        std::stringstream linestream(line);
        int s1, s2;
        double J;
        linestream >> s1 >> s2 >> J;

        if (linestream.fail())
        {
            std::cout << "***ERROR***: Unexpected formatting of .edge file at line " << linenum << ". Expected 3 columns with value types; int, int, float.\n";
            in.close();
            exit(-1);
        }
        
        if (s1 >= sites.size() || s2 >= sites.size()) 
        { 
            std::cout << "***ERROR***: Site index out of range in .edge file at line " << linenum << "\n";
            in.close();
            exit(-1);
        }

        sites[s1].addNeighbour(&sites[s2], J);

    }

    in.close();

    return sites;

}