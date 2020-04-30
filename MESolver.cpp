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

// Output formatting options
bool verbose = false;
bool precondition = false;
bool rescale = false;
bool tol_auto = true;
double tolerance = 0.0;
double transE = 0.0;

int main(int argc, char* argv[])
{
    //Parse command line parameters.
    if (argc < 4) {
        std::cout << "*** ERROR ***: Expect at least three input files: .sim, .xyz, .edge\n";
        exit(-1);
    }
    char sim[128], xyz[128], edge[128], occ[128];
    for (int i = 1; i < argc; i++) {
        
        // Input files
        if (strstr(argv[i], ".sim"))  strcpy_s(sim, argv[i]);
        if (strstr(argv[i], ".xyz"))  strcpy_s(xyz, argv[i]);
        if (strstr(argv[i], ".edge")) strcpy_s(edge, argv[i]);
        if (strstr(argv[i], ".occ")) strcpy_s(occ, argv[i]);

        // Options
        if (strcmp(argv[i], "-v") == 0) verbose = true;
        if (strcmp(argv[i], "--precondition") == 0) precondition = true;
        if (strcmp(argv[i], "--rescale") == 0) rescale = true;
        if (strstr(argv[i], "--tol"))
        {
            char* substr = strchr(argv[i], '='); // Extract substring of all characters after and including '='
            tolerance = atof(++substr); // Increment to get rid of '=' char, then attempt to convert to double
            if (tolerance) tol_auto = false; // If parsed value is non-zero and interpretable, then use as tolerance
        }
        if (strstr(argv[i], "--transE"))
        {
            char* substr = strchr(argv[i], '='); // Extract substring of all characters after and including '='
            transE = atof(++substr); // Increment to get rid of '=' char, then attempt to convert to double. If it is not interpretable then E0 is set to 0.
        }
    }

    // Print options
    std::cout << "Taking input from " << sim << ", " << xyz << ", " << edge << " ...\n";
    std::cout << "Preconditioning "; if (precondition) std::cout << "on\n"; else std::cout << "off\n";
    std::cout << "Rescaling "; if (rescale) std::cout << "on\n"; else std::cout << "off\n";
    std::cout << "Singular value threshold = "; if (tol_auto) std::cout << "auto\n"; else std::cout << tolerance << "\n";
    std::cout << "Verbosity "; if (verbose) std::cout << "high\n"; else std::cout << "low\n";


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

    if (precondition) std::cout << "\nCreating preconditioned rate matrix A...\n";
    else std::cout << "\nCreating rate matrix A...\n";

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
                        sum += allSites[i].Rate(&allSites[k], F_z, kBT, reorg) * allSites[i].PrecondFactor(F_z, kBT, transE, precondition);
                el = -sum;
            }
            else
            {
                el = allSites[f].Rate(&allSites[i], F_z, kBT, reorg) * allSites[f].PrecondFactor(F_z, kBT, transE, precondition); // May need to switch i and f?
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


 
    if (rescale)
    {
        std::cout << "\nTo reduce precision errors, rescale A by 1e-" << highestO << "\n";
        gsl_matrix_scale(A, pow(10, -highestO));

        if (verbose)
        {
            std::cout << "Reduced A = \n";
            printMatrix(A);
        }
    }

    std::cout << "\nSolving ME using SVD...\n";
    gsl_matrix* U = gsl_matrix_alloc(M, M);
    gsl_matrix* V = gsl_matrix_alloc(M, M);
    gsl_matrix* VT = gsl_matrix_alloc(M, M);
    gsl_vector* S = gsl_vector_alloc(M);
    gsl_matrix* Sigma = gsl_matrix_alloc(M, M);
    gsl_vector* work = gsl_vector_alloc(M);
    gsl_vector* Q = gsl_vector_alloc(M);
    gsl_matrix* UxSigma = gsl_matrix_alloc(M, M);
    gsl_matrix* UxSigmaxVT = gsl_matrix_alloc(M, M);
    gsl_vector* AxP = gsl_vector_alloc(M);

    gsl_matrix_memcpy(U, A);
    gsl_linalg_SV_decomp(U, V, S, work);
    gsl_matrix_transpose_memcpy(VT, V);

    //Create Sigma matrix
    gsl_matrix_set_zero(Sigma);
    for (int i = 0; i < M; i++)
        gsl_matrix_set(Sigma, i, i, gsl_vector_get(S, i));

    if (verbose)
    {
        std::cout << "\nU = \n";
        printMatrix(U);
    }

    std::cout << "\nSingular values = \n";
    printVector(S, false);

    if (verbose)
    {
        std::cout << "\nV = \n";
        printMatrix(V);


        std::cout << "\nV^T = \n";
        printMatrix(VT);


        std::cout << "\nCHECK: does U x Sigma x VT = A?\nU x Sigma x VT =\n";
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, U, Sigma, 0.0, UxSigma);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, UxSigma, VT, 0.0, UxSigmaxVT);
        printMatrix(UxSigmaxVT);
    }

    if (tol_auto)
        tolerance = std::numeric_limits<double>::epsilon() * 
                    std::max(std::max(gsl_matrix_max(U),std::abs(gsl_matrix_min(U))),
                             std::max(gsl_matrix_max(V), std::abs(gsl_matrix_min(V))));

    std::cout << "\n\nDisregarding singular values greater than threshold = " << tolerance << "\n";
    std::cout << "Printing possible solutions\n";
    int solnum = 0;
    for (int i = 0; i < S->size; i++)
    {
        double sval = gsl_vector_get(S, i);
        if (sval <= tolerance)
        {
            ++solnum;
            std::cout << "\n\nPossible solution " << solnum << " : singular value = " << sval << "\n";

            gsl_matrix_get_col(Q, V, i);
            if (precondition) {
                std::cout << "\nConditioned densities\n";
                printVector(Q);

                // Reverse preconditioning
                for (int j = 0; j < Q->size; j++)
                    gsl_vector_set(Q, j, gsl_vector_get(Q, j) * allSites[j].PrecondFactor(F_z, kBT,transE, precondition));

                // Renormalise so squared values add to 1
                normalise(Q);
            }

            // If largest value is negative then flip all signs.
            if (abs(gsl_vector_min(Q)) > gsl_vector_max(Q)) gsl_vector_scale(Q, -1.0);

            std::cout << "\nOccupation densities\n";
            for (int j = 0; j < Q->size; j++)
                allSites[j].occProb = gsl_vector_get(Q, j);
            printOccProbs(allSites, 2);

            //if (verbose)
            //{
            //    std::cout << "\nCHECK: is P a solution?\nA x P =\n";
            //    gsl_blas_dgemv(CblasNoTrans, 1.0, A, P, 0.0, AxP);
            //    printVector(AxP);
            //}

        }
    }

    // Free memory
    gsl_matrix_free(A);
    gsl_matrix_free(U);
    gsl_matrix_free(V);
    gsl_matrix_free(VT);
    gsl_vector_free(S);
    gsl_matrix_free(Sigma);
    gsl_vector_free(work);
    gsl_vector_free(Q);
    gsl_matrix_free(UxSigma);
    gsl_matrix_free(UxSigmaxVT);
    gsl_vector_free(AxP);

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