// This program should be capable of using SVD to solve the master equation describing the steady state transport of charge carriers in an array of molecules

#include "pch.h"
#include "consts.h"
#include "utility.h"
#include "IO.h"
#include "site.h"
#include "transporter.h"


// Simulation parameter labels
const char* label_F_z = "fieldZ"; // Electric field strength in Z direction.
const char* label_T = "temp"; // Temperature
const char* label_reorg = "reorg"; // Reorganisation energy
const char* label_periodic = "periodicZ"; // Information for periodic boundary conditions 

// Options
bool verbose = false;
bool periodic = false;
bool rescale = false;
double tolerance = 0.0;
double transE = 0.0;
std::vector<double> propagate;
transporter::PrecondForm form = transporter::PrecondForm::off;

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
        if (strstr(argv[i], "--precondition"))
        {
            form = transporter::PrecondForm::boltzmann;
            if (strstr(argv[i], "="))
            {
                char* substr = strchr(argv[i], '='); // Extract substring of all characters after and including '='
                ++substr; // Increment to get rid of '=' char
                if (strcmp(substr, "rateSum") == 0) form = transporter::PrecondForm::rateSum;
                else if (strcmp(substr, "boltzmannSquared") == 0) form = transporter::PrecondForm::boltzmannSquared;
            }
        }
        if (strcmp(argv[i], "--rescale") == 0) rescale = true;
        if (strstr(argv[i], "--tol="))
        {
            char* substr = strchr(argv[i], '=');
            tolerance = atof(++substr); // If not interpretable then atof will return 0.
        }
        if (strstr(argv[i], "--transE="))
        {
            char* substr = strchr(argv[i], '=');
            transE = atof(++substr); // If not interpretable then atof will return 0.
        }
        if (strstr(argv[i], "--propagate="))
        {
            char* substr = strchr(argv[i], '=');
            ++substr;
            char* token;
            char* next_token;
            double val;
            token = strtok_s(substr, ",", &next_token);
            while (token)
            {
                val = atof(token);
                if (val > 0.0) propagate.push_back(val);
                token = strtok_s(NULL, ",",&next_token);
            }
        }
    }

    // Print options
    std::cout << "Taking input from " << sim << ", " << xyz << ", " << edge << " ...\n";
    std::cout << "Preconditioning ";
    switch (form)
    {
        case transporter::PrecondForm::off: std::cout << "off\n"; break;
        case transporter::PrecondForm::boltzmann: std::cout << "on, form = boltzmann\n"; break;
        case transporter::PrecondForm::boltzmannSquared: std::cout << "on, form = boltzmannSquared\n"; break;
        case transporter::PrecondForm::rateSum: std::cout << "on, form = rateSum\n"; break;
    }
    std::cout << "Rescaling "; if (rescale) std::cout << "on\n"; else std::cout << "off\n";
    std::cout << "Singular value threshold = "; if (tolerance == 0.0) std::cout << "auto\n"; else std::cout << tolerance << "\n";
    if (!propagate.empty()) std::cout << "Testing time propagation of "; for (int i = 0; i < propagate.size(); i++) { std::cout << propagate[i] << "s "; }; std::cout << "\n";
    std::cout << "Verbosity "; if (verbose) std::cout << "high\n"; else std::cout << "low\n";


    std::cout << "\nReading simulation parameters...\n";
    const double F_z = ReadParameter(sim, label_F_z); // V/Ang
    const double temp = ReadParameter(sim, label_T); // K
    const double kBT = kB * temp; // eV
    const double reorg = ReadParameter(sim, label_reorg); // eV
    const double zsize = ReadParameterDefaultValue(sim, label_periodic, -1.0); // Ang
    if (zsize != -1.0) periodic = true;

    std::cout << "fieldZ (V/Ang) = " << F_z
        << "\ntemp (K) = " << temp
        << "\nreorg (eV) = " << reorg;
    if (periodic) std::cout << "\nPeriodic in z, zsize (Ang) = " << zsize;
    std::cout << "\n";

    std::cout << "\nCreating sites...\n";
    std::vector<site> allSites = CreateSites(xyz, edge);
    const size_t M = allSites.size(); // # sites
    for (int i = 0; i < M; i++)
        std::cout << allSites[i] << std::endl;

    // Create transporter object
    transporter transport(kBT, F_z, reorg, transE, periodic, zsize);

    if (form != transporter::PrecondForm::off) std::cout << "\nCreating preconditioned rate matrix A...\n";
    else std::cout << "\nCreating rate matrix A...\n";
    gsl_matrix* A = transport.CreateRateMatrix(allSites, form, rescale, verbose);

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

    if (tolerance == 0.0) // If tolerance has not been set manually, calculate it using machine epsilon.
        tolerance = std::numeric_limits<double>::epsilon() * 
                    std::max(std::max(gsl_matrix_max(U),std::abs(gsl_matrix_min(U))),
                             std::max(gsl_matrix_max(V), std::abs(gsl_matrix_min(V))));

    std::cout << "\n\nDisregarding singular values greater than threshold = " << tolerance << "\n";
    std::cout << "Printing possible solutions\n";
    gsl_matrix* cleanA = transport.CreateRateMatrix(allSites, transporter::PrecondForm::off, false, false); // Need the non-conditioned, non-scaled rate matrix for time propogation and velocity calculations.
    int solnum = 0;
    for (int i = 0; i < S->size; i++)
    {
        double sval = gsl_vector_get(S, i);
        if (sval <= tolerance)
        {
            ++solnum;
            std::cout << "\n\nPossible solution " << solnum << " : singular value = " << sval << "\n";

            gsl_matrix_get_col(Q, V, i);
            if (form != transporter::PrecondForm::off) {
                std::cout << "\nConditioned densities\n";
                printVector(Q);

                // Reverse preconditioning
                for (int j = 0; j < Q->size; j++)      
                    gsl_vector_set(Q, j, gsl_vector_get(Q, j) * transport.PrecondFactor(&allSites[j], form));

                // Renormalise so squared values add to 1
                normalise(Q);
            }

            // If largest value is negative then flip all signs.
            if (abs(gsl_vector_min(Q)) > gsl_vector_max(Q)) gsl_vector_scale(Q, -1.0);

            std::cout << "\nOccupation densities\n";
            for (int j = 0; j < Q->size; j++)
                allSites[j].occProb = gsl_vector_get(Q, j);
            printOccProbs(allSites, 2);

            // Propagate densities in time (Can be useful to check if the solution is steady state).
            if (!propagate.empty())
            {
                std::cout << "\nTime propagation\n";
                gsl_matrix* Axt = gsl_matrix_alloc(M,M);
                gsl_matrix* expAxt = gsl_matrix_alloc(M, M);
                gsl_vector* Qt = gsl_vector_alloc(M);
                for (int i = 0; i < propagate.size(); i++)
                {
                    gsl_matrix_memcpy(Axt, cleanA);
                    gsl_matrix_scale(Axt, propagate[i]);
                    gsl_linalg_exponential_ss(Axt, expAxt, GSL_PREC_DOUBLE);
                    gsl_blas_dgemv(CblasNoTrans, 1.0, expAxt, Q, 0.0, Qt);
                    std::cout << "\nP( " << propagate[i] << "s ) = \n";
                    printVector(Qt);
                }
                gsl_matrix_free(Axt);
                gsl_matrix_free(expAxt);
                gsl_vector_free(Qt);
            }

            double v_z = transport.velocity_z(allSites, cleanA);
            std::cout << "\nvelocity_z = " << v_z << "\n";
            if (F_z != 0.0) std::cout << "\nmobility = " << v_z / F_z << "\n";
            
            //if (verbose)
            //{
            //    std::cout << "\nCHECK: is P a solution?\nA x P =\n";
            //    gsl_blas_dgemv(CblasNoTrans, 1.0, A, P, 0.0, AxP);
            //    printVector(AxP);
            //}

        }
    }

    // Free memory
    gsl_matrix_free(cleanA);
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