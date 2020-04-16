#include "pch.h"
#include "utility.h"

void printMatrix(gsl_matrix* m)
{
    std::stringstream sstream;
    sstream.setf(std::ios::scientific);
    sstream.precision(2);

    for (int i = 0; i < m->size1; i++)
    {
        for (int j = 0; j < m->size2; j++)
        {
            double el = gsl_matrix_get(m, i, j);
            if (el)
                sstream << std::setw(10) << std::left << el;
            else
                sstream << std::setw(10) << std::left << "0";
        }
        std::cout << sstream.str() << "\n";
        sstream.str("");
    }
}

void printVector(gsl_vector* v, bool horizontal)
{
    std::stringstream sstream;
    sstream.setf(std::ios::scientific);
    sstream.precision(2);

    if (horizontal)
    {
        sstream << "( ";
        for (size_t i = 0; i < v->size; i++)
        {
            double el = gsl_vector_get(v, i);
            if (el) sstream << el; else sstream << "0";
            if (i < v->size - 1) sstream << ", ";
        }
        sstream << " )\n";
    }
    else 
    {
        for (size_t i = 0; i < v->size; i++)
        {
            double el = gsl_vector_get(v, i);
            if (el) sstream << el; else sstream << "0";
            sstream << "\n";
        }
    }

    std::cout << sstream.str();
}

void printDiagonal(gsl_matrix* m, bool horizontal)
{
    size_t minSize = std::min(m->size1, m->size2);
    gsl_vector* V = gsl_vector_alloc(minSize);

    for (size_t i = 0; i < minSize; i++)
        gsl_vector_set(V, i, gsl_matrix_get(m, i, i));

    printVector(V, horizontal);
    gsl_vector_free(V);
}

void printOccProbs(std::vector<site>& sites, double period)
{
    std::stringstream sstream;
    sstream.setf(std::ios::scientific);
    sstream.precision(2);

    //for (size_t i = 0; i < v->size; i++)
    //{
    //    double el = gsl_vector_get(v, i);
    //    if (el) sstream << el; else sstream << "0";
    //    if (i < v->size - 1) sstream << ", ";
    //}
    //sstream << "\n";


    double el;
    el = sites[0].occProb;
    if (el) sstream << el; else sstream << "0";

    for (int i = 1; i < sites.size(); i++)
    {
        if (std::abs(sites[i].pos.Z - sites[i - 1].pos.Z) > period / 2)
        {
            sstream << "\n";
            std::cout << sstream.str();
            sstream.str("");

        } else
        if (std::abs(sites[i].pos.Y - sites[i - 1].pos.Y) > period / 2)
        {
            sstream << "          ";
        }


        sstream << std::setw(10) << std::left;
        el = sites[i].occProb;
        if (el) sstream << el; else sstream << "0";

    }

    std::cout << sstream.str();

}