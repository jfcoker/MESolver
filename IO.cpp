#include "pch.h"
#include "IO.h"

void open(char* filename, std::ifstream& in) {
    in.open(filename);
    if (!in) {
        std::cout << "***ERROR***: Unable to open " << filename << std::endl;
        exit(-1);
    }
}

double ReadParameter(char* filename, std::string name) {
    std::ifstream in;
    open(filename, in);

    std::string word;
    double value;
    while (in) {
        in >> word;
        if (word == name) 
        {
            in >> value;
            if (in.fail()) 
            {
                std::cout << "***ERROR***: Could not convert parameter value to double.\n";
                exit(-1);
            }
            else 
            {
                return value;
            }
        }
    }
    std::cout << "***ERROR***: Reached end-of-file. Could not find " << name << " in " << filename << "\n";
    in.close();
    exit(-1);
}