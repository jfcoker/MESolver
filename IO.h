#pragma once
#include "pch.h"

double ReadParameter(char* filename, std::string name);
double ReadParameterDefaultValue(char* filename, std::string name, double def);
void open(char* filename, std::ifstream& filestream);