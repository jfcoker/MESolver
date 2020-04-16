#pragma once

// Physical constansts
const double hbar = 6.58211899e-16; // eV, reduced plank constant
const double pi = 3.14159265; // dimensionless
const double e = 1.60217646e-19; // C, elementary charge

// Simulation parameters
const double F_Z = -1e-2; // V/Ang, Electric field strength in Z direction.
const double kBT = 0.02585; // eV, for T=300K
const double reorg = 0.01; // eV, reorganisation energy
const double threshold = 1; // disregard null spaces with associated singular values above this threshold. 