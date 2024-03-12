#ifndef PARAMETERS
#define PARAMETERS

extern const double Omega;      //total carrying capacity
extern const double h;          // lattice size
extern const int N;             // number of voxels
extern const int R;             // number of reactions

extern const double Tmax;       // maximum simulation time

extern const double Du;         // diffusion coefficient of u; Du>Dv
extern const double Dv;         // diffusion coefficient of v
extern const double F;          // Gray-Scott model parameter
extern const double k;          // Gray-Scott model parameter
extern const double rho;        // Gray-Scott model parameter

#endif