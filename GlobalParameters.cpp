#include "Declarations.h"
#include"Headers/GlobalParameters.h"

extern const double Omega(25000);         //total system size
extern const int N(50);                 // number of voxels
extern const double h(1.0/(double)N);   // lattice size
extern const int R(6);                  // number of reactions

extern const double Tmax(500.0);          // maximum simulation time

extern const double Du(1.0e-5);        // diffusion coefficient of u; Du>Dv
extern const double Dv(1.0e-6);         // diffusion coefficient of v
extern const double F(0.053);           // Gray-Scott model parameter
extern const double k(0.0516);          // Gray-Scott model parameter
extern const double rho(1.0);           // Gray-Scott model parameter
