#include "fargo3d.h"

void _CondInit(real epsilon) {
  
  int i,j,k;
  real r, omega;
  
  real *rho  = Density->field_cpu;
  real *cs   = Energy->field_cpu;
  real *vphi = Vx->field_cpu;
  real *vr   = Vy->field_cpu;
  
  real rhog, rhod;
  real vk;
  
  i = j = k = 0;
  
  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx+2*NGHX; i++) {
	
	r     = Ymed(j);
	omega = sqrt(G*MSTAR/r/r/r);                       //Keplerian frequency
	rhog  = SIGMA0*pow(r/R0,-SIGMASLOPE);              //Gas surface density
        rhod  = rhog*epsilon;                              //Dust surface density

	if (Fluidtype == GAS) {
	  rho[l]   = rhog;
	  vphi[l]  = omega*r*sqrt(1.0 + pow(ASPECTRATIO,2)*pow(r/R0,2*FLARINGINDEX)*
				  (2.0*FLARINGINDEX - 1.0 - SIGMASLOPE));
	  vr[l]    = 0.0;
	  cs[l]    = ASPECTRATIO*pow(r/R0,FLARINGINDEX)*omega*r;
	}
	
	if (Fluidtype == DUST) {
	  rho[l]  = rhod;
	  vphi[l] = omega*r;
	  vr[l]   = 0.0;
	  cs[l]   = 0.0;
	}
	
	vphi[l] -= OMEGAFRAME*r;
	
      }
    }
  }
}

void CondInit() {
  float InvStokes[NFLUIDS] = { 0.0, INVSTOKES1, INVSTOKES2, INVSTOKES3 };
  float Epsilons[NFLUIDS] = { 1.0, EPSILON1, EPSILON2, EPSILON3 };
  
  int id_gas = 0;
  //We first create the gaseous fluid and store it in the array Fluids[]
  Fluids[id_gas] = CreateFluid("gas",GAS);

  //We now select the fluid
  SelectFluid(id_gas);

  //and fill its fields
  _CondInit(Epsilons[id_gas]);

  //We repeat the process for the dust fluids
  char dust_name[MAXNAMELENGTH];
  int id_dust;

  for(id_dust = 1; id_dust<NFLUIDS; id_dust++) {
    sprintf(dust_name,"dust%d",id_dust); //We assign different names to the dust fluids

    Fluids[id_dust] = CreateFluid(dust_name, DUST);
    SelectFluid(id_dust);
    _CondInit(Epsilons[id_dust]);

    /*We now fill the collision matrix (Feedback from dust included)
     Note: ColRate() moves the collision matrix to the device.
     If feedback=NO, gas does not feel the drag force.*/
  
    // Normalise the inverse stokes parameter to be in units of
    // the initial gas density at R=R0.  This will then be used to
    // rescale to the current cell gas density in collision_kernel.h
    ColRate(InvStokes[id_dust]/SIGMA0, id_gas, 1, GASFEELSDUSTDRAG);
  }
}
