#include "idefix.hpp"
#include "setup.hpp"
#include <cmath>

// Default constructor
real v_from_lM(real lMx, real c2x);
real ShockCompressionFactor (real M1x);
real FindEntropy(real zx, real vx, real c2x);
real Phi(real lMx, real c2x, real zx, real Sx, real Ah);
real HeatingCooling (real zx, real px, real Ah);
real CoolingCutoff(real Sx, real Sx_rpns);
void static ComputeBondi(real radx, real& rho, real& pgas, real& vx);
real LinearInterpolation(std::vector<real>& x, std::vector<real>& y, real x_interp);
real *RK4_RHS (real t, int n, real u[], real Ah);
real *RK4Update (real t0, int m, real u0[], real dt, real *f (real t, int m, real u[], real Ah), real Ah); 
real CalcStationaryShock (real M_rsh0x, real c2_rsh0x, real r_shx, real Ah);
void UserdefBoundary(Hydro *hydro, int dir, BoundarySide side, real t);
void MySourceTerm(Hydro *hydro, const real t, const real dtin);
void CoarsenFunction(DataBlock &data);

real Aheat, Acool, rsh; 
real gammaGlob, gm1, gm1i;

real M1       = 5.0; // Pre-shock Mach number
real rsh0     = ONE_F; // Initial shock position for solution with cooling
real rpns     = 0.4; // initial PNS radius
real epsbar   = 0.0; // Dimensionaless dissociation energy (eps = epsbar*v1^2 = 0.5*epsbar*vff^2)
real kappa    = ShockCompressionFactor(M1);  
real v1       = -pow((1.0 / rsh0 + 3.0 * 0.006431) / (0.5 + 3.0 / pow(M1,2)), 0.5); 
real c1       = -v1 / M1;    // pre-shock sound speed
real v_rsh0   = v1 / kappa;  // post-shock velocity
real c2_rsh0  = (pow(v1, 2) - pow(v_rsh0, 2) + 6 * pow(c1, 2) - 2 * epsbar) / 6.0; 
real M_rsh0   = -v_rsh0 / pow(c2_rsh0, 0.5);
real rho_rsh0 = v1 / v_rsh0; // post-shock density
real p_rsh0   = c2_rsh0 * rho_rsh0 / (4.0/3.0);      // post-shock pressure

real v_rsh    = v_rsh0;
real c2_rsh   = c2_rsh0;
real M_rsh    = M_rsh0;
real rho_rsh  = rho_rsh0;
real p_rsh    = p_rsh0;

real Mmin    = 0.0; // Minimum Mach number at PNS surface r*
real s_min   = 0.0; // Use s_rpns = 0.0 
real s_rpns  = 0.0; // Use s_rpns = 0.0 

real v_rpns, rho_rpns, p_rpns, c_rpns, M_rpns;

int const Nmax = 300; // number of grid cells for initial data calcution using RK4.

std::vector<real> z_arr, M_arr, c_arr;

//----------------------------------------------------------------------------------------

// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {

  gammaGlob = data.hydro->eos->GetGamma();
  gm1  = gammaGlob - ONE_F; 
  gm1i = ONE_F / gm1; 

  s_min  = input.Get<real>("Setup", "s_min",  ZERO_F);
  Mmin  = input.Get<real>("Setup", "Mmin",   ZERO_F);
  Acool = input.Get<real>("Setup", "Acool",  ZERO_F);
  Aheat = input.Get<real>("Setup", "Aheat",  ZERO_F);
  rsh   = input.Get<real>("Setup", "Rshock", ZERO_F);

  data.hydro->EnrollUserDefBoundary(&UserdefBoundary);
  data.hydro->EnrollUserSourceTerm(&MySourceTerm);
  if(data.haveGridCoarsening) {
    data.EnrollGridCoarseningLevels(&CoarsenFunction);
  }
}
//----------------------------------------------------------------------------------------

void Setup::InitFlow(DataBlock &data) {

  DataBlockHost d(data); // Create a host copy

  rpns = CalcStationaryShock (M_rsh0, c2_rsh0, rsh0, 0.0); // assuming s_rpns = 0.0

  v_rpns   = - M_arr[Nmax-1] * c_arr[Nmax-1]; // velocity at r=rpns
  c_rpns   = c_arr[Nmax-1]; // sound speed at r=rpns;
  rho_rpns = v1 / v_rpns / pow(rpns, 2); // density at r=rpns
  p_rpns   = rho_rpns * pow(c_rpns, 2) / gammaGlob; // pressure at r=rpns

  s_rpns   = gm1i * log(p_rpns / pow(rho_rpns, gammaGlob)); // entropy at r=rpns
  
  std::cout << "PNS radius: " << rpns << std::endl;
  std::cout << "Entropy at PNS surface: " << s_rpns << std::endl;
  std::cout << "Mach at PNS: " << M_arr[Nmax-1] << std::endl;
  std::cout << "vr at PNS: " << v_rpns << std::endl;
  std::cout << "cs at PNS: " << c_rpns << std::endl;

  s_rpns   = s_min; // assign entropy

  // heating solution
  if (FABS(rsh - rsh0) > 1.0e-12) {
    // heating/cooling solution
    real v_tmp, rho_tmp, p_tmp, c_tmp;
    ComputeBondi(rsh, rho_tmp, p_tmp, v_tmp);
    c_tmp   = std::sqrt(gammaGlob * p_tmp / rho_tmp);
    kappa   = ShockCompressionFactor(- v_tmp / c_tmp);  
    v_rsh   = v_tmp / kappa;  // post-shock velocity
    c2_rsh  = (pow(v_tmp, 2) - pow(v_rsh, 2) + 6 * pow(c_tmp, 2) - 2 * epsbar) / 6.0; 
    M_rsh   = - v_rsh / pow(c2_rsh, 0.5);
    rho_rsh = v1 / (v_rsh * pow(rsh,2)); // post-shock density
    p_rsh   = c2_rsh * rho_rsh0 / gammaGlob;  // post-shock pressure

    std::cout << "Post-shock Mach number: " << M_rsh << std::endl;
    std::cout << "Post-shock density: " << rho_rsh << std::endl;
    std::cout << "Post-shock pressure: " << p_rsh << std::endl;
    std::cout << "Post-shock sound speed: " << std::sqrt(c2_rsh) << std::endl;
    std::cout << "Post-shock velocity: " << v_rsh << std::endl;
  
    std::cout << "Aheat = " << Aheat << std::endl;

    z_arr.clear(); M_arr.clear(); c_arr.clear();
    rpns = CalcStationaryShock (M_rsh, c2_rsh, rsh, Aheat); 

    v_rpns   = - M_arr[Nmax-1] * c_arr[Nmax-1]; // velocity at r=rpns
    c_rpns   = c_arr[Nmax-1]; // sound speed at r=rpns;
    rho_rpns = v1 / v_rpns / pow(rpns, 2); // density at r=rpns
    p_rpns   = rho_rpns * pow(c_rpns, 2) / gammaGlob; // pressure at r=rpns

    std::cout << "PNS radius: " << rpns << std::endl;
    std::cout << "Mach at PNS: " << M_arr[Nmax-1] << std::endl;
  }

  real dEdt = ZERO_F; // total energy change due to cooling
  for(int i = 0; i < d.np_tot[IDIR]; i++) {
    real r  = d.x[IDIR](i); // cell center
    real ri = r - HALF_F * d.dx[IDIR](i); // left cell boundary
    if (r < rsh && r > rpns) {
      //std::cout << i << " " << d.dx[IDIR](i) / r << " " << d.dx[IDIR](i) << std::endl;
      // Calculate cooling in the post-shock region
      real cx   = LinearInterpolation(z_arr, c_arr, ri);
      real vx   = - cx * LinearInterpolation(z_arr, M_arr, ri);
      real rhox = FABS(v1 / (vx * pow(ri, 2)));
      real px   = rhox * pow(cx, 2) / gammaGlob;
      real sx   = gm1i * log(FABS(px / pow(rhox, gammaGlob)));
      dEdt += HeatingCooling(ri, px, 0.0) *
              rhox * 4.0 * M_PI * pow(r, 2) * d.dx[IDIR](i);
    }
  }
  // Adjust cooling rate
  real Acool_correction = (4.0 * M_PI * v1 / rpns) / dEdt;
  std::cout << "Cooling correction factor: " << Acool_correction << std::endl;
  Acool *= Acool_correction; 

  // interpolate the initial solution to the grid
  for(int i = 0; i < d.np_tot[IDIR] ; i++) {
    real r  = d.x[IDIR](i); // cell center
    real ri = r - HALF_F * d.dx[IDIR](i); // left cell boundary
    real ro = r + HALF_F * d.dx[IDIR](i); // left cell boundary
    real drx = 0.1 * d.dx[IDIR](i); // size of sub-grid

    real dVx = FOUR_F / THREE_F * M_PI * (pow(ro,3) - pow(ri,3)); // cell volume

    real px, rhox;
    real cx  = ZERO_F;
    real vx = ZERO_F;
    real dMx = ZERO_F;

    if (r < rsh && r > rpns) {
      // divide into subgrids for accurate interpolation
      for (int ii = 0; ii < 10; ii++){

        real rx  = ri + (ii + 0.05) * drx; // subgrid cell center
        real rxi = ri + (ii)        * drx; // subgrid cell left boundary
        real rxo = ri + (ii + 1)    * drx; // subgrid cell right boundary
        real dVxx = FOUR_F / THREE_F * M_PI * (pow(rxo,3) - pow(rxi,3)); // sub-grid cell volume

        real cxx   = LinearInterpolation(z_arr, c_arr, rx);
        real vxx   = - cxx * LinearInterpolation(z_arr, M_arr, rx);
        real rhoxx = FABS(v1 / (vxx * pow(rx, 2)));

        dMx  +=       rhoxx * dVxx; // mass in cell
        vx   += vxx * rhoxx * dVxx; // momentum in subgrid
        cx   += cxx * rhoxx * dVxx; // spped of sound "momentum" in subgrid
      }
      cx   /= dMx;
      vx   /= dMx;
      rhox = dMx / dVx;
      px   = rhox * pow(cx, 2) / gammaGlob;
    } else if (r > rsh) {
      ComputeBondi(r, rhox, px, vx);
    }
    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
      for(int j = 0; j < d.np_tot[JDIR] ; j++) {
        real theta = d.x[JDIR](j);
        real phi = d.x[KDIR](k);
        d.Vc(PRS,k,j,i) = px;
        d.Vc(RHO,k,j,i) = rhox;
        if(r >= 5 && r <= 7){
          d.Vc(VX1,k,j,i) = vx;// - (0.01/d.Vc(RHO,k,j,i))*sqrt(3/M_PI)*(5*cos(theta)*pow(sin(theta), 3/2)*sin(M_PI*(5 - r)/(5 - 7)))/(4*r*r*sqrt(sin(theta)*sin(theta))); //l=1
          d.Vc(VX2,k,j,i) = ZERO_F;// - (0.01/d.Vc(RHO,k,j,i))*(sqrt(3*M_PI*sin(theta)*sin(theta)*sin(theta))*cos(M_PI*(5 - r)/(5 - 7)))/(2*r*(5 - 7));
          //d.Vc(VX1,k,j,i) = ZERO_F - (0.015/d.Vc(RHO,k,j,i))*(sqrt((15/(2*M_PI)))*pow(sin(theta), 3/2)*sin((M_PI*(5 - r))/(5 - 7))*(7*cos(2*theta) + 3))/(8*r*r*sqrt(sin(theta)*sin(theta))); //l=2
          //d.Vc(VX2,k,j,i) = ZERO_F - (0.015/d.Vc(RHO,k,j,i))*(sqrt(15*M_PI*sin(theta)*sin(theta)*sin(theta)/2)*cos(theta)*cos(M_PI*(5 - r)/(5 - 7)))/(2*r*(5 - 7));
          //d.Vc(VX1,k,j,i) = vx - (0.015/d.Vc(RHO,k,j,i))*sqrt(5/M_PI)*3*(27 + 56*cos(2*theta) + 77*cos(4*theta))*pow(sin(theta), 3/2)*sin(M_PI*(5 - r)/(5 - 7))/(128*r*r*sqrt(pow(sin(theta), 2))); //l=4
          //d.Vc(VX2,k,j,i) = ZERO_F + (0.015/d.Vc(RHO,k,j,i))*sqrt(5*M_PI*pow(sin(theta), 3))*3*(3 - 7*pow(cos(theta), 2))*cos(theta)*cos(M_PI*(5 - r)/(5 - 7))/(8*(5 - 7)*r);
        }
        else{
          d.Vc(VX1,k,j,i) = vx;
          d.Vc(VX2,k,j,i) = ZERO_F;
        }
        d.Vc(VX3,k,j,i) = ZERO_F;
      }
    }
  }
  z_arr.clear(); M_arr.clear(); c_arr.clear();

  d.SyncToDevice(); // Send it all, if needed
}
//----------------------------------------------------------------------------------------

void MySourceTerm(Hydro *hydro, const real t, const real dtin) {
  IdefixArray4D<real> Vc = data->hydro->Vc;
  IdefixArray4D<real> Uc = data->hydro->Uc;
  IdefixArray1D<real> x1 = data->x[IDIR];

  idefix_for("MySourceTerm",0,data->np_tot[KDIR],0,data->np_tot[JDIR],0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int k, int j, int i) {

      real c2 = gammaGlob * Vc(PRS,k,j,i) / Vc(RHO,k,j,i);
      real v2 = pow(Vc(VX1,k,j,i), 2) + pow(Vc(VX2,k,j,i), 2) + pow(Vc(VX3,k,j,i), 2);

      if (v2 / c2 < 1.0 || x1(i) < 0.5) { 
        real sx = gm1i * log(FABS(Vc(PRS,k,j,i) / pow(Vc(RHO,k,j,i), gammaGlob)));
        Uc(ENG,k,j,i) += dtin * Vc(RHO,k,j,i) * 
          HeatingCooling (x1(i), Vc(PRS,k,j,i), Aheat) * 
          exp(- pow(sx / s_rpns, 2));
      } 
    });
}
//----------------------------------------------------------------------------------------

void UserdefBoundary(Hydro *hydro, int dir, BoundarySide side, real t) {
  IdefixArray4D<real> Vc = data->hydro->Vc;
  IdefixArray1D<real> x1 = data->x[IDIR];
  const int nxi    = data->np_int[IDIR];
  const int ighost = data->nghost[IDIR];

  if(dir==IDIR) {
    data->hydro->boundary->BoundaryFor("UserDefBoundary", dir, side,
      KOKKOS_LAMBDA (int k, int j, int i) {
        if (side == right) {

          ComputeBondi(x1(i), Vc(RHO,k,j,i), Vc(PRS,k,j,i), Vc(VX1,k,j,i));
          Vc(VX2,k,j,i) = ZERO_F;
          Vc(VX3,k,j,i) = ZERO_F;

        } else if (side == left) {

          // const int iref = ighost + side*(nxi-1); // outflow          
          // Vc(RHO,k,j,i) =   Vc(RHO,k,j,iref);
          // Vc(PRS,k,j,i) =   Vc(PRS,k,j,iref);
          // Vc(VX1,k,j,i) =   - Vc(VX1,k,j,iref); // v_rpns; // fixed velocity
          // Vc(VX2,k,j,i) =   Vc(VX2,k,j,iref);
          // Vc(VX3,k,j,i) =   Vc(VX3,k,j,iref);

          // reflection
          const int iref = 2*(ighost + side*(nxi-1)) - i - 1; // reflection
          Vc(RHO,k,j,i)  =  Vc(RHO,k,j,iref);
          Vc(PRS,k,j,i)  =  Vc(PRS,k,j,iref);
          Vc(VX1,k,j,i)  =  - Vc(VX1,k,j,iref); // reflection
          Vc(VX2,k,j,i)  =  Vc(VX2,k,j,iref);
          Vc(VX3,k,j,i)  =  Vc(VX3,k,j,iref);
        }
      });   
  }
}
//----------------------------------------------------------------------------------------

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {
}
//----------------------------------------------------------------------------------------

real ShockCompressionFactor (real M1x) {
  return (7.0 / 3.0) *  M1x * M1x / ((1.0 / 3.0) *  M1x * M1x + 2.0);
}
//----------------------------------------------------------------------------------------

real v_from_lM(real lMx, real c2x) {
  return - exp(lMx) * pow(FABS(c2x), 0.5);
}
//----------------------------------------------------------------------------------------

real FindEntropy(real zx, real vx, real c2x) {
  real rhox = FABS(v1 / vx) / (zx*zx);
  real px   = c2x * rhox / gammaGlob;
  return gm1i * log(px / pow(rhox, gammaGlob));
}
//----------------------------------------------------------------------------------------

real Phi (real lMx, real c2x, real zx, real Ah) {
  real vx   = v_from_lM(lMx, c2x);
  real Sx   = FindEntropy(zx, vx, c2x);
  real rhox = FABS(v1 / vx) / (zx*zx);
  real px   = rhox * c2x / gammaGlob; 
  return HeatingCooling(zx, px, Ah) * CoolingCutoff(Sx, s_rpns) / v_from_lM(lMx, c2x);
}
//----------------------------------------------------------------------------------------

real HeatingCooling (real zx, real px, real Ah) {
  return (Ah / (zx * zx) - Acool * pow(px, 1.5)); 
}
//----------------------------------------------------------------------------------------

real CoolingCutoff(real Sx, real Sx_rpns) {
  // the specific heating/cooling finction.
  if (fabs(Sx_rpns) < 1.0e-12) {
    return 1.0; // no cooling cutoff on the 1st calculation
  } else {
    return exp(- pow(Sx / Sx_rpns, 2));
  }
}
//----------------------------------------------------------------------------------------

real CalcStationaryShock (real M_rshx, real c2_rshx, real r_shx, real Ah) {
  int i, j;
  int n = 2;
  real t0, t1, dt; 
  real *u0, *u1;
  real lMmin = log(Mmin);   // log of minimum Mach number
  real lMmax = log(M_rshx); // log of post-shock Mach number

  t0 = lMmax; 
  dt = (lMmin - lMmax) / fabs(Nmax);

  u0 = new real[n];
  u0[0] = c2_rshx; // post-shock sound speed
  u0[1] = r_shx;   // shock position

  z_arr.push_back(r_shx);
  M_arr.push_back(M_rshx);
  c_arr.push_back(std::sqrt(c2_rshx));

  while (FABS(t0 - lMmin) > 1.0e-12) {
    t1 = t0 + dt;
    u1 = RK4Update (t0, n, u0, dt, RK4_RHS, Ah);
    t0 = t1;
    for ( i = 0; i < n; i++ ) {
      u0[i] = u1[i];
    }
    z_arr.push_back(u0[1]);
    M_arr.push_back(exp(t0));
    c_arr.push_back(std::sqrt(u0[0]));
    delete [] u1;
  }
  return z_arr[z_arr.size() - 1];
}

//----------------------------------------------------------------------------------------
//  Evaluates the right hand side of a vector ODE.
//  Parameters:
//    Input, t = the current time.
//    Input, n = the dimension of the system.
//    Input, U[m] = the current solution value.
//    Output, Real the value of the derivative, dU/dT.

real *RK4_RHS (real t, int n, real u[], real Ah) {
  real *u_dot;
  u_dot = new real[n];

  real Phix = Phi(t, u[0], u[1], Ah);

  real b1   = -u[0]; 
  real b2   = -exp(2*t) * u[0];
  real a11  = 0.5 + 3.0;
  real a12  = 2 * u[0] / u[1] - (4.0/3.0) * Phix;
  real a21  = 0.5 * exp(2*t) + 3.0; 
  real a22  = 1.0 / (u[1]*u[1]) - Phix;

  u_dot[0] = (a22*b1-a12*b2) / (a11*a22-a12*a21); // dc2_dlM;
  u_dot[1] = (a11*b2-a21*b1) / (a11*a22-a12*a21); //dr_dlM;
 
  return u_dot;
}

//----------------------------------------------------------------------------------------
//  Purpose: takes one Runge-Kutta step for a vector ODE
//      du/dt = f (t, u), u(t0) = u0
//  Parameters:
//    Input, Real T0, the current time.
//    Input, int M, the spatial dimension.
//    Input, Real U0[M], the solution estimate at the current time.
//    Input, Real DT, the time step.
//    Input, Real *F ( Real T, int M, Real U[] ), a function which evaluates
//    the derivative, or right hand side of the problem.
//    Output, Real RK4VEC[M], the estimate at time T0+DT.

real *RK4Update ( real t0, int m, real u0[], real dt, 
  real *f ( real t, int m, real u[], real Ah), real Ah) {
  real *f0;
  real *f1;
  real *f2;
  real *f3;
  real t1, t2, t3;
  real *u;
  real *u1;
  real *u2;
  real *u3;
  int i;

  f0 = f ( t0, m, u0, Ah);

  t1 = t0 + dt / 2.0;
  u1 = new real[m];
  for ( i = 0; i < m; i++ )
  {
    u1[i] = u0[i] + dt * f0[i] / 2.0;
  }
  f1 = f ( t1, m, u1, Ah);

  t2 = t0 + dt / 2.0;
  u2 = new real[m];
  for ( i = 0; i < m; i++ )
  {
    u2[i] = u0[i] + dt * f1[i] / 2.0;
  }
  f2 = f ( t2, m, u2, Ah);

  t3 = t0 + dt;
  u3 = new real[m];
  for ( i = 0; i < m; i++ )
  {
     u3[i] = u0[i] + dt * f2[i];
  }
  f3 = f ( t3, m, u3, Ah);

  u = new real[m];
  for ( i = 0; i < m; i++ )
  {
     u[i] = u0[i] + dt * ( f0[i] + 2.0 * f1[i] + 2.0 * f2[i] + f3[i] ) / 6.0;
  }
  delete [] f0;
  delete [] f1;
  delete [] f2;
  delete [] f3;
  delete [] u1;
  delete [] u2;
  delete [] u3;

  return u;
}
//----------------------------------------------------------------------------------------

real LinearInterpolation(std::vector<real>& x, std::vector<real>& y, real x_interp) {
    int i_int = 0;
    for (int i = 0; i < x.size() - 1; i++) {
        if (x_interp <= x[i] && x_interp >= x[i + 1]) {
            i_int = i;
            break;
        }
    }
    return y[i_int] + (y[i_int+1]-y[i_int]) * (x_interp-x[i_int]) / (x[i_int+1]-x[i_int]);;
}
//----------------------------------------------------------------------------------------

void CoarsenFunction(DataBlock &data) {
  IdefixArray2D<int> coarseningLevel = data.coarseningLevel[KDIR];
  IdefixArray1D<real> th = data.x[JDIR];
  idefix_for("set_coarsening", 0, data.np_tot[JDIR], 0, data.np_tot[IDIR],
      KOKKOS_LAMBDA(int j,int i) {
        coarseningLevel(j,i) = 1;
        if(th(j) < 0.3 || M_PI - th(j) < 0.3) {
          coarseningLevel(j,i) = 2;
          if(th(j)<0.1 || M_PI - th(j) < 0.1) {
            coarseningLevel(j,i) = 3;
          }
        }
      });
}
//----------------------------------------------------------------------------------------

real MachNumberResidual(real Mx, real rx) {
  return 2.0 * pow(rx, 4) * pow(Mx, 2) * pow((2.0 + 6.0 * rx) / (rx * (pow(Mx, 2) + 6.0)), 7) - 1.0;
}
//----------------------------------------------------------------------------------------

real MachNumberBisect(real r, real t_min, real t_max) {
  // Parameters
  const int max_iterations = 20;
  const real tol_residual = 1.0e-6;
  const real tol_temperature = 1.0e-6;

  // Find initial residuals
  real res_min = MachNumberResidual(t_min, r);
  real res_max = MachNumberResidual(t_max, r);
  if (std::abs(res_min) < tol_residual) {
    return t_min;
  }
  if (std::abs(res_max) < tol_residual) {
    return t_max;
  }
  if ((res_min < 0.0 && res_max < 0.0) || (res_min > 0.0 && res_max > 0.0)) {
    return NAN;
  }

  // Iterate to find root
  real t_mid;
  for (int i = 0; i < max_iterations; ++i) {
    t_mid = (t_min + t_max) / 2.0;
    if (t_max - t_min < tol_temperature) {
      return t_mid;
    }
    real res_mid = MachNumberResidual(t_mid, r);
    if (std::abs(res_mid) < tol_residual) {
      return t_mid;
    }
    if ((res_mid < 0.0 && res_min < 0.0) || (res_mid > 0.0 && res_min > 0.0)) {
      t_min = t_mid;
      res_min = res_mid;
    } else {
      t_max = t_mid;
      res_max = res_mid;
    }
  }
  return t_mid;
}
//----------------------------------------------------------------------------------------

void static ComputeBondi(real rx, real& rho, real& pgas, real& vx) {
  
  // works only for gamma=4/3
  // Change units:
  // v1           -> c_infty
  // R1 = GR/v1^2 -> R_B = GM/cinfty^2

  real r1_over_rb   = 0.006431; // radius where M=5 in units of R_B
  real cinf_over_v1 = 0.0801935;
  real vr_bondi     = -15.9887 * cinf_over_v1; // radial velocity at r1 in unit of v1;;

  real radx = r1_over_rb * rx;
  real rsonic = ONE_FOURTH_F; // sonic radius in units of R_B
  real Mx, cx;

  if (radx > rsonic) {
    Mx   = MachNumberBisect(radx, 0.001, 0.999999999999);
  } else {
    Mx   = MachNumberBisect(radx, 1.00000000001, 10.0);
  }

  // convert back to units v1 and R1
  cx   = cinf_over_v1 * std::sqrt((1.0 / radx + 3.0) / (0.5 * Mx*Mx + 3.0));
  vx  =  - Mx * cx;
  rho  = vr_bondi / (vx * pow(rx, 2)); 
  pgas = rho * pow(cx, 2) * 0.75;

  return;
}
//----------------------------------------------------------------------------------------
