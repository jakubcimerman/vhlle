#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <cfloat>

#include "fld.h"
#include "eos.h"
#include "ic3F.h"
#include "rmn.h"
#include "inc.h"
#include "nucleon.h"

using namespace std;

IC3F::IC3F(Fluid *f_p, Fluid *f_t, int _nevents, double _snn, double _b_min, double _b_max, int _projA, int _targA, int _projZ, int _targZ, double _Rg, double _tau0, EoS* eos) {
 nx = f_p->getNX();
 ny = f_p->getNY();
 nz = f_p->getNZ();
 dx = f_p->getDx();
 dy = f_p->getDy();
 dz = f_p->getDz();
 xmin = f_p->getX(0);
 xmax = f_p->getX(nx - 1);
 ymin = f_p->getY(0);
 ymax = f_p->getY(ny - 1);
 zmin = f_p->getZ(0);
 zmax = f_p->getZ(nz - 1);

 snn = _snn;
 b_min = _b_min;
 b_max = _b_max;
 nevents = _nevents;
 projA = _projA;
 targA = _targA;
 projZ = _projZ;
 targZ = _targZ;

 const double vcoll = sqrt(1-pow(2*nucleon_mass/snn, 2)); // velocity of the beam
 const double gamma = 1/sqrt(1 - vcoll * vcoll); // gamma factor of nuclei
 const double rap_beam = 0.5 * log((1 + vcoll) / (1-vcoll)); // beam rapidity
 const double Rproj = 1.1 * pow(projA, 0.333) - 0.656 * pow(projA, -0.333); // projectile nucleus radius
 const double Rtarg = 1.1 * pow(targA, 0.333) - 0.656 * pow(targA, -0.333); // target nucleus radius
 const double z0_proj = - (Rproj + 2 * WSdelta) / gamma;
 const double z0_targ = (Rtarg + 2 * WSdelta) / gamma;
 tau0 = _tau0;
 Rgx = _Rg;
 Rgy = _Rg;
 Rgz = _Rg;
 nsmoothx = (int)(3.0 * Rgx / dx);
 nsmoothy = (int)(3.0 * Rgy / dy);
 nsmoothz = (int)(3.0 * Rgz / dz);

 // allocate tensor field for projectile
 T00_p = new double**[nx];
 T0z_p = new double**[nx];
 QB_p = new double**[nx];
 QE_p = new double**[nx];
 for (int ix = 0; ix < nx; ix++) {
  T00_p[ix] = new double*[ny];
  T0z_p[ix] = new double*[ny];
  QB_p[ix] = new double*[ny];
  QE_p[ix] = new double*[ny];
  for (int iy = 0; iy < ny; iy++) {
   T00_p[ix][iy] = new double[nz];
   T0z_p[ix][iy] = new double[nz];
   QB_p[ix][iy] = new double[nz];
   QE_p[ix][iy] = new double[nz];
   for (int iz = 0; iz < nz; iz++) {
    T00_p[ix][iy][iz] = 0.0;
    T0z_p[ix][iy][iz] = 0.0;
    QB_p[ix][iy][iz] = 0.0;
    QE_p[ix][iy][iz] = 0.0;
   }
  }
 }

  // allocate tensor field for target
 T00_t = new double**[nx];
 T0z_t = new double**[nx];
 QB_t = new double**[nx];
 QE_t = new double**[nx];
 for (int ix = 0; ix < nx; ix++) {
  T00_t[ix] = new double*[ny];
  T0z_t[ix] = new double*[ny];
  QB_t[ix] = new double*[ny];
  QE_t[ix] = new double*[ny];
  for (int iy = 0; iy < ny; iy++) {
   T00_t[ix][iy] = new double[nz];
   T0z_t[ix][iy] = new double[nz];
   QB_t[ix][iy] = new double[nz];
   QE_t[ix][iy] = new double[nz];
   for (int iz = 0; iz < nz; iz++) {
    T00_t[ix][iy][iz] = 0.0;
    T0z_t[ix][iy][iz] = 0.0;
    QB_t[ix][iy][iz] = 0.0;
    QE_t[ix][iy][iz] = 0.0;
   }
  }
 }
 nucleons.clear();

 // random generators
 random_device rd;
 mt19937 gen(rd());
 uniform_real_distribution<> dis_proj(- Rproj - 4.0, Rproj + 4.0);
 uniform_real_distribution<> dis_targ(- Rtarg - 4.0, Rtarg + 4.0);
 uniform_real_distribution<> dis_uniform(0.0, 1.0);

 if (nevents < 0) {
  cout << "nevents has to be positive" << endl;
  exit(1);
 } else if (nevents == 0) {
  double ep_norm = 0, nbp_norm = 0, nqp_norm = 0;
  double et_norm = 0, nbt_norm = 0, nqt_norm = 0;
  for (int ix = 0; ix < f_p->getNX(); ix++)
   for (int iy = 0; iy < f_p->getNY(); iy++)
    for (int iz = 0; iz < f_p->getNZ(); iz++) {
     if (b_max != b_min)
      cout << "symmetric initial state does not generate random impact parameter, value b_min is used" << endl;
     double b = b_min;
     double x = f_p->getX(ix);
     double y = f_p->getY(iy);
     double eta = f_p->getZ(iz);
     double z_p = tau0 * (sinh(eta - rap_beam) + sinh(rap_beam)) / (cosh(rap_beam));
     double z_t = tau0 * (sinh(eta + rap_beam) - sinh(rap_beam)) / (cosh(rap_beam));

     double r_p = sqrt((x - b/2) * (x - b/2) + y * y + gamma * gamma * (z_p - z0_proj) * (z_p - z0_proj));
     double r_t = sqrt((x + b/2) * (x + b/2) + y * y + gamma * gamma * (z_t - z0_targ) * (z_t - z0_targ));

     // calculate energy-momentum tensor and charges in Cartesian frame at ix,iy,iz
     double ep = 0.17 * nucleon_mass / (1 + exp((r_p - Rproj) / WSdelta));
     double et = 0.17 * nucleon_mass / (1 + exp((r_t - Rtarg) / WSdelta));
     double nbp = ep / nucleon_mass;
     double nbt = et / nucleon_mass;
     double nqp = nbp * projZ / projA;
     double nqt = nbt * targZ / targA;
     double pp = eos->p(ep, nbp, 0, nqp);
     double pt = eos->p(et, nbt, 0, nqt);

     double Ttt_p = (ep + pp) * gamma * gamma - pp;
     double Ttz_p = (ep + pp) * gamma * gamma * vcoll;
     double Tzz_p = (ep + pp) * gamma * gamma * vcoll * vcoll + pp;
     double Ttt_t = (et + pt) * gamma * gamma - pt;
     double Ttz_t = (et + pt) * gamma * gamma * (-vcoll);
     double Tzz_t = (et + pt) * gamma * gamma * vcoll * vcoll + pt;
     double Nbt_p = nbp * gamma;
     double Nbz_p = nbp * gamma * vcoll;
     double Nbt_t = nbt * gamma;
     double Nbz_t = nbt * gamma * (-vcoll);
     double Nqt_p = nqp * gamma;
     double Nqz_p = nqp * gamma * vcoll;
     double Nqt_t = nqt * gamma;
     double Nqz_t = nqt * gamma * (-vcoll);

     // transform energy-momentum tensor and charges to hyperbolic frame
     double z = tau0 * sinh(eta);
     double t = tau0 * cosh(eta);
     T00_p[ix][iy][iz] = (t*t*Ttt_p - 2*z*t*Ttz_p + z*z*Tzz_p)/pow(tau0,2);
     T00_t[ix][iy][iz] = (t*t*Ttt_t - 2*z*t*Ttz_t + z*z*Tzz_t)/pow(tau0,2);
     T0z_p[ix][iy][iz] = (-z*t*Ttt_p + (t*t+z*z)*Ttz_p - z*t*Tzz_p)/pow(tau0,2);
     T0z_t[ix][iy][iz] = (-z*t*Ttt_t + (t*t+z*z)*Ttz_t - z*t*Tzz_t)/pow(tau0,2);
     QB_p[ix][iy][iz] = (t*Nbt_p - z*Nbz_p) / tau0;
     QB_t[ix][iy][iz] = (t*Nbt_t - z*Nbz_t) / tau0;
     QE_p[ix][iy][iz] = (t*Nqt_p - z*Nqz_p) / tau0;
     QE_t[ix][iy][iz] = (t*Nqt_t - z*Nqz_t) / tau0;

     // calculate normalization factors
     double nsp, vxp, vyp, vzp, mubp, muqp, musp, tp;
     double nst, vxt, vyt, vzt, mubt, muqt, must, tt;
     double Qp[7] = {T00_p[ix][iy][iz], 0, 0, T0z_p[ix][iy][iz], QB_p[ix][iy][iz], QE_p[ix][iy][iz], 0};
     double Qt[7] = {T00_t[ix][iy][iz], 0, 0, T0z_t[ix][iy][iz], QB_t[ix][iy][iz], QE_t[ix][iy][iz], 0};
     transformPV(eos, Qp, ep, pp, nbp, nqp, nsp, vxp, vyp, vzp);
     transformPV(eos, Qt, et, pt, nbt, nqt, nst, vxt, vyt, vzt);
     vzp = eta + 1. / 2. * log((1. + vzp) / (1. - vzp));
     vxp = vxp * cosh(vzp - eta) / cosh(vzp);
     vyp = vyp * cosh(vzp - eta) / cosh(vzp);
     vzt = eta + 1. / 2. * log((1. + vzt) / (1. - vzt));
     vxt = vxt * cosh(vzt - eta) / cosh(vzt);
     vyt = vyt * cosh(vzt - eta) / cosh(vzt);
     eos->eos(ep, nbp, nqp, nsp, tp, mubp, muqp, musp, pp);
     eos->eos(et, nbt, nqt, nst, tt, mubt, muqt, must, pt);
     const double cosh_int = (sinh(eta + 0.5 * dz) - sinh(eta - 0.5 * dz)) / dz;
     const double sinh_int = (cosh(eta + 0.5 * dz) - cosh(eta - 0.5 * dz)) / dz;
     ep_norm += tau0 * (ep + pp) / (1. - vxp * vxp - vyp * vyp - tanh(vzp) * tanh(vzp)) *
             (cosh_int - tanh(vzp) * sinh_int) -
         tau0 * pp * cosh_int;
     et_norm += tau0 * (et + pt) / (1. - vxt * vxt - vyt * vyt - tanh(vzt) * tanh(vzt)) *
             (cosh_int - tanh(vzt) * sinh_int) -
         tau0 * pt * cosh_int;
     nbp_norm += QB_p[ix][iy][iz];
     nbt_norm += QB_t[ix][iy][iz];
     nqp_norm += QE_p[ix][iy][iz];
     nqt_norm += QE_t[ix][iy][iz];
  }
  ep_norm *= dx * dy * dz;
  et_norm *= dx * dy * dz;
  nbp_norm *= tau0 * dx * dy * dz;
  nbt_norm *= tau0 * dx * dy * dz;
  nqp_norm *= tau0 * dx * dy * dz;
  nqt_norm *= tau0 * dx * dy * dz;

  // normalize to E=sqrt(sNN)*A/2, Nb=A, Nq=Z
  cout << ep_norm << " " << nbp_norm << " " << nqp_norm << endl;
  cout << et_norm << " " << nbt_norm << " " << nqt_norm << endl;
  double norm_p = 0;
  double norm_t = 0;
  for (int ix = 0; ix < f_p->getNX(); ix++)
   for (int iy = 0; iy < f_p->getNY(); iy++)
    for (int iz = 0; iz < f_p->getNZ(); iz++) {
     QB_p[ix][iy][iz] *= projA / nbp_norm;
     QB_t[ix][iy][iz] *= targA / nbt_norm;
     QE_p[ix][iy][iz] *= projZ / nqp_norm;
     QE_t[ix][iy][iz] *= targZ / nqt_norm;
     T00_p[ix][iy][iz] *= snn * projA / (2 * ep_norm);
     T00_t[ix][iy][iz] *= snn * targA / (2 * et_norm);
     T0z_p[ix][iy][iz] *= snn * projA / (2 * ep_norm);
     T0z_t[ix][iy][iz] *= snn * targA / (2 * et_norm);
  }
 } else {
  for (int iev = 0; iev < nevents; iev++) {
   // generating impact parameter
   bool generated = false;
   double bx, by, b;
   if (b_max > Rproj + Rtarg) b_max = Rproj + Rtarg;
   if (b_min < 0) b_min = 0;
   if (b_min > b_max) {
    cout << "b_min has to be smaller than b_max" << endl;
    exit(1);
   }
   if (b_min == b_max) {
    b = b_min;
    generated = true;
   }
   uniform_real_distribution<> dis_impact(- b_max, b_max);
   while (!generated) {
    bx = dis_impact(gen);
    by = dis_impact(gen);
    b = sqrt(bx*bx + by*by);
    if (b < b_max && b > b_min) generated = true;
   }

   cout << "generated impact parameter: " << b << endl;

   // output for debugging
   //ofstream fout("nucleons.dat");
   std::vector<Nucleon> nucl1;
   nucleons.push_back(nucl1);
   nucleons[iev].clear();

   // generating nucleons of projectile nucleus
   for (int i = 0; i < projA; i++) {
    generated = false;
    double x, y, z, r;
    while (!generated) {
     x = dis_proj(gen);
     y = dis_proj(gen);
     z = dis_proj(gen);
     r = sqrt(x*x + y*y + z*z);
     if (dis_uniform(gen) < 1/(1 + exp((r - Rproj) / WSdelta))) generated = true;
    }
    x += b / 2;
    z = z / gamma + z0_proj;
    double eta = asinh(z * cosh(rap_beam) / tau0 - sinh(rap_beam)) + rap_beam;
    int charge = i < projZ ? 1 : 0;
    nucleons[iev].push_back(Nucleon(x, y, eta, rap_beam, charge));
    makeSmoothPart(x, y, eta, charge, rap_beam, true);
    //fout << x << " " << y << " " << z << " " << r << " " << eta << endl;
   }

   // generating nucleons of target nucleus
   for (int i = 0; i < targA; i++) {
    generated = false;
    double x, y, z, r;
    while (!generated) {
     x = dis_targ(gen);
     y = dis_targ(gen);
     z = dis_targ(gen);
     r = sqrt(x*x + y*y + z*z);
     if (dis_uniform(gen) < 1/(1 + exp((r - Rtarg) / WSdelta))) generated = true;
    }
    x -= b / 2;
    z = z / gamma + z0_targ;
    double eta = asinh(z * cosh(rap_beam) / tau0 + sinh(rap_beam)) - rap_beam;
    int charge = i < targZ ? 1 : 0;
    nucleons[iev].push_back(Nucleon(x, y, eta, -rap_beam, charge));
    makeSmoothPart(x, y, eta, charge, -rap_beam, false);
    //fout << x << " " << y << " " << z << " " << r << " " << eta << endl;
   }
  }
 }
 //fout.close();
}

IC3F::~IC3F() {
 for (int ix = 0; ix < nx; ix++) {
  for (int iy = 0; iy < ny; iy++) {
   delete[] T00_p[ix][iy];
   delete[] T0z_p[ix][iy];
   delete[] QB_p[ix][iy];
   delete[] QE_p[ix][iy];
   delete[] T00_t[ix][iy];
   delete[] T0z_t[ix][iy];
   delete[] QB_t[ix][iy];
   delete[] QE_t[ix][iy];
  }
  delete[] T00_p[ix];
  delete[] T0z_p[ix];
  delete[] QB_p[ix];
  delete[] QE_p[ix];
  delete[] T00_t[ix];
  delete[] T0z_t[ix];
  delete[] QB_t[ix];
  delete[] QE_t[ix];
 }
 delete[] T00_p;
 delete[] T0z_p;
 delete[] QB_p;
 delete[] QE_p;
 delete[] T00_t;
 delete[] T0z_t;
 delete[] QB_t;
 delete[] QE_t;
}

std::vector<std::vector<Nucleon>> IC3F::getNucleons(void) {
 return nucleons;
}

void IC3F::makeSmoothPart(double x, double y, double eta, int Charge, double rap, bool isProjectile) {
 int ixc = (int)round((x - xmin) / dx);
 int iyc = (int)round((y - ymin) / dy);
 int izc = (int)round((eta - zmin) / dz);
 double m = nucleon_mass;
 double norm_gauss = 0.0;
 for (int ix = ixc - nsmoothx; ix < ixc + nsmoothx + 1; ix++)
  for (int iy = iyc - nsmoothy; iy < iyc + nsmoothy + 1; iy++)
   for (int iz = izc - nsmoothz; iz < izc + nsmoothz + 1; iz++)
    if (ix > 0 && ix < nx && iy > 0 && iy < ny && iz > 0 && iz < nz) {
     const double xdiff = x - (xmin + ix * dx);
     const double ydiff = y - (ymin + iy * dy);
     const double zdiff = eta - (zmin + iz * dz);
     norm_gauss +=
         exp(-xdiff * xdiff / 2 / Rgx / Rgx - ydiff * ydiff / 2 / Rgy / Rgy -
             zdiff * zdiff / 2 / Rgz / Rgz * tau0 * tau0 * cosh(eta) * cosh(eta) *
                 cosh(rap) * cosh(rap));
    }

 for (int ix = ixc - nsmoothx; ix < ixc + nsmoothx + 1; ix++)
  for (int iy = iyc - nsmoothy; iy < iyc + nsmoothy + 1; iy++)
   for (int iz = izc - nsmoothz; iz < izc + nsmoothz + 1; iz++)
    if (ix > 0 && ix < nx && iy > 0 && iy < ny && iz > 0 && iz < nz) {
     const double xdiff = x - (xmin + ix * dx);
     const double ydiff = y - (ymin + iy * dy);
     const double zdiff = eta - (zmin + iz * dz);
     double weight;
      weight = 1.0 / norm_gauss *
          exp(-xdiff * xdiff / 2 / Rgx / Rgx - ydiff * ydiff / 2 / Rgy / Rgy -
             zdiff * zdiff / 2 / Rgz / Rgz * tau0 * tau0 * cosh(eta) * cosh(eta) *
                 cosh(rap) * cosh(rap));
     if (weight != weight || fabs(weight) > DBL_MAX) {
      weight = 0.0;
     }
     if (isProjectile) {
      T00_p[ix][iy][iz] += weight * m * cosh(rap - eta + zdiff);
      T0z_p[ix][iy][iz] += weight * m * sinh(rap - eta + zdiff);
      QB_p[ix][iy][iz] += weight;
      QE_p[ix][iy][iz] += Charge * weight;
     } else {
      T00_t[ix][iy][iz] += weight * m * cosh(rap - eta + zdiff);
      T0z_t[ix][iy][iz] += weight * m * sinh(rap - eta + zdiff);
      QB_t[ix][iy][iz] += weight;
      QE_t[ix][iy][iz] += Charge * weight;
     }
     //}
    }
}

void IC3F::setIC(Fluid* f_p, Fluid* f_t, EoS* eos) {
 double E = 0.0, Px = 0.0, Py = 0.0, Pz = 0.0, Nb = 0.0, S = 0.0;
 double Q_p[7], e_p, p_p, nb_p, nq_p, ns_p, vx_p, vy_p, vz_p;
 double Q_t[7], e_t, p_t, nb_t, nq_t, ns_t, vx_t, vy_t, vz_t;
 // output for debugging
 //ofstream feout("energy.dat");
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    if (nevents == 0) {
     Q_p[T_] = T00_p[ix][iy][iz];  // /tau for Milne
     Q_p[X_] = 0.0;
     Q_p[Y_] = 0.0;
     Q_p[Z_] = T0z_p[ix][iy][iz];
     Q_p[NB_] = QB_p[ix][iy][iz];
     Q_p[NQ_] = QE_p[ix][iy][iz];
     Q_p[NS_] = 0.0;

     Q_t[T_] = T00_t[ix][iy][iz];  // /tau for Milne
     Q_t[X_] = 0.0;
     Q_t[Y_] = 0.0;
     Q_t[Z_] = T0z_t[ix][iy][iz];
     Q_t[NB_] = QB_t[ix][iy][iz];
     Q_t[NQ_] = QE_t[ix][iy][iz];
     Q_t[NS_] = 0.0;
    } else {
     Q_p[T_] = T00_p[ix][iy][iz] / nevents / dx / dy / dz / tau0;  // /tau for Milne
     Q_p[X_] = 0.0;
     Q_p[Y_] = 0.0;
     Q_p[Z_] = T0z_p[ix][iy][iz] / nevents / dx / dy / dz / tau0;
     Q_p[NB_] = QB_p[ix][iy][iz] / nevents / dx / dy / dz / tau0;
     Q_p[NQ_] = QE_p[ix][iy][iz] / nevents / dx / dy / dz / tau0;
     Q_p[NS_] = 0.0;

     Q_t[T_] = T00_t[ix][iy][iz] / nevents / dx / dy / dz / tau0;  // /tau for Milne
     Q_t[X_] = 0.0;
     Q_t[Y_] = 0.0;
     Q_t[Z_] = T0z_t[ix][iy][iz] / nevents / dx / dy / dz / tau0;
     Q_t[NB_] = QB_t[ix][iy][iz] / nevents / dx / dy / dz / tau0;
     Q_t[NQ_] = QE_t[ix][iy][iz] / nevents / dx / dy / dz / tau0;
     Q_t[NS_] = 0.0;
    }

    transformPV(eos, Q_p, e_p, p_p, nb_p, nq_p, ns_p, vx_p, vy_p, vz_p);
    transformPV(eos, Q_t, e_t, p_t, nb_t, nq_t, ns_t, vx_t, vy_t, vz_t);
    if (e_p < 1e-7) {
     e_p = 1e-10; nb_p = nq_p = 0.0;
     vx_p = vy_p = vz_p = 0.0;
    }
    if (e_t < 1e-7) {
     e_t = 1e-10; nb_t = nq_t = 0.0;
     vx_t = vy_t = vz_t = 0.0;
    }
    Cell* c_p = f_p->getCell(ix, iy, iz);
    c_p->setPrimVar(eos, tau0, e_p, nb_p, nq_p, ns_p, vx_p, vy_p, vz_p);

    Cell* c_t = f_t->getCell(ix, iy, iz);
    c_t->setPrimVar(eos, tau0, e_t, nb_t, nq_t, ns_t, vx_t, vy_t, vz_t);

    //feout << f_p->getX(ix) << " " << f_p->getY(iy) << " " << f_p->getZ(iz) << " " << e_p+e_t << endl;

    if (e_p > 0.) c_p->setAllM(1.);
    if (e_t > 0.) c_t->setAllM(1.);
    /*const double gamma = 1.0 / sqrt(1.0 - vx * vx - vy * vy - vz * vz);
    double u[4] = {gamma, gamma * vx, gamma * vy, gamma * vz};
    double eta = zmin + iz * dz;
    E += tau0 * ((e + p) * u[0] * (u[0] * cosh(eta) + u[3] * sinh(eta)) -
                 p * cosh(eta)) *
         dx * dy * dz;
    // if(zmin+iz*dz>0.)
    Pz += tau0 * ((e + p) * u[0] * (u[0] * sinh(eta) + u[3] * cosh(eta)) -
                  p * sinh(eta)) *
          dx * dy * dz;
    Px += tau0 * (e + p) * u[1] * u[0] * dx * dy * dz;
    Py += tau0 * (e + p) * u[2] * u[0] * dx * dy * dz;
    Nb += tau0 * nb * u[0] * dx * dy * dz;
    S += tau0 * eos->s(e, nb, nq, ns) * u[0] * dx * dy * dz;

    if (iz == (int)nz/2 and iy == (int)ny/2) {
      cout << "x " << xmin + ix * dx << "  " << vx << "  " << e << endl;
    }
    if (iz == (int)nz/2 and ix == (int)ny/2) {
      cout << "y " << ymin + iy * dy << "  " << vy << "  " << e << endl;
    }
   }
 cout << "hydrodynamic E = " << E << "  Pz = " << Pz << "  Nbar = " << Nb
      << endl
      << "  Px = " << Px << "  Py = " << Py << endl;
 cout << "initial_entropy S_ini = " << S << endl;*/
   }
 //feout.close();
 //exit(1);
}
