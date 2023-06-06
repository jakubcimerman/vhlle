#include <cmath>
#include <fstream>
#include <iostream>
#include <TMatrixDEigen.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>
#include <iomanip>

#include "multiHydro.h"
#include "hdo.h"
#include "fld.h"
#include "eos.h"
#include "rmn.h"
#include "trancoeff.h"
#include "cll.h"
#include "xsect.h"
#include "cornelius.h"
#include "nucleon.h"

#define OUTPI

using namespace std;

MultiHydro::MultiHydro(Fluid *_f_p, Fluid *_f_t, Fluid *_f_f, Hydro *_h_p,
 Hydro *_h_t, Hydro *_h_f, EoS *_eos, TransportCoeff *_trcoeff, double _dtau,
 double eCrit, double _sNN, double _xi_fa, double _lambda, double _formationTime,
 int _frictionModel, int _decreasingFormTime, double _xi_q, double _xi_h, std::vector<std::vector<Nucleon>> nucl)
{
 f_p = _f_p;
 f_t = _f_t;
 f_f = _f_f;
 h_p = _h_p;
 h_t = _h_t;
 h_f = _h_f;
 eos = _eos;
 trcoeff = _trcoeff;
 xsect = new CrossSections;
 nx = f_p->getNX();
 ny = f_p->getNY();
 nz = f_p->getNZ();
 dx = f_p->getDx();
 dy = f_p->getDy();
 dz = f_p->getDz();
 dtau = _dtau;
 tau0 = h_p->getTau();
 sNN = _sNN;
 xi_fa = _xi_fa;
 lambda = _lambda;
 formationTime = _formationTime;
 frictionModel = _frictionModel;
 decreasingFormTime = _decreasingFormTime;
 xi_q = _xi_q;
 xi_h = _xi_h;
 dtauf = formationTime / 10.0;
 nucleons = nucl;

 //---- Cornelius init
 double arrayDx[4] = {h_p->getDtau(), f_p->getDx(), f_p->getDy(), f_p->getDz()};
 cornelius = new Cornelius;
 cornelius->init(4, eCrit, arrayDx);
 ecrit = eCrit;
 vEff = 0.;
 vEff_p = 0.;
 vEff_t = 0.;
 vEff_f = 0.;

 // allocate field for oveall energy density
 MHeps = new double**[nx];
 MHepsPrev = new double**[nx];
 for (int ix = 0; ix < nx; ix++) {
  MHeps[ix] = new double*[ny];
  MHepsPrev[ix] = new double*[ny];
  for (int iy = 0; iy < ny; iy++) {
   MHeps[ix][iy] = new double[nz];
   MHepsPrev[ix][iy] = new double[nz];
   for (int iz = 0; iz < nz; iz++) {
    MHeps[ix][iy][iz] = 0.0;
    MHepsPrev[ix][iy][iz] = 0.0;
   }
  }
 }

 /* Debug of scattering rates calculation
 double u[4] = {1, 0, 0, 0};
 for (double p = 0.1; p <= 2; p+=0.1) {
  cout << p << " " << totalScatRate(p, 0.12, 0.3, u) << endl;
 }
 exit(1);*/
}

MultiHydro::~MultiHydro() {
 for (int ix = 0; ix < nx; ix++) {
  for (int iy = 0; iy < ny; iy++) {
   delete[] MHeps[ix][iy];
   delete[] MHepsPrev[ix][iy];
  }
  delete[] MHeps[ix];
  delete[] MHepsPrev[ix];
 }
 delete[] MHeps;
 delete[] MHepsPrev;
}

void MultiHydro::setFluids(Fluid *_f_p, Fluid *_f_t, Fluid *_f_f, Hydro *_h_p,
 Hydro *_h_t, Hydro *_h_f) {
 f_p = _f_p;
 f_t = _f_t;
 f_f = _f_f;
 h_p = _h_p;
 h_t = _h_t;
 h_f = _h_f;
 nx = f_p->getNX();
 ny = f_p->getNY();
 nz = f_p->getNZ();
 dx = f_p->getDx();
 dy = f_p->getDy();
 dz = f_p->getDz();
 dtau = h_p->getDtau();

 //---- Cornelius init
 double arrayDx[4] = {h_p->getDtau(), f_p->getDx(), f_p->getDy(), f_p->getDz()};
 cornelius = new Cornelius;
 cornelius->init(4, ecrit, arrayDx);

 resizeMHeps();
}

void MultiHydro::initOutput(const char *dir) {
 string outfreeze_p = dir, outfreeze_f = dir, outfreeze_t = dir, outfreeze_all = dir;
 outfreeze_p.append("/freezeout_p.dat");
 outfreeze_t.append("/freezeout_t.dat");
 outfreeze_f.append("/freezeout_f.dat");
 outfreeze_all.append("/freezeout_all.dat");
 fmhfreeze_p.open(outfreeze_p.c_str());
 fmhfreeze_t.open(outfreeze_t.c_str());
 fmhfreeze_f.open(outfreeze_f.c_str());
 fmhfreeze_all.open(outfreeze_all.c_str());
}

void MultiHydro::performStep()
{
 h_p->performStep();
 h_t->performStep();
 h_f->performStep();
 frictionSubstep();
}

void MultiHydro::frictionSubstep()
{
 // here it is assumed that projectile and target grids
 // have same dimensions and physical sizes
 for (int iy = 0; iy < f_p->getNY(); iy++)
  for (int iz = 0; iz < f_p->getNZ(); iz++)
   for (int ix = 0; ix < f_p->getNX(); ix++) {
    Cell *c_p = f_p->getCell(ix, iy, iz);
    Cell *c_t = f_t->getCell(ix, iy, iz);
    Cell *c_f = f_f->getCell(ix, iy, iz);
    double ep, pp, nbp, nqp, nsp, vxp, vyp, vzp;
    double et, pt, nbt, nqt, nst, vxt, vyt, vzt;
    double ef, pf, nbf, nqf, nsf, vxf, vyf, vzf;
    c_p->getPrimVar(eos, h_p->getTau(), ep, pp, nbp, nqp, nsp, vxp, vyp, vzp);
    c_t->getPrimVar(eos, h_t->getTau(), et, pt, nbt, nqt, nst, vxt, vyt, vzt);
    c_f->getPrimVar(eos, h_f->getTau(), ef, pf, nbf, nqf, nsf, vxf, vyf, vzf);
    double TCp, mubCp, muqCp, musCp, pCp;
    double TCt, mubCt, muqCt, musCt, pCt;
    eos->eos(ep, nbp, nqp, nsp, TCp, mubCp, muqCp, musCp, pCp);
    eos->eos(et, nbt, nqt, nst, TCt, mubCt, muqCt, musCt, pCt);
    // 4-velocities, u_p and u_t
    double gammap = 1.0/sqrt(1.0-vxp*vxp-vyp*vyp-vzp*vzp);
    double up [4] = {gammap,gammap*vxp,gammap*vyp,gammap*vzp};
    double gammat = 1.0/sqrt(1.0-vxt*vxt-vyt*vyt-vzt*vzt);
    double ut [4] = {gammat,gammat*vxt,gammat*vyt,gammat*vzt};
    double gammaf = 1.0/sqrt(1.0-vxf*vxf-vyf*vyf-vzf*vzf);
    double uf [4] = {gammaf,gammaf*vxf,gammaf*vyf,gammaf*vzf};
    TLorentzVector upLV(up[1], up[2], up[3], up[0]);
    TLorentzVector utLV(ut[1], ut[2], ut[3], ut[0]);
    TLorentzVector ufLV(uf[1], uf[2], uf[3], uf[0]);
    double upt_abs = sqrt(pow(up[0]+ut[0],2)-pow(up[1]+ut[1],2)-pow(up[2]+ut[2],2)-pow(up[3]+ut[3],2));
    double U_F[4] = {(up[0]+ut[0])/upt_abs, (up[1]+ut[1])/upt_abs, (up[2]+ut[2])/upt_abs, (up[3]+ut[3])/upt_abs};
    double flux_p [4] = {0.}, flux_t [4] = {0.}, flux_f [4] = {0.},
           flux_pf [4] = {0.}, flux_tf [4] = {0.};
     // 1. projectile-target friction
    if (ep>0. && et>0.) {
    // u_p^\mu u_t_\mu
    double uput = gammap*gammat*(1.0 - vxp*vxt - vyp*vyt - vzp*vzt);
    // s (Mandelstam) variable
    double s = 2.0*mN*mN*(1.0 + uput);
    double Ekin = s/(2.0*mN) - 2.0*mN;
    double sigmaT, sigmaE, sigmaP;
    if (Ekin <= 0.) {
     Ekin = 0;
     sigmaT = 0;
     sigmaE = 0;
     sigmaP = 0;
    } else {
     xsect->NN(Ekin, sigmaT, sigmaE, sigmaP);
    }
    // Moeller factor
    double Vrel = sqrt(uput*uput - 1.0);
    // friction coefficient
    double D_P = mN*Vrel*sigmaP;
    double D_E = mN*Vrel*sigmaE; // SAME cross section moment for testing
    double maxVal = max(max(Fermi(nbp)/mN, 3*TCp/(2*mN)), max(Fermi(nbt)/mN, 3*TCt/(2*mN)));
    double deltaV = sqrt(1-pow(maxVal+1,-2));
    double theta = 1 - exp(-pow(Vrel/deltaV,4));
    // qgb densities
    double zeta3 = 1.20205690315959;
    double nquark_p = 18 * zeta3 * pow(TCp, 3) / pow(M_PI, 2) + 2*pow(mubCp/3, 3);
    double ngluon_p = 16 * zeta3 * pow(TCp, 3) / pow(M_PI, 2);
    double nquark_t = 18 * zeta3 * pow(TCt, 3) / pow(M_PI, 2) + 2*pow(mubCt/3, 3);
    double ngluon_t = 16 * zeta3 * pow(TCt, 3) / pow(M_PI, 2);
    double dens_p, dens_t;
    //if (ep < 0.7) {
     dens_p = xi_h*pow(2*mN/sqrt(s), 0.5)*nbp;
    /*} else {
     dens_p = xi_q*pow(2*mN/sqrt(s), 0.5)*(nquark_p+ngluon_p)/3;
    }*/
    //if (et < 0.7) {
     dens_t = xi_h*pow(2*mN/sqrt(s), 0.5)*nbt;
    /*} else {
     dens_t = xi_q*pow(2*mN/sqrt(s), 0.5)*(nquark_t+ngluon_t)/3;
    }*/
    //double sigma_gg = 0.3; // fm^2
    //double D_QGP = mN*Vrel*sigma_gg;
    upLV.Boost(-vxt, -vyt, -vzt);
    utLV.Boost(-vxp, -vyp, -vzp);
    for(int i=0; i<4; i++){
     if (frictionModel == 1) {
      // Ivanov's friction terms
      flux_p[i] += -theta*dens_p*dens_t*(D_P*(up[i] - ut[i]) + D_E*(up[i] + ut[i]))*h_p->getDtau();
      flux_t[i] += -theta*dens_p*dens_t*(D_P*(ut[i] - up[i]) + D_E*(up[i] + ut[i]))*h_p->getDtau();
     } else {
      // simplified and parametrized friction terms
      double vpAbs = sqrt(pow(upLV[0],2) + pow(upLV[1],2) + pow(upLV[2],2))/upLV[3];
      double vtAbs = sqrt(pow(utLV[0],2) + pow(utLV[1],2) + pow(utLV[2],2))/utLV[3];
      double vLimit = sqrt(1-pow(2*0.938/sNN,2));
      if (vpAbs > vLimit) {
        upLV[3] = 1.0/sqrt(1.0-vLimit*vLimit);
        upLV[0] = upLV[0]*(vLimit/vpAbs)*sqrt((1.0-vpAbs*vpAbs)/(1.0-vLimit*vLimit));
        upLV[1] = upLV[1]*(vLimit/vpAbs)*sqrt((1.0-vpAbs*vpAbs)/(1.0-vLimit*vLimit));
        upLV[2] = upLV[2]*(vLimit/vpAbs)*sqrt((1.0-vpAbs*vpAbs)/(1.0-vLimit*vLimit));
      }
      if (vtAbs > vLimit) {
        utLV[3] = 1.0/sqrt(1.0-vLimit*vLimit);
        utLV[0] = utLV[0]*(vLimit/vtAbs)*sqrt((1.0-vtAbs*vtAbs)/(1.0-vLimit*vLimit));
        utLV[1] = utLV[1]*(vLimit/vtAbs)*sqrt((1.0-vtAbs*vtAbs)/(1.0-vLimit*vLimit));
        utLV[2] = utLV[2]*(vLimit/vtAbs)*sqrt((1.0-vtAbs*vtAbs)/(1.0-vLimit*vLimit));
      }
      flux_p[i] += -xi_h*upLV[(i+3)%4]*sqrt(ep*et)*h_p->getDtau()/lambda;
      flux_t[i] += -xi_h*utLV[(i+3)%4]*sqrt(ep*et)*h_p->getDtau()/lambda;
     }
     if (formationTime > 0) {
      addRetardedFriction((-flux_p[i]-flux_t[i])*f_p->getDx()*f_p->getDy()*f_p->getDz()*h_p->getTau(),
        f_p->getX(ix)+U_F[1]*formationTime, f_p->getY(iy)+U_F[2]*formationTime,
        f_p->getZ(iz)+U_F[3]*formationTime, h_p->getTau()+U_F[0]*formationTime, i);
      flux_f[i] += calculateRetardedFriction(f_p->getX(ix), f_p->getY(iy), f_p->getZ(iz),
                  h_p->getTau(), i)/(f_p->getDx()*f_p->getDy()*f_p->getDz()*h_p->getTau());
     } else {
      flux_f[i] += -flux_p[i]-flux_t[i];
     }
    }
   }
   // 2. projectile-fireball friction
   if(ep>0. && ef>0.) {
    double s = mpi*mpi + mN*mN + 2.*mpi*mN*gammaf*gammap*
     (1.0 - vxf*vxp - vyf*vyp - vzf*vzp);
    double Vrel = 0.5/(mN*mpi)*sqrt(pow(s - mN*mN - mpi*mpi,2)
     - 4.*mN*mN*mpi*mpi);
    double D = xi_fa*pow(mN/sqrt(s), 2)*Vrel*xsect->piN(s);
    for(int i=0; i<4; i++){
     flux_pf[i] += D*nbp*(ef + pf)*uf[i]*h_p->getDtau();
    }
    flux_pf[0] += -D*nbp*pf/uf[0]*h_p->getDtau();
   }
   if(et>0. && ef>0.) { // target-fireball friction
    double s = mpi*mpi + mN*mN + 2.*mpi*mN*gammaf*gammat*
     (1.0 - vxf*vxt - vyf*vyt - vzf*vzt);
    double Vrel = 0.5/(mN*mpi)*sqrt(pow(s - mN*mN - mpi*mpi,2)
     - 4.*mN*mN*mpi*mpi);
    double D = xi_fa*pow(mN/sqrt(s), 2)*Vrel*xsect->piN(s);
    for(int i=0; i<4; i++){
     flux_tf[i] += D*nbt*(ef + pf)*uf[i]*h_p->getDtau();
    }
    flux_tf[0] += -D*nbt*pf/uf[0]*h_p->getDtau();
   }
   double taup = h_p->getTau();
   double taut = h_t->getTau();
   double tauf = h_f->getTau();
   // check thermodynamic variables of cell after update
   double _Q_p[7], _Q_t[7], _Q_f[7], Q_p_new[7]={0}, Q_t_new[7]={0}, Q_f_new[7]={0};
   c_p->getQ(_Q_p);
   c_t->getQ(_Q_t);
   c_f->getQ(_Q_f);
   for(int i=0;i<4;i++){
    Q_p_new[i]=_Q_p[i]+(flux_p[i]+flux_pf[i])*taup;
    Q_t_new[i]=_Q_t[i]+(flux_t[i]+flux_tf[i])*taut;
    Q_f_new[i]=_Q_f[i]+(flux_f[i]-flux_pf[i]-flux_tf[i])*tauf;
   }
   double e_p_new, e_t_new, e_f_new, nb_p_new, nb_t_new, nb_f_new, p_new, nq_new, ns_new, vx_new, vy_new, vz_new;
   transformPV(eos,Q_p_new,e_p_new,p_new,nb_p_new,nq_new,ns_new,vx_new,vy_new,vz_new,false);
   transformPV(eos,Q_t_new,e_t_new,p_new,nb_t_new,nq_new,ns_new,vx_new,vy_new,vz_new,false);
   transformPV(eos,Q_f_new,e_f_new,p_new,nb_f_new,nq_new,ns_new,vx_new,vy_new,vz_new,false);
   // compute largest relative change to e-mN*nb across the fluids
   double friction_rel=max((ep-e_p_new-mN*nb_p_new)/(ep-mN*nb_p_new),
                        max((et-e_t_new-mN*nb_t_new)/(et-mN*nb_t_new),
                        (ef-e_f_new-mN*nb_f_new)/(ef-mN*nb_f_new)));
   // overall rescaling factor for all friction terms:
   // if friction_rel<<1, the factor is 1, but it smoothly transitions to a rescaling that
   // will keep the relative change to e-mN*nb below the value specified by MaxRelFriction
   double rescaling=pow(pow(friction_rel/MaxRelFriction,10.0)+1.0,-0.1);
   for(int i=0;i<4;i++){
    flux_f[i]*=rescaling;
    flux_p[i]*=rescaling;
    flux_t[i]*=rescaling;
    flux_pf[i]*=rescaling;
    flux_tf[i]*=rescaling;
   }
   if (_Q_p[0] + (flux_p[0]+flux_pf[0])*taup >= 0 &&
       _Q_t[0] + (flux_t[0]+flux_tf[0])*taut >= 0 &&
       _Q_f[0] + (-flux_pf[0]-flux_tf[0]+flux_f[0])*tauf >= 0) {
    c_p->addFlux((flux_p[0]+flux_pf[0])*taup, (flux_p[1]+flux_pf[1])*taup,
     (flux_p[2]+flux_pf[2])*taup, (flux_p[3]+flux_pf[3])*taup, 0., 0., 0.);
    c_t->addFlux((flux_t[0]+flux_tf[0])*taut, (flux_t[1]+flux_tf[1])*taut,
     (flux_t[2]+flux_tf[2])*taut, (flux_t[3]+flux_tf[3])*taut, 0., 0., 0.);
    c_f->addFlux((-flux_pf[0]-flux_tf[0]+flux_f[0])*tauf, (-flux_pf[1]-flux_tf[1]+flux_f[1])*tauf,
     (-flux_pf[2]-flux_tf[2]+flux_f[2])*tauf, (-flux_pf[3]-flux_tf[3]+flux_f[3])*tauf, 0., 0., 0.);
    c_p->updateByFlux();
    c_t->updateByFlux();
    c_f->updateByFlux();
    c_p->clearFlux();
    c_t->clearFlux();
    c_f->clearFlux();
   }
   if(-flux_p[0]-flux_t[0] > 0. && c_f->getMaxM()<0.01)
    c_f->setAllM(1.0);
   } // end cell loop
 clearRetardedFriction();
 if (decreasingFormTime == 1) {
  formationTime -= dtau * dtauf;
  if (formationTime < 0) formationTime = 0;
 }
}

void MultiHydro::addRetardedFriction(double flux, double x, double y, double z, double t, int i)
{
 vector<double> v = {flux, x, y, z, t, (double)i};
 if (abs(flux) > 1e-3) retardedFriction.push_back(v);
}


double MultiHydro::calculateRetardedFriction(double x, double y, double z, double t, int i)
{
 double total_flux = 0.;
 for (int it = 0; it < retardedFriction.size(); it++)
 {
  if ((abs(retardedFriction[it][1]-x) < dx) && (abs(retardedFriction[it][2]-y) < dy) &&
   (abs(retardedFriction[it][3]-z) < dz) && (abs(retardedFriction[it][4]-t) < dtau) &&
   (retardedFriction[it][5] == (double)i))
   {
    double wCenX = 1-abs(retardedFriction[it][1]-x)/dx;
    double wCenY = 1-abs(retardedFriction[it][2]-y)/dy;
    double wCenZ = 1-abs(retardedFriction[it][3]-z)/dz;
    double wCenT = 1-abs(retardedFriction[it][4]-t)/dtau;
    total_flux += retardedFriction[it][0] * wCenX * wCenY * wCenZ * wCenT;
   }
 }
 return total_flux;
}

void MultiHydro::clearRetardedFriction()
{
 retardedFriction.erase(
  remove_if(
   retardedFriction.begin(),
   retardedFriction.end(),
   [&](vector<double> element){
    if (element[4] < h_p->getTau()) return true;
    else return false;
   }
  ),
  retardedFriction.end()
 );
 // next part is only for printing energy stored in the retardedFriction vector
 /*double storedEnergy = 0.;
 for (int it = 0; it < retardedFriction.size(); it++)
 {
  if (abs(retardedFriction[it][4]-h_p->getTau()) < dtau)
  {
   storedEnergy += retardedFriction[it][0]*(abs(retardedFriction[it][4]-h_p->getTau())/dtau);
  }
  else
  {
   storedEnergy += retardedFriction[it][0];
  }
 }
 cout << "Stored energy: " << storedEnergy << endl;*/
}

void MultiHydro::getEnergyMomentumTensor(double (&T)[4][4], double Q_p[7], double Q_f[7], double Q_t[7])
{
 const double delta[4][4] = {
     {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
 double ep, pp, nbp, nqp, nsp, vxp, vyp, vzp;
 double et, pt, nbt, nqt, nst, vxt, vyt, vzt;
 double ef, pf, nbf, nqf, nsf, vxf, vyf, vzf;
 transformPV(eos, Q_p, ep, pp, nbp, nqp, nsp, vxp, vyp, vzp);
 transformPV(eos, Q_t, et, pt, nbt, nqt, nst, vxt, vyt, vzt);
 transformPV(eos, Q_f, ef, pf, nbf, nqf, nsf, vxf, vyf, vzf);
 // 4-velocities, u_p and u_t
 double gammap = 1.0/sqrt(1.0-vxp*vxp-vyp*vyp-vzp*vzp);
 double up [4] = {gammap,gammap*vxp,gammap*vyp,gammap*vzp};
 double gammat = 1.0/sqrt(1.0-vxt*vxt-vyt*vyt-vzt*vzt);
 double ut [4] = {gammat,gammat*vxt,gammat*vyt,gammat*vzt};
 double gammaf = 1.0/sqrt(1.0-vxf*vxf-vyf*vyf-vzf*vzf);
 double uf [4] = {gammaf,gammaf*vxf,gammaf*vyf,gammaf*vzf};

 // calculation of the energy-momentum tensor
 for (int i=0; i<4; i++){
  for (int j=0; j<4; j++){
   T[i][j] = (ep + pp) * up[i] * up[j] - pp * gmunu[i][j]
    + (et + pt) * ut[i] * ut[j] - pt * gmunu[i][j]
    + (ef + pf) * uf[i] * uf[j] - pf * gmunu[i][j];
  }
 }
}

void MultiHydro::getDiagonalizedEnergyDensity()
{
 double Q_p[7], Q_f[7], Q_t[7];
 double Ttemp[4][4];
 for (int iy = 0; iy < f_p->getNY(); iy++)
  for (int iz = 0; iz < f_p->getNZ(); iz++)
   for (int ix = 0; ix < f_p->getNX(); ix++) {
    Cell *c_p = f_p->getCell(ix, iy, iz);
    Cell *c_t = f_t->getCell(ix, iy, iz);
    Cell *c_f = f_f->getCell(ix, iy, iz);
    c_p->getQ(Q_p);
    c_f->getQ(Q_f);
    c_t->getQ(Q_t);
    for (int i = 0; i < 7; i++) {
     Q_p[i] = Q_p[i]/h_p->getTau();
     Q_t[i] = Q_t[i]/h_t->getTau();
     Q_f[i] = Q_f[i]/h_f->getTau();
    }
    getEnergyMomentumTensor(Ttemp, Q_p, Q_f, Q_t);

    // calculation of the energy-momentum tensor
    TMatrixD T(4,4);
    for (int i=0; i<4; i++)
     for (int j=0; j<4; j++){
      T[i][j] = Ttemp[i][j]*gmunu[j][j];
    }
    if (T[0][0] == 0 && T[1][1] == 0 && T[2][2] == 0 && T[3][3] == 0)
     MHeps[ix][iy][iz] = 0;
    else
    {
     // diagonalization of the energy-momentum tensor
     TMatrixDEigen Te(T);
     TMatrixD eigenValues = Te.GetEigenValues();
     TMatrixD eigenVectors = Te.GetEigenVectors();

     double energyDensity;
     TVectorD v(4);
     for (int i=0; i<4; i++) {
      double vmuvmu = 0;
      energyDensity = eigenValues[i][i];
      v = TMatrixDColumn(eigenVectors,i);
      for (int j=0; j<4; j++) {
       vmuvmu += v[j]*v[j]*gmunu[j][j];
      }
      if (vmuvmu > 0 && energyDensity >= 0) {
       break;
      }
      else if (i == 3) {
       cout << "Multihydro: None of the eigenvectors is time-like, ";
       cout << "using largest eigenvalue for energy density." << endl;
       energyDensity = eigenValues[0][0];
       v = TMatrixDColumn(eigenVectors,0);
       break;
      }
     }

     // save computed energy density into private field
     MHeps[ix][iy][iz] = energyDensity;
    }
   }
}

void MultiHydro::getMaxEnergyDensity()
{
 double Q_p[7], Q_f[7], Q_t[7];
 double Ttemp[4][4];
 for (int iy = 0; iy < f_p->getNY(); iy++)
  for (int iz = 0; iz < f_p->getNZ(); iz++)
   for (int ix = 0; ix < f_p->getNX(); ix++) {
    Cell *c_p = f_p->getCell(ix, iy, iz);
    Cell *c_t = f_t->getCell(ix, iy, iz);
    Cell *c_f = f_f->getCell(ix, iy, iz);
    c_p->getQ(Q_p);
    c_f->getQ(Q_f);
    c_t->getQ(Q_t);
    for (int i = 0; i < 7; i++) {
     Q_p[i] = Q_p[i]/h_p->getTau();
     Q_t[i] = Q_t[i]/h_t->getTau();
     Q_f[i] = Q_f[i]/h_f->getTau();
    }
    double ep, pp, nbp, nqp, nsp, vxp, vyp, vzp;
    double et, pt, nbt, nqt, nst, vxt, vyt, vzt;
    double ef, pf, nbf, nqf, nsf, vxf, vyf, vzf;
    transformPV(eos, Q_p, ep, pp, nbp, nqp, nsp, vxp, vyp, vzp);
    transformPV(eos, Q_t, et, pt, nbt, nqt, nst, vxt, vyt, vzt);
    transformPV(eos, Q_f, ef, pf, nbf, nqf, nsf, vxf, vyf, vzf);

    MHeps[ix][iy][iz] = max(ep, max(et, ef));
   }
}

void MultiHydro::getSumEnergyDensity()
{
 double Q_p[7], Q_f[7], Q_t[7];
 double Ttemp[4][4];
 for (int iy = 0; iy < f_p->getNY(); iy++)
  for (int iz = 0; iz < f_p->getNZ(); iz++)
   for (int ix = 0; ix < f_p->getNX(); ix++) {
    Cell *c_p = f_p->getCell(ix, iy, iz);
    Cell *c_t = f_t->getCell(ix, iy, iz);
    Cell *c_f = f_f->getCell(ix, iy, iz);
    c_p->getQ(Q_p);
    c_f->getQ(Q_f);
    c_t->getQ(Q_t);
    for (int i = 0; i < 7; i++) {
     Q_p[i] = Q_p[i]/h_p->getTau();
     Q_t[i] = Q_t[i]/h_t->getTau();
     Q_f[i] = Q_f[i]/h_f->getTau();
    }
    double ep, pp, nbp, nqp, nsp, vxp, vyp, vzp;
    double et, pt, nbt, nqt, nst, vxt, vyt, vzt;
    double ef, pf, nbf, nqf, nsf, vxf, vyf, vzf;
    transformPV(eos, Q_p, ep, pp, nbp, nqp, nsp, vxp, vyp, vzp);
    transformPV(eos, Q_t, et, pt, nbt, nqt, nst, vxt, vyt, vzt);
    transformPV(eos, Q_f, ef, pf, nbf, nqf, nsf, vxf, vyf, vzf);

    MHeps[ix][iy][iz] = ep + et + ef;
   }
}

void MultiHydro::getEnergyDensity()
{
 // In this function we decide which energy density will be used
 // to find the freezeout hypersurface.
 // - use getMaxEnergyDensity() to use maximum energy of the three fluids
 // - use getDiagonalizedEnergyDensity() to use energy density obtained
 //   by diagonalization of sum of the three energy-momentum tensors
 // getMaxEnergyDensity();
 getDiagonalizedEnergyDensity();
 // getSumEnergyDensity();
}

void MultiHydro::updateEnergyDensity()
{
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    MHepsPrev[ix][iy][iz] = MHeps[ix][iy][iz];
    MHeps[ix][iy][iz] = 0.0;
   }
}

void MultiHydro::outputEnergyDensity()
{
 ofstream fenergy;
 fenergy.open("energy.dat");
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    fenergy << h_p->getTau() << "\t" << f_p->getX(ix) << "\t" << f_p->getY(iy) << "\t" << f_p->getZ(iz) << "\t" << MHeps[ix][iy][iz] << endl;
   }
 fenergy.close();
}

void MultiHydro::resizeMHeps()
{
 double temp[nx][ny][nz];
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    if (2*(ix-(nx-1)/2)+(nx-1)/2 >= 0 && 2*(ix-(nx-1)/2)+(nx-1)/2 < nx && 2*(iy-(ny-1)/2)+(ny-1)/2 >=0 && 2*(iy-(ny-1)/2)+(ny-1)/2 < ny) {
     temp[ix][iy][iz] = MHeps[2*(ix-(nx-1)/2)+(nx-1)/2][2*(iy-(ny-1)/2)+(ny-1)/2][iz];
    } else {
     temp[ix][iy][iz] = 0.;
    }
   }
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    MHeps[ix][iy][iz] = temp[ix][iy][iz];
   }
}

int MultiHydro::findFreezeout(EoS* eosH)
{
 updateEnergyDensity();
 getEnergyDensity();
 swap(eos, eosH);
 int nelements = 0;
 int ne_pos = 0;
 double E=0., Efull = 0.;
 double Ttemp[4][4];

 // allocating corner points for Cornelius
 double ****ccube = new double ***[2];
 for (int i1 = 0; i1 < 2; i1++) {
  ccube[i1] = new double **[2];
  for (int i2 = 0; i2 < 2; i2++) {
   ccube[i1][i2] = new double *[2];
   for (int i3 = 0; i3 < 2; i3++) {
    ccube[i1][i2][i3] = new double[2];
   }
  }
 }

 for (int ix = 2; ix < nx - 2; ix++)
  for (int iy = 2; iy < ny - 2; iy++)
   for (int iz = 2; iz < nz - 2; iz++) {
    double QCube_p[2][2][2][2][7], QCube_f[2][2][2][2][7], QCube_t[2][2][2][2][7];
    // array for storing full energy-momentum tensor of all three fluids at corners
    double TCube[2][2][2][2][4][4];
    double piSquare_p[2][2][2][10], PiSquare_p[2][2][2];
    double piSquare_t[2][2][2][10], PiSquare_t[2][2][2];
    double piSquare_f[2][2][2][10], PiSquare_f[2][2][2];

    // fill all corner cell with energy-momentum tensor
    for (int jx = 0; jx < 2; jx++)
     for (int jy = 0; jy < 2; jy++)
      for (int jz = 0; jz < 2; jz++) {
       ccube[0][jx][jy][jz] = MHepsPrev[ix + jx][iy + jy][iz + jz];
       ccube[1][jx][jy][jz] = MHeps[ix + jx][iy + jy][iz + jz];
       Cell *cc_p = f_p->getCell(ix + jx, iy + jy, iz + jz);
       Cell *cc_t = f_t->getCell(ix + jx, iy + jy, iz + jz);
       Cell *cc_f = f_f->getCell(ix + jx, iy + jy, iz + jz);
       double Qc_p[7], Qc_f[7], Qc_t[7];
       cc_p->getQ(QCube_p[1][jx][jy][jz]);
       cc_f->getQ(QCube_f[1][jx][jy][jz]);
       cc_t->getQ(QCube_t[1][jx][jy][jz]);
       for (int i = 0; i < 7; i++) {
        QCube_p[1][jx][jy][jz][i] = QCube_p[1][jx][jy][jz][i]/h_p->getTau();
        QCube_t[1][jx][jy][jz][i] = QCube_t[1][jx][jy][jz][i]/h_t->getTau();
        QCube_f[1][jx][jy][jz][i] = QCube_f[1][jx][jy][jz][i]/h_f->getTau();
       }
       getEnergyMomentumTensor(Ttemp, QCube_p[1][jx][jy][jz], QCube_f[1][jx][jy][jz], QCube_t[1][jx][jy][jz]);
       for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
         TCube[1][jx][jy][jz][i][j] = Ttemp[i][j];
       }
       cc_p->getQprev(QCube_p[0][jx][jy][jz]);
       cc_f->getQprev(QCube_f[0][jx][jy][jz]);
       cc_t->getQprev(QCube_t[0][jx][jy][jz]);
       for (int i = 0; i < 7; i++) {
        QCube_p[0][jx][jy][jz][i] = QCube_p[0][jx][jy][jz][i]/h_p->getTau();
        QCube_t[0][jx][jy][jz][i] = QCube_t[0][jx][jy][jz][i]/h_t->getTau();
        QCube_f[0][jx][jy][jz][i] = QCube_f[0][jx][jy][jz][i]/h_f->getTau();
       }
       getEnergyMomentumTensor(Ttemp, QCube_p[0][jx][jy][jz], QCube_f[0][jx][jy][jz], QCube_t[0][jx][jy][jz]);
       for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
         TCube[0][jx][jy][jz][i][j] = Ttemp[i][j];
       }
       for (int ii = 0; ii < 4; ii++)
        for (int jj = 0; jj <= ii; jj++) {
         piSquare_p[jx][jy][jz][index44(ii, jj)] = cc_p->getpi(ii, jj);
         piSquare_t[jx][jy][jz][index44(ii, jj)] = cc_t->getpi(ii, jj);
         piSquare_f[jx][jy][jz][index44(ii, jj)] = cc_f->getpi(ii, jj);
       }
       PiSquare_p[jx][jy][jz] = cc_p->getPi();
       PiSquare_t[jx][jy][jz] = cc_t->getPi();
       PiSquare_f[jx][jy][jz] = cc_f->getPi();
    }

    // cornelius
    cornelius->find_surface_4d(ccube);
    const int Nsegm = cornelius->get_Nelements();
    for (int isegm = 0; isegm < Nsegm; isegm++) {
     nelements++;

     // interpolation procedure
     double vxC = 0., vyC = 0., vzC = 0., TC = 0., mubC = 0., muqC = 0.,
            musC = 0., PiC_p = 0., PiC_t = 0., PiC_f = 0., nbC = 0., nqC = 0.;
     double piC_p[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
     double piC_t[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
     double piC_f[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
     double QC_p[7] = {0., 0., 0., 0., 0., 0., 0.};
     double QC_t[7] = {0., 0., 0., 0., 0., 0., 0.};
     double QC_f[7] = {0., 0., 0., 0., 0., 0., 0.};
     double TmunuC[4][4];
     for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
       TmunuC[i][j] = 0;
     }
     double eC = 0., pC = 0.;
     for (int ii = 0; ii < 10; ii++) {
      piC_p[ii] = 0.0;
      piC_t[ii] = 0.0;
      piC_f[ii] = 0.0;
     }
     double wCenT[2] = {1. - cornelius->get_centroid_elem(isegm, 0) / h_p->getDtau(),
                        cornelius->get_centroid_elem(isegm, 0) / h_p->getDtau()};
     double wCenX[2] = {1. - cornelius->get_centroid_elem(isegm, 1) / dx,
                        cornelius->get_centroid_elem(isegm, 1) / dx};
     double wCenY[2] = {1. - cornelius->get_centroid_elem(isegm, 2) / dy,
                        cornelius->get_centroid_elem(isegm, 2) / dy};
     double wCenZ[2] = {1. - cornelius->get_centroid_elem(isegm, 3) / dz,
                        cornelius->get_centroid_elem(isegm, 3) / dz};

     for (int jt = 0; jt < 2; jt++)
      for (int jx = 0; jx < 2; jx++)
       for (int jy = 0; jy < 2; jy++)
        for (int jz = 0; jz < 2; jz++) {
         for (int i = 0; i < 4; i++)
          for (int j = 0; j < 4; j++) {
           TmunuC[i][j] += TCube[jt][jx][jy][jz][i][j] * wCenT[jt] * wCenX[jx] *
                   wCenY[jy] * wCenZ[jz];
         }
         for (int i = 0; i < 7; i++) {
          QC_p[i] += QCube_p[jt][jx][jy][jz][i] * wCenT[jt] * wCenX[jx] *
                   wCenY[jy] * wCenZ[jz];
          QC_t[i] += QCube_t[jt][jx][jy][jz][i] * wCenT[jt] * wCenX[jx] *
                   wCenY[jy] * wCenZ[jz];
          QC_f[i] += QCube_f[jt][jx][jy][jz][i] * wCenT[jt] * wCenX[jx] *
                   wCenY[jy] * wCenZ[jz];
         }
     }

     /*for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
       TmunuC[i][j] = TmunuC[i][j] / (h_p->getTau() - h_p->getDtau() + cornelius->get_centroid_elem(isegm, 0));
     for (int i = 0; i < 7; i++) {
      QC_p[i] = QC_p[i] / (h_p->getTau() - h_p->getDtau() + cornelius->get_centroid_elem(isegm, 0));
      QC_t[i] = QC_t[i] / (h_t->getTau() - h_t->getDtau() + cornelius->get_centroid_elem(isegm, 0));
      QC_f[i] = QC_f[i] / (h_f->getTau() - h_f->getDtau() + cornelius->get_centroid_elem(isegm, 0));
     }*/
     double _ns = 0.0;
     double ep, pp, nbp, nqp, nsp, vxp, vyp, vzp;
     double et, pt, nbt, nqt, nst, vxt, vyt, vzt;
     double ef, pf, nbf, nqf, nsf, vxf, vyf, vzf;
     transformPV(eos, QC_p, ep, pp, nbp, nqp, nsp, vxp, vyp, vzp);
     transformPV(eos, QC_t, et, pt, nbt, nqt, nst, vxt, vyt, vzt);
     transformPV(eos, QC_f, ef, pf, nbf, nqf, nsf, vxf, vyf, vzf);

     TMatrixD T(4,4);
     for (int i=0; i<4; i++)
      for (int j=0; j<4; j++){
       T[i][j] = TmunuC[i][j]*gmunu[j][j];
     }
     // diagonalization of the energy-momentum tensor
     TVectorD v(4);
     if (T[0][0] == 0 && T[1][1] == 0 && T[2][2] == 0 && T[3][3] == 0)
     {
      eC = 0;
      v[0] = 1;
      v[1] = 0;
      v[2] = 0;
      v[3] = 0;
     }
     else
     {
      TMatrixDEigen Te(T);
      TMatrixD eigenValues = Te.GetEigenValues();
      TMatrixD eigenVectors = Te.GetEigenVectors();

      eC = eigenValues[0][0];
      v = TMatrixDColumn(eigenVectors,0);
     }
     vxC = v[1]/v[0];
     vyC = v[2]/v[0];
     vzC = v[3]/v[0];

     //transformPV(eos, QC, eC, pC, nbC, nqC, _ns, vxC, vyC, vzC);
     nbC = nbp + nbt + nbf;
     nqC = nqp + nqt + nqf;
     _ns = nsp + nst + nsf;
     eos->eos(eC, nbC, nqC, _ns, TC, mubC, muqC, musC, pC);
     double TCp, mubCp, muqCp, musCp, pCp;
     double TCt, mubCt, muqCt, musCt, pCt;
     double TCf, mubCf, muqCf, musCf, pCf;
     eos->eos(ep, nbp, nqp, nsp, TCp, mubCp, muqCp, musCp, pCp);
     eos->eos(et, nbt, nqt, nst, TCt, mubCt, muqCt, musCt, pCt);
     eos->eos(ef, nbf, nqf, nsf, TCf, mubCf, muqCf, musCf, pCf);
     if (TC > 0.4 || fabs(mubC) > 0.85) {
      cout << "#### Error (multifluid surface): high T/mu_b (T=" << TC << "/mu_b=" << mubC << ") ####\n";
     }
     for (int jx = 0; jx < 2; jx++)
      for (int jy = 0; jy < 2; jy++)
       for (int jz = 0; jz < 2; jz++) {
        for (int ii = 0; ii < 10; ii++) {
         piC_p[ii] += piSquare_p[jx][jy][jz][ii] * wCenX[jx] * wCenY[jy] * wCenZ[jz];
         piC_t[ii] += piSquare_t[jx][jy][jz][ii] * wCenX[jx] * wCenY[jy] * wCenZ[jz];
         piC_f[ii] += piSquare_f[jx][jy][jz][ii] * wCenX[jx] * wCenY[jy] * wCenZ[jz];
        }
        PiC_p += PiSquare_p[jx][jy][jz] * wCenX[jx] * wCenY[jy] * wCenZ[jz];
        PiC_t += PiSquare_t[jx][jy][jz] * wCenX[jx] * wCenY[jy] * wCenZ[jz];
        PiC_f += PiSquare_f[jx][jy][jz] * wCenX[jx] * wCenY[jy] * wCenZ[jz];
       }
     double v2C = vxC * vxC + vyC * vyC + vzC * vzC;
     if (v2C > 1.) {
      vxC *= sqrt(0.99 / v2C);
      vyC *= sqrt(0.99 / v2C);
      vzC *= sqrt(0.99 / v2C);
      v2C = 0.99;
     }
     double etaC = f_p->getZ(iz) + cornelius->get_centroid_elem(isegm, 3);
     transformToLab(etaC, vxC, vyC, vzC);  // viC is now in lab.frame!
     transformToLab(etaC, vxp, vyp, vzp);
     transformToLab(etaC, vxt, vyt, vzt);
     transformToLab(etaC, vxf, vyf, vzf);
     double gammaC = 1. / sqrt(1. - vxC * vxC - vyC * vyC - vzC * vzC);
     double gammaC_p = 1. / sqrt(1. - vxp * vxp - vyp * vyp - vzp * vzp);
     double gammaC_t = 1. / sqrt(1. - vxt * vxt - vyt * vyt - vzt * vzt);
     double gammaC_f = 1. / sqrt(1. - vxf * vxf - vyf * vyf - vzf * vzf);
     double uC[4] = {gammaC, gammaC * vxC, gammaC * vyC, gammaC * vzC};
     double uC_p[4] = {gammaC_p, gammaC_p * vxp, gammaC_p * vyp, gammaC_p * vzp};
     double uC_t[4] = {gammaC_t, gammaC_t * vxt, gammaC_t * vyt, gammaC_t * vzt};
     double uC_f[4] = {gammaC_f, gammaC_f * vxf, gammaC_f * vyf, gammaC_f * vzf};
     const double tauC = h_p->getTau() - h_p->getDtau() + cornelius->get_centroid_elem(isegm, 0);
     double dsigma[4], dsds;
     // ---- transform dsigma to lab.frame :
     const double ch = cosh(etaC);
     const double sh = sinh(etaC);
     dsigma[0] = tauC * (ch * cornelius->get_normal_elem(0, 0) -
                         sh / tauC * cornelius->get_normal_elem(0, 3));
     dsigma[3] = tauC * (-sh * cornelius->get_normal_elem(0, 0) +
                         ch / tauC * cornelius->get_normal_elem(0, 3));
     dsigma[1] = tauC * cornelius->get_normal_elem(0, 1);
     dsigma[2] = tauC * cornelius->get_normal_elem(0, 2);
     dsds = dsigma[0]*dsigma[0] - dsigma[1]*dsigma[1] - dsigma[2]*dsigma[2] - dsigma[3]*dsigma[3];
     double dVEff_p = 0.0, dVEff_t = 0.0, dVEff_f = 0.0, dVEff = 0.0;
     for (int ii = 0; ii < 4; ii++) {
      dVEff += dsigma[ii] * uC[ii];
      dVEff_p += dsigma[ii] * uC_p[ii];  // normalize for Delta eta=1
      dVEff_t += dsigma[ii] * uC_t[ii];
      dVEff_f += dsigma[ii] * uC_f[ii];
     }
     if (dVEff_p > 0) ne_pos++;
     vEff += dVEff;
     vEff_p += dVEff_p;
     vEff_t += dVEff_t;
     vEff_f += dVEff_f;

     double picart_p[10];
     double picart_t[10];
     double picart_f[10];
#ifdef OUTPI
     // projectile fluid
     /*pi00*/ picart_p[index44(0, 0)] = ch * ch * piC_p[index44(0, 0)] +
                                      2. * ch * sh * piC_p[index44(0, 3)] +
                                      sh * sh * piC_p[index44(3, 3)];
     /*pi01*/ picart_p[index44(0, 1)] =
         ch * piC_p[index44(0, 1)] + sh * piC_p[index44(3, 1)];
     /*pi02*/ picart_p[index44(0, 2)] =
         ch * piC_p[index44(0, 2)] + sh * piC_p[index44(3, 2)];
     /*pi03*/ picart_p[index44(0, 3)] =
         ch * sh * (piC_p[index44(0, 0)] + piC_p[index44(3, 3)]) +
         (ch * ch + sh * sh) * piC_p[index44(0, 3)];
     /*pi11*/ picart_p[index44(1, 1)] = piC_p[index44(1, 1)];
     /*pi12*/ picart_p[index44(1, 2)] = piC_p[index44(1, 2)];
     /*pi13*/ picart_p[index44(1, 3)] =
         sh * piC_p[index44(0, 1)] + ch * piC_p[index44(3, 1)];
     /*pi22*/ picart_p[index44(2, 2)] = piC_p[index44(2, 2)];
     /*pi23*/ picart_p[index44(2, 3)] =
         sh * piC_p[index44(0, 2)] + ch * piC_p[index44(3, 2)];
     /*pi33*/ picart_p[index44(3, 3)] = sh * sh * piC_p[index44(0, 0)] +
                                      ch * ch * piC_p[index44(3, 3)] +
                                      2. * sh * ch * piC_p[index44(0, 3)];
     // target fluid
     /*pi00*/ picart_t[index44(0, 0)] = ch * ch * piC_t[index44(0, 0)] +
                                      2. * ch * sh * piC_t[index44(0, 3)] +
                                      sh * sh * piC_t[index44(3, 3)];
     /*pi01*/ picart_t[index44(0, 1)] =
         ch * piC_t[index44(0, 1)] + sh * piC_t[index44(3, 1)];
     /*pi02*/ picart_t[index44(0, 2)] =
         ch * piC_t[index44(0, 2)] + sh * piC_t[index44(3, 2)];
     /*pi03*/ picart_t[index44(0, 3)] =
         ch * sh * (piC_t[index44(0, 0)] + piC_t[index44(3, 3)]) +
         (ch * ch + sh * sh) * piC_t[index44(0, 3)];
     /*pi11*/ picart_t[index44(1, 1)] = piC_t[index44(1, 1)];
     /*pi12*/ picart_t[index44(1, 2)] = piC_t[index44(1, 2)];
     /*pi13*/ picart_t[index44(1, 3)] =
         sh * piC_t[index44(0, 1)] + ch * piC_t[index44(3, 1)];
     /*pi22*/ picart_t[index44(2, 2)] = piC_t[index44(2, 2)];
     /*pi23*/ picart_t[index44(2, 3)] =
         sh * piC_t[index44(0, 2)] + ch * piC_t[index44(3, 2)];
     /*pi33*/ picart_t[index44(3, 3)] = sh * sh * piC_t[index44(0, 0)] +
                                      ch * ch * piC_t[index44(3, 3)] +
                                      2. * sh * ch * piC_t[index44(0, 3)];
     // projectile fluid
     /*pi00*/ picart_f[index44(0, 0)] = ch * ch * piC_f[index44(0, 0)] +
                                      2. * ch * sh * piC_f[index44(0, 3)] +
                                      sh * sh * piC_f[index44(3, 3)];
     /*pi01*/ picart_f[index44(0, 1)] =
         ch * piC_f[index44(0, 1)] + sh * piC_f[index44(3, 1)];
     /*pi02*/ picart_f[index44(0, 2)] =
         ch * piC_f[index44(0, 2)] + sh * piC_f[index44(3, 2)];
     /*pi03*/ picart_f[index44(0, 3)] =
         ch * sh * (piC_f[index44(0, 0)] + piC_f[index44(3, 3)]) +
         (ch * ch + sh * sh) * piC_f[index44(0, 3)];
     /*pi11*/ picart_f[index44(1, 1)] = piC_f[index44(1, 1)];
     /*pi12*/ picart_f[index44(1, 2)] = piC_f[index44(1, 2)];
     /*pi13*/ picart_f[index44(1, 3)] =
         sh * piC_f[index44(0, 1)] + ch * piC_f[index44(3, 1)];
     /*pi22*/ picart_f[index44(2, 2)] = piC_f[index44(2, 2)];
     /*pi23*/ picart_f[index44(2, 3)] =
         sh * piC_f[index44(0, 2)] + ch * piC_f[index44(3, 2)];
     /*pi33*/ picart_f[index44(3, 3)] = sh * sh * piC_f[index44(0, 0)] +
                                      ch * ch * piC_f[index44(3, 3)] +
                                      2. * sh * ch * piC_f[index44(0, 3)];
#endif

     double dEtotSurf[3] = {0., 0., 0.};
     dEtotSurf[0] = (ep + pCp) * uC_p[0] * dVEff_p - pCp * dsigma[0]; // projectile
     dEtotSurf[1] = (et + pCt) * uC_t[0] * dVEff_t - pCt * dsigma[0]; // target
     dEtotSurf[2] = (ef + pCf) * uC_f[0] * dVEff_f - pCf * dsigma[0]; // fireball
     EtotSurf[0] += dEtotSurf[0];
     EtotSurf[1] += dEtotSurf[1];
     EtotSurf[2] += dEtotSurf[2];
     if (dEtotSurf[0] > 0) EtotSurf_positive[0] += dEtotSurf[0];
     else EtotSurf_negative[0] += dEtotSurf[0];
     if (dEtotSurf[1] > 0) EtotSurf_positive[1] += dEtotSurf[1];
     else EtotSurf_negative[1] += dEtotSurf[1];
     if (dEtotSurf[2] > 0) EtotSurf_positive[2] += dEtotSurf[2];
     else EtotSurf_negative[2] += dEtotSurf[2];

     // exclude segments which fulfills dSigma_0 < 0 & dSigma^2 > 0 - those cells have energy flow into >
     if (dsds > 0) {
      if (dsigma[0] > 0) {
       printFreezeout(
        fmhfreeze_p,
        h_p->getTau() - h_p->getDtau() + cornelius->get_centroid_elem(isegm, 0),
        f_p->getX(ix) + cornelius->get_centroid_elem(isegm, 1),
        f_p->getY(iy) + cornelius->get_centroid_elem(isegm, 2),
        f_p->getZ(iz) + cornelius->get_centroid_elem(isegm, 3),
        dsigma, uC_p, TCp, mubCp, muqCp, musCp, picart_p, PiC_p, dVEff_p
       );
       printFreezeout(
        fmhfreeze_t,
        h_t->getTau() - h_t->getDtau() + cornelius->get_centroid_elem(isegm, 0),
        f_t->getX(ix) + cornelius->get_centroid_elem(isegm, 1),
        f_t->getY(iy) + cornelius->get_centroid_elem(isegm, 2),
        f_t->getZ(iz) + cornelius->get_centroid_elem(isegm, 3),
        dsigma, uC_t, TCt, mubCt, muqCt, musCt, picart_t, PiC_t, dVEff_t
       );
       printFreezeout(
        fmhfreeze_f,
        h_f->getTau() - h_f->getDtau() + cornelius->get_centroid_elem(isegm, 0),
        f_f->getX(ix) + cornelius->get_centroid_elem(isegm, 1),
        f_f->getY(iy) + cornelius->get_centroid_elem(isegm, 2),
        f_f->getZ(iz) + cornelius->get_centroid_elem(isegm, 3),
        dsigma, uC_f, TCf, mubCf, muqCf, musCf, picart_f, PiC_f, dVEff_f
       );
       /*printFreezeout(
        fmhfreeze_all,
        h_f->getTau() - h_f->getDtau() + cornelius->get_centroid_elem(isegm, 0),
        f_f->getX(ix) + cornelius->get_centroid_elem(isegm, 1),
        f_f->getY(iy) + cornelius->get_centroid_elem(isegm, 2),
        f_f->getZ(iz) + cornelius->get_centroid_elem(isegm, 3),
        dsigma, uC, TC, mubC, muqC, musC, picart, PiC, dVEff
       );*/
      }
     } else {
      if (dEtotSurf[0] > 0 && dVEff_p > 0) printFreezeout(
       fmhfreeze_p,
       h_p->getTau() - h_p->getDtau() + cornelius->get_centroid_elem(isegm, 0),
       f_p->getX(ix) + cornelius->get_centroid_elem(isegm, 1),
       f_p->getY(iy) + cornelius->get_centroid_elem(isegm, 2),
       f_p->getZ(iz) + cornelius->get_centroid_elem(isegm, 3),
       dsigma, uC_p, TCp, mubCp, muqCp, musCp, picart_p, PiC_p, dVEff_p
      );
      if (dEtotSurf[1] > 0 && dVEff_t > 0) printFreezeout(
       fmhfreeze_t,
       h_t->getTau() - h_t->getDtau() + cornelius->get_centroid_elem(isegm, 0),
       f_t->getX(ix) + cornelius->get_centroid_elem(isegm, 1),
       f_t->getY(iy) + cornelius->get_centroid_elem(isegm, 2),
       f_t->getZ(iz) + cornelius->get_centroid_elem(isegm, 3),
       dsigma, uC_t, TCt, mubCt, muqCt, musCt, picart_t, PiC_t, dVEff_t
      );
      if (dEtotSurf[2] > 0 && dVEff_f > 0) printFreezeout(
       fmhfreeze_f,
       h_f->getTau() - h_f->getDtau() + cornelius->get_centroid_elem(isegm, 0),
       f_f->getX(ix) + cornelius->get_centroid_elem(isegm, 1),
       f_f->getY(iy) + cornelius->get_centroid_elem(isegm, 2),
       f_f->getZ(iz) + cornelius->get_centroid_elem(isegm, 3),
       dsigma, uC_f, TCf, mubCf, muqCf, musCf, picart_f, PiC_f, dVEff_f
      );
      /*if (dEtotSurf[2] > 0 && dVEff > 0) printFreezeout(
        fmhfreeze_all,
        h_f->getTau() - h_f->getDtau() + cornelius->get_centroid_elem(isegm, 0),
        f_f->getX(ix) + cornelius->get_centroid_elem(isegm, 1),
        f_f->getY(iy) + cornelius->get_centroid_elem(isegm, 2),
        f_f->getZ(iz) + cornelius->get_centroid_elem(isegm, 3),
        dsigma, uC, TC, mubC, muqC, musC, picart, PiC, dVEff
       );*/
     }

    }
 }

 cout << setw(10) << h_p->getTau() << setw(10) << nelements << "\t" << ne_pos << "\t"
      << EtotSurf[0] << "\t" << EtotSurf_positive[0] << "\t" << EtotSurf_negative[0] << "\t"
      << EtotSurf[1] << "\t" << EtotSurf_positive[1] << "\t" << EtotSurf_negative[1] << "\t"
      << EtotSurf[2] << "\t" << EtotSurf_positive[2] << "\t" << EtotSurf_negative[2] << endl;
 swap(eos, eosH); // get back to the hydrodynamic EoS
 for (int i1 = 0; i1 < 2; i1++) {
  for (int i2 = 0; i2 < 2; i2++) {
   for (int i3 = 0; i3 < 2; i3++) {
    delete[] ccube[i1][i2][i3];
   }
   delete[] ccube[i1][i2];
  }
  delete[] ccube[i1];
 }
 delete[] ccube;
 if (nelements == 0 && h_p->getTau() > 5 + tau0)
  return 1;   // particlization surface ended - return 1 for the evolution to stop
 else
  return 0;   // return 0 for the evolution to continue
}

void MultiHydro::printFreezeout(std::ofstream &fout, double t, double x, double y, double z, double dsigma[4], double uC[4], double TC, double mub, double muq, double mus, double picart[10], double PiC, double dVEff)
{
 fout.precision(15);
 fout << setw(24) << t << setw(24) << x << setw(24) << y << setw(24) << z;
 for (int ii = 0; ii < 4; ii++) {
  fout << setw(24) << dsigma[ii];
 }
 for (int ii = 0; ii < 4; ii++) {
  fout << setw(24) << uC[ii];
 }
 fout << setw(24) << TC << setw(24) << mub << setw(24) << muq << setw(24) << mus;
#ifdef OUTPI
 for (int ii = 0; ii < 10; ii++) {
  fout << setw(24) << picart[ii];
 }
 fout << setw(24) << PiC;
#else
 fout << setw(24) << dVEff;
#endif
 fout << endl;
}


double MultiHydro::totalScatRate(double px, double T, double mu, double u[4])
{
 return 2*calculateScatRates(px, T, mu, u, 1) + 3*calculateScatRates(px, T, mu, u, 0);
}


double MultiHydro::calculateScatRates(double px, double T, double mu, double u[4], int particle)
{
 /*
   This function calculates scattering rates of nucleon on nucleons in fluids
 */
 double m1 = 0.938; // scattered particle
 double m2 = particle == 0 ? 0.1396 : 0.938; // scatterers in cloud
 TLorentzVector pLV;
 pLV.SetPxPyPzE(px, 0, 0, sqrt(px*px+m1*m1));
 pLV.Boost(-u[1]/u[0], -u[2]/u[0], -u[3]/u[0]);
 double p = sqrt(pow(pLV.Px(),2) + pow(pLV.Py(),2) + pow(pLV.Pz(),2));
 double gamma = u[0];
 double R;
 double Epi = sqrt(m1*m1+p*p);
 double scrit = m1*m1 + m2*m2 + 2*Epi*m2;
 double ds = 0.1;
 double smin = pow(m1 + m2,2);
 double smax = 10;
 double g = 2;
 for (double s = smin + ds/2; s < smax; s += ds)
 {
  //double Ekin = s/(2.0*mN) - 2.0*mN;
  //double pLab = sqrt(s*(s-4*mN*mN)) / (2.0*mN);
  double sa = pow(m1+m2,2);
  double sb = pow(m1-m2,2);
  double sprime = s - m1*m1 - m2*m2;
  double ebar, deltaE;
  if (s <= scrit) {
   ebar = 0.5*((Epi*sprime+p*sqrt(sprime*sprime-4*m1*m1*m2*m2))/(2*m1*m1)+2*Epi*m2*m2/sprime);
   deltaE = 0.5*((Epi*sprime+p*sqrt(sprime*sprime-4*m1*m1*m2*m2))/(2*m1*m1)-2*Epi*m2*m2/sprime);
  } else {
   ebar = Epi*sprime/(2*m1*m1);
   deltaE = p*sqrt(sprime*sprime-4*m1*m1*m2*m2)/(2*m1*m1);
  }
  double dR = sqrt(s-sa)*sqrt(s-sb)*exp(-ebar/T)*sinh(deltaE/T)/pow(0.197,3);
  // case of proton-proton scattering
  if (particle == 1) dR *= 0.1 * pp_total(s);
  // case of pion-proton scattering
  if (particle == 0) dR *= xsect->piN(s);
  //cout << s << " " << exp(ebar/T) << " " << exp(-sqrt(m1*m1+p*p)*(s - m1*m1 - m2*m2)/(2*T*m1*m1)) << " " << dR << endl;
  if (dR > 1e-20 && dR < 1e20) R += dR * ds;
 }
 R *= g*T*exp(mu/T)/(8*pow(M_PI,2)*p*Epi);
 R *= gamma;
 //exit(1);
 return R;
}

double MultiHydro::pp_total(double mandelstam_s) {
  double mN = 0.938;
  const double p_lab = sqrt(mandelstam_s*(mandelstam_s-4*mN*mN)) / (2.0*mN);
  if (p_lab < 0.4) {
    return 34 * pow(p_lab / 0.4, -2.104);
  } else if (p_lab < 0.8) {
    return 23.5 + 1000 * pow(p_lab - 0.7, 4);
  } else if (p_lab < 1.5) {
    return 23.5 + 24.6 / (1 + exp(-(p_lab - 1.2) / 0.1));
  } else if (p_lab < 5.0) {
    return 41 + 60 * (p_lab - 0.9) * exp(-1.2 * p_lab);
  } else {
    const auto logp = log(p_lab);
    return 48.0 + 0.522 * logp * logp - 4.51 * logp;
  }
}

double MultiHydro::Fermi(double nb) {
 return pow(0.197,2)*pow(6*nb/M_PI,2/3)*M_PI*M_PI/(2*0.938);
}

double MultiHydro::getLocalEnergyDensity(double x, double y, double eta) {
// this function returns effective energy density, pre-computed in the MHeps[][][] array,
// at an arbitrary position using tri-linear interpolation over the table(array)
 int ix = (int)((x - f_p->getX(0))/dx);
 int iy = (int)((y - f_p->getY(0))/dy);
 int iz = (int)((eta - f_p->getZ(0))/dz);
 if(ix<0 || ix>nx-1 || iy<0 || iy>ny-1 || iz<0  || iz>nz-1) {
  cout << "MultiHydro::getLocalEnergyDensity: out of bounds\n";
  exit(77);
 }
 double xm = (x - f_p->getX(0) - ix*dx) / dx;
 double ym = (y - f_p->getY(0) - iy*dy) / dy;
 double zm = (eta - f_p->getZ(0) - iz*dz) / dz;
 double wx [2] = {1. - xm, xm};
 double wy [2] = {1. - ym, ym};
 double wz [2] = {1. - zm, zm};
 double interpolatedEps = 0.;
 for(int i=0; i<2; i++)
 for(int j=0; j<2; j++)
 for(int k=0; k<2; k++)
  interpolatedEps += wx[i]*wy[j]*wz[k]*MHeps[ix+i][iy+j][iz+k];
 return interpolatedEps;
}

void MultiHydro::evolveSpectators(void) {
 for(int iev=0; iev<nucleons.size(); iev++) {
 for(int i=0; i<nucleons[iev].size(); i++) {
  nucleons[iev][i].eta = nucleons[iev][i].rap +
   asinh(h_f->getTau()/(h_f->getTau()+dtau)*sinh(nucleons[iev][i].eta - nucleons[iev][i].rap));
 }
 }
 for(int iev=0; iev<nucleons.size(); iev++) {
  for (auto it = nucleons[iev].begin(); it != nucleons[iev].end();) {
   if(getLocalEnergyDensity((*it).x, (*it).y, (*it).eta)>ecrit) {
    it = nucleons[iev].erase(it);
   } else
   it++;
  }
 }
}

void MultiHydro::printSpectators(std::ofstream &fout) {
 for(int iev=0; iev<nucleons.size(); iev++) {
  fout << "#event " << iev << endl;
  for(int i=0; i<nucleons[iev].size(); i++) {
  // id offset for the output is big enough so that there is no overlap
  // with ids of hadrons, generated by smash-hadron-sampler
  const int id = i + 50000;
  Nucleon &n = nucleons[iev][i];
  const int pdg = n.charge==1 ? 2212 : 2112;
  const double mass = 0.938;
  fout << tau0*cosh(n.getEta0()) << " " << n.x
   << " " << n.y << " " << tau0*sinh(n.getEta0())
   << " " << mass << " " << mass*cosh(n.rap)
   << " " << 0 << " " << 0 << " " << mass*sinh(n.rap)
   << " " << pdg << " " << id << " " << n.charge << endl;
  }
 }
 fout << "#event end" << endl;
}
