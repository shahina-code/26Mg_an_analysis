/*************************************************************************
 *
 *  Filename: Reactions.h
 *
 *  Description:
 *		Class for nuclear reactions
 *
 *	Author(s):
 *     Michael T. Febbraro
 *
 *  Creation Date: 7/2/2016
 *  Last modified: 3/31/2019
 *
 * -----------------------------------------------------
 * 	Physics Divsion - Oak Ridge National Laboratory
 *	Oak Ridge, TN, USA
 *
 */
 
#include <iostream>
#include <fstream>
#include <ctime>
#include <signal.h>
#include <string>
#include <cstring>
//#include "PhysicsConstants.h"
//#include "Nucleus.h"
//#include "UncertainNumber.h"
using namespace std;

class Reaction {
private :

  double pi = 4.0*TMath::ATan(1.0);

	double	calRutherford();
	double	cal2body();
	double	calCoulBarrier();


protected :

public :

	Nucleus		projectile;
	Nucleus		target;
	Nucleus		ejectile;
	Nucleus		recoil;

	UncertainNumber	qval;
	UncertainNumber coulBar;
	UncertainNumber ebeam, etot, ecm, ex;
	UncertainNumber e3_cm, e4_cm, p3_cm, p4_cm;
	UncertainNumber e3, e4, p3, p4;

	UncertainNumber A, B, C, D;

	Reaction (Nucleus A, Nucleus b, Nucleus c, Nucleus D) : target (A), projectile (b), ejectile (c), recoil(D)
  {
    // Check for charge and nucleon conservation
    if (projectile.Z + target.Z != ejectile.Z + recoil.Z) {cout << "WARNING - charge not conserved!" << endl;}
    else if (projectile.A + target.A != ejectile.A + recoil.A) {cout << "WARNING - number of nucleons not conserved!" << endl;}

  	// Calculate q-value
    qval = projectile.m + target.m - ejectile.m - recoil.m;

  	// Calculate coul barrier
    coulBar.value = 1.4399643929*((double)(projectile.Z + target.Z))/(projectile.r + target.r);

    // Calculate kinematics
    ebeam.set(0.0, 0.0);
    etot = ebeam + target.m + projectile.m;

  	// Calculate distance of closest approach

  };

	~Reaction(){};

  void printSummary()
  { cout << "________________________________________________________" << endl;
    cout << " Reaction: " << target.A << target.symbol << "(" << projectile.A << projectile.symbol << "," << ejectile.A << ejectile.symbol << ")" << recoil.A << recoil.symbol;
    cout << "Q-value: " << qval.value << " MeV" << endl;
    cout << "Coulomb barrier: " << coulBar.value << " MeV" << endl;
    cout << "Beam energy: " << ebeam.value << " MeV" << endl;
    cout << "CM energy: " << ecm.value << " MeV" << endl;
    cout << "________________________________________________________" << endl;
  };

	UncertainNumber	getQval() {return qval;};
	UncertainNumber getCoulBar() {return coulBar;};

  void setExcitation(double e, double err)
  {
    ex.set(e, err);
    qval = qval - ex;
  }

	void	setEbeam(double e, double err)
  {
    ebeam.set(e, err);
    etot = ebeam + qval;
    ecm = (target.m / (projectile.m + target.m))*ebeam;

    A = (projectile.m*recoil.m*(ebeam/etot))/((projectile.m + target.m)*(ejectile.m + recoil.m));
    B = (projectile.m*ejectile.m*(ebeam/etot))/((projectile.m + target.m)*(ejectile.m + recoil.m));
    C = ((target.m*ejectile.m)/((projectile.m + target.m)*(ejectile.m + recoil.m)))*((UncertainNumber)1 + ((projectile.m*qval)/(target.m*etot)));
    D = ((target.m*recoil.m)/((projectile.m + target.m)*(ejectile.m + recoil.m)))*((UncertainNumber)1 + ((projectile.m*qval)/(target.m*etot)));

  };

  void  setEbeam(UncertainNumber e)
  {
    ebeam = e;
    etot = ebeam + qval;
    ecm = (target.m / (projectile.m + target.m))*ebeam;

    A = (projectile.m*recoil.m*(ebeam/etot))/((projectile.m + target.m)*(ejectile.m + recoil.m));
    B = (projectile.m*ejectile.m*(ebeam/etot))/((projectile.m + target.m)*(ejectile.m + recoil.m));
    C = ((target.m*ejectile.m)/((projectile.m + target.m)*(ejectile.m + recoil.m)))*((UncertainNumber)1 + ((projectile.m*qval)/(target.m*etot)));
    D = ((target.m*recoil.m)/((projectile.m + target.m)*(ejectile.m + recoil.m)))*((UncertainNumber)1 + ((projectile.m*qval)/(target.m*etot)));
  };

  double getEjectileElab_thcm(double thetaCM)
  {
    double cth_cm = TMath::Cos(TMath::DegToRad()*thetaCM);

    e3 = etot*(B + D + (UncertainNumber)2*((A*C)^0.5)*(UncertainNumber)cth_cm);
    return e3.getValue();
  };

  double getEjectileElab_thlab(double thetaLAB)
  {
    double cth_lab = TMath::Cos(TMath::DegToRad()*thetaLAB);
    double sth_lab2 = TMath::Power(TMath::Sin(TMath::DegToRad()*thetaLAB), 2);

    if (B.getValue() < D.getValue())
      e3 = etot*B*(((UncertainNumber)cth_lab + (((D/B) - (UncertainNumber)sth_lab2)^0.5))^2.);

    return e3.getValue();
  };

  double getEjectileThcm_thlab(double thetaLAB)
  {
    double cth_lab = TMath::Cos(TMath::DegToRad()*thetaLAB);
    double sth_lab2 = TMath::Power(TMath::Sin(TMath::DegToRad()*thetaLAB), 2);

    if (B.getValue() < D.getValue())
      e3 = etot*B*(((UncertainNumber)cth_lab + (((D/B) - (UncertainNumber)sth_lab2)^0.5))^2.);

    return TMath::RadToDeg()*TMath::ASin(((e3/etot)/D).getValue()*TMath::Sin(TMath::DegToRad()*thetaLAB));
  };


};
