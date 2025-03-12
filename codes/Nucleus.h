/*************************************************************************
 *
 *  Filename: Nucleus.h
 *
 *		Class for creating a Nucleus
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
//#include "UncertainNumber.h"

// Need mass table

using namespace std;

class Nucleus {
private :


 // 931.4940090

	// Liquid drop model (J.W. Rohlf 94)
	double av		=	15.75;
	double as		=	17.8;
	double ac		=	0.711;
	double aa		=	23.7;
	double del_ee	= +11.18;
	double del_eo	=	-11.18;

protected :

public :

	UncertainNumber mamu;	// Mass in AMU
	UncertainNumber m;	// Mass in MeV/c2

	int		  A;				  // Mass number
	int		  Z;				  // Atomic number
	double	Ex;				  // Excitation energy (MeV)
	double  r;					// Nuclear radius (fm)
	double	massExcess;	// Mass excess (MeV)
	double	J;				  // Spin
	int		  P;				  // Parity
	string 	symbol;			// Element symbol

	Nucleus(int z , int a) : A(a), Z(z)
	{
		// ------------------------------------
		// Find nucleus in NIST tables
		// ------------------------------------
		string line, sline;
		int Ai, Zi, index = 0;
		ifstream fip("AtomicWeight.txt");
		for (int i = 0; i < 6; i++) {getline(fip, line);}		// Header

		while(getline(fip, line))
		{
			if (line == "_________________________________________________________________________")
			{
					getline(fip, line);
					if (strlen(line.c_str()) > 11)
					{
						Zi = atoi(line.substr(0,3).c_str());
						Ai = atoi(line.substr(6,11).c_str());
						if (Zi == Z)
						{
							symbol = line.substr(3,5);
							while (line != "_________________________________________________________________________")
							{
								if (Ai == A)
								{
									sline = line.substr(11,32);
									mamu.value = atof(sline.substr(0,line.find(".")).c_str());
								}
								getline(fip, line);
								Ai = atoi(line.substr(6,11).c_str());
							}
						}
					}
			}
		}

		fip.close();
		m = mamu*AMU;
		r = 1.25*pow(((double)A),(1./3.));
	};

	~Nucleus() {};

	void 	printSummary()
	{
		cout << A << symbol << Z << endl;
		cout << mamu.toString() << endl;
		cout << m.toString() << endl;
	};

	/* --------------------------------
		Semiemperical mass formula
	   --------------------------------
	*/

	//double SEM_volume() const {return av*A;}
	//double SEM_surface() const {return as*pow(A, 2./3.);}
	//double SEM_coulomb() const {return ac*(pow(Z,2)/pow(A,1/3);}
	//double getLD_Eb

};
