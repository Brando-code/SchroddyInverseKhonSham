/*
 * Schroddy.cpp
 *
 *  Created on: 17 feb 2018
 *      Author: brando
 */
//Class for Schrodinger equation solver
//# include <math.h>
//# include <stdlib.h>
//# include <stdio.h>
# include <iostream>
# include "Potentials.h"
# include "Parameters.h"
# include "Schroddy.h"
//# include <cmath>
# include <vector>
# include <fstream>
# define MAX(a, b) (((a) > (b)) ? (a) : (b))


Eigenvalues::~Eigenvalues() {}

/*=======================================================================
 * HO Eigenvalues generator
 *=====================================================================*/

HarmonicEigenvalues::HarmonicEigenvalues(unsigned int n, int l): m_n(n), m_l(l) {}
double HarmonicEigenvalues::eigenvalue() const
{
    return Parameters::hbar_omega*(2*m_n+m_l+(3/2));
}
Eigenvalues* HarmonicEigenvalues::clone() const
{
    return new HarmonicEigenvalues(*this);
}

/*========================================================================
 * Schrodinger equation solver by Runge-Kutta and Shooting method
 *======================================================================*/
Schroddy::Schroddy(const InitialPot& pot, unsigned int nQuant): m_pot(pot.clone()), m_nQuant(nQuant) {}

Schroddy::~Schroddy()
{
    delete m_pot;
}

double Schroddy::solveSchroddyByRK(double x0, double x1, double psi0, double psiPrime0) const
{
	double eigenValLower = 0.0;
	double eigenValUpper = 100.0;
	double err = 10e-8;
	int End_sign = -1;
	bool Limits_are_defined = false;
	double trialEigenvalue;
	double normalCoef;
	unsigned int nodes, step = 1000 /*Parameters::stepH*/, Nstep = (x1-x0)/step;
	const double factor = 2*Parameters::mn/(Parameters::hbarc*Parameters::hbarc);
	double runningX = x0, runningPsi = psi0, runningPsiPrime = psiPrime0, psiR = 0.0 /*Parameters::psiRight*/;
	std::vector<double> psiArray;
	std::vector<double> eigenEnergies;
	std::vector<double> normalPsi;
	psiArray.push_back(psi0);

	for(int i = 1; i <= m_nQuant; ++i)
	{
		Limits_are_defined = false;
		while (Limits_are_defined == false)
		{
			/* First, determine an upper limit for energy, so that the wave-
			function Psi[i] has one node more than physically needed.
			 */
			nodes = 0;
			trialEigenvalue = eigenValUpper;
			for (int i=1; i <= Nstep; ++i)
			{
				//compute decoupled RK factors
				const double k1 = runningPsiPrime;
				const double l1 = (m_pot -> potential(runningX) - trialEigenvalue)*runningPsi;
				const double k2 = runningPsiPrime + step/2.*l1;
				const double l2 = (m_pot -> potential(runningX + step/2.) - trialEigenvalue)*(runningPsi + step/2.*k1);
				const double k3 = runningPsiPrime + step/2.*l2;
				const double l3 = (m_pot -> potential(runningX + step/2.) - trialEigenvalue)*(runningPsi + step/2.*k2);
				const double k4 = runningPsiPrime + step*l3;
				const double l4 = (m_pot -> potential(runningX + step) - trialEigenvalue)*(runningPsi + step*k3);

				//advance running variables and store intermediate results
				runningPsi += step/6.*factor*(k1 + 2*k2 + 2*k3 + k4);
				runningPsiPrime += step/6.*factor*(l1 + 2*l2 + 2*l3 + l4);
				runningX += step;
				psiArray.push_back(runningPsi/runningX);

				if (psiArray[i]*psiArray[i-1] < 0)
					++nodes;
			}

			/* If one runs into the following condition, the modification
			of the upper limit was too aggressive. */
			if (eigenValUpper < eigenValLower)
				eigenValUpper = MAX(2*eigenValUpper, -2*eigenValUpper);
			if (nodes > m_nQuant) eigenValUpper *= 0.7;
				else if (nodes < m_nQuant) eigenValUpper *= 2.0;
			else Limits_are_defined = true; // At least one node should appear.
		}// End of the loop: while (Limits_are_defined == false)

		// Refine the energy by satisfying the right boundary condition.
		End_sign = -End_sign;
		while ((eigenValUpper - eigenValLower) > err)
		{
			trialEigenvalue = (eigenValUpper + eigenValLower) / 2.0;
			for (int i=2; i <= Nstep; ++i)
			{
				//compute decoupled RK factors
				const double k1 = runningPsiPrime;
				const double l1 = (m_pot -> potential(runningX) - trialEigenvalue)*runningPsi;
				const double k2 = runningPsiPrime + step/2.*l1;
				const double l2 = (m_pot -> potential(runningX + step/2.) - trialEigenvalue)*(runningPsi + step/2.*k1);
				const double k3 = runningPsiPrime + step/2.*l2;
				const double l3 = (m_pot -> potential(runningX + step/2.) - trialEigenvalue)*(runningPsi + step/2.*k2);
				const double k4 = runningPsiPrime + step*l3;
				const double l4 = (m_pot -> potential(runningX + step) - trialEigenvalue)*(runningPsi + step*k3);

				//advance running variables and store intermediate results
				runningPsi += step/6.*factor*(k1 + 2*k2 + 2*k3 + k4);
				runningPsiPrime += step/6.*factor*(l1 + 2*l2 + 2*l3 + l4);
				runningX += step;
				psiArray.push_back(runningPsi/runningX);
			}
			if (End_sign*psiArray[Nstep] > psiR)
				eigenValLower = trialEigenvalue;
			else eigenValUpper = trialEigenvalue;
		}// End of loop: while ((E_upperLimit - E_lowerLimit) > Epsilon)

		// Initialization for the next iteration in main loop
		trialEigenvalue = (eigenValUpper+eigenValLower)/2;
		eigenEnergies.push_back(trialEigenvalue);
		eigenValUpper = trialEigenvalue;
		eigenValLower = trialEigenvalue;

		// Now find the normalization coefficient
		double Integral = 0.0;
		for (int i=1; i <= Nstep; ++i)
		{
			// Simple integration
			Integral += 0.5*Nstep*(psiArray[i-1]*psiArray[i-1]+psiArray[i]*psiArray[i]);
		}

		normalCoef = sqrt(1.0/Integral);
		//double normalPsi;
		for (int i=0; i <= Nstep; ++i)
		{
			normalPsi.push_back(normalCoef*psiArray[i]);
		}
	} return 0;//return normalPsi;

    std::vector <double>::iterator walk1 = eigenEnergies.begin();
    std::vector <double>::iterator walk2 = normalPsi.begin();
    while (walk1 != eigenEnergies.end() && walk2 != normalPsi.end())
    {
    	//file << *walk1 << "\t\t" << *walk2 << std::endl;
    	std::cout << *walk1 << "\t\t" << *walk2 << std::endl;
    	walk1++;
    	walk2++;
    }
}





