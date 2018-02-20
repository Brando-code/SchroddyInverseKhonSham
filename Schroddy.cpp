/*
 * Schroddy.cpp
 *
 *  Created on: 17 feb 2018
 *      Author: brando
 */
//Class for Schrodinger equation solver
# include <iostream>
# include "Potentials.h"
# include "Parameters.h"
# include "Schroddy.h"
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
    return Parameters::hbar_omega*(2*(m_n-1)+m_l+(3/2));
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

std::vector<double> Schroddy::solveSchroddyByRK(double x0, double x1, double psi0, double psiPrime0) const
{
	double eigenValLower = 0.0;
	double eigenValUpper = 100.0;
	double err = 10e-8;
	int End_sign = -1;
	bool Limits_are_defined = false;
	double trialEigenvalue;
	double normalCoef;
	unsigned int nodes;
	int step = 1000 /*Parameters::stepH*/,Nstep = (x1-x0)/step;
	const double factor = 2*Parameters::mn/(Parameters::hbarc*Parameters::hbarc);
	double runningX = x0, runningPsi = psi0, runningPsiPrime = psiPrime0, psiR = 0.0 /*Parameters::psiRight*/;
	std::vector<double> psiArray;
	std::vector<double> eigenEnergies;
	std::vector<double> normalPsi;
	psiArray.push_back(psi0);


	for(unsigned int i = 1; i <= m_nQuant; ++i)
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
				const double k1 = step*runningPsiPrime;
				const double l1 = factor*(m_pot -> potential(runningX) - trialEigenvalue)*runningPsi;
				const double k2 = step*(runningPsiPrime + 1./2.*l1);
				const double l2 = factor*(m_pot -> potential(runningX + step/2.) - trialEigenvalue)*(runningPsi + 1./2.*k1);
				const double k3 = step*(runningPsiPrime + 1./2.*l2);
				const double l3 = factor*(m_pot -> potential(runningX + step/2.) - trialEigenvalue)*(runningPsi + 1./2.*k2);
				const double k4 = step*(runningPsiPrime + l3);
				const double l4 = factor*(m_pot -> potential(runningX + step) - trialEigenvalue)*(runningPsi + k3);

				//advance running variables and store intermediate results
				//runningPsi += step/6.*(k1 + 2*k2 + 2*k3 + k4);
				runningPsi += 1./6.*(k1 + 2*k2 + 2*k3 + k4);
				//runningPsiPrime += step/6.*(l1 + 2*l2 + 2*l3 + l4);
				runningPsiPrime += 1./6.*(l1 + 2*l2 + 2*l3 + l4);
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
				const double k1 = step*runningPsiPrime;
				const double l1 = factor*(m_pot -> potential(runningX) - trialEigenvalue)*runningPsi;
				const double k2 = step*(runningPsiPrime + 1./2.*l1);
				const double l2 = factor*(m_pot -> potential(runningX + step/2.) - trialEigenvalue)*(runningPsi + 1./2.*k1);
				const double k3 = step*(runningPsiPrime + 1./2.*l2);
				const double l3 = factor*(m_pot -> potential(runningX + step/2.) - trialEigenvalue)*(runningPsi + 1./2.*k2);
				const double k4 = step*(runningPsiPrime + l3);
				const double l4 = factor*(m_pot -> potential(runningX + step) - trialEigenvalue)*(runningPsi + k3);

				//advance running variables and store intermediate results
				//runningPsi += step/6.*(k1 + 2*k2 + 2*k3 + k4);
				runningPsi += 1./6.*(k1 + 2*k2 + 2*k3 + k4);
				//runningPsiPrime += step/6.*(l1 + 2*l2 + 2*l3 + l4);
				runningPsiPrime += 1./6.*(l1 + 2*l2 + 2*l3 + l4);
				runningX += step;
				psiArray.push_back(runningPsi/runningX);
			}

			if (End_sign*psiArray[Nstep-1] > psiR)
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
		for (int i=1; i <= step; ++i)
		{
			// Simple integration
			Integral += 0.5*Nstep*(psiArray[i-1]*psiArray[i-1]+psiArray[i]*psiArray[i]);
		}

		normalCoef = sqrt(1.0/Integral);
		//double normalPsi;
		for (int i=0; i <= step; ++i)
		{
			normalPsi.push_back(normalCoef*psiArray[i]);
		}
	}

    std::vector <double>::iterator walk1 = eigenEnergies.begin();
    std::vector <double>::iterator walk2 = normalPsi.begin();
    while (walk1 != eigenEnergies.end() && walk2 != normalPsi.end())
    {
    	//file << *walk1 << "\t\t" << *walk2 << std::endl;
    	std::cout << *walk1 << "\t\t" << *walk2 << std::endl;
    	walk1++;
    	walk2++;
    }return normalPsi;
}

std::vector<double> Schroddy::solveSchroddyByNV(double psiL, double psiR) const
{
	double eigenValLower = 0.0;
	double eigenValUpper = 10.0;
	double err = 1e-10;
	int End_sign = -1;
	bool Limits_are_defined = false;
	double trialEigenvalue;
	double normalCoef;
	unsigned int nodes;
	double Ksquare;
	int step = 1000 /*Parameters::stepH*/;
	//const double factor = 2*Parameters::mn/(Parameters::hbarc*Parameters::hbarc);
	std::vector<double> eigenEnergies(m_nQuant+1);
	std::vector<double> normalPsi(step+1);
	std::vector<double> psi (step+1);
	std::vector<double> potential(step+1);
	psi[0]=psiL;
	psi[1]=psiL+1.e-3;
	const static double xRange = 12;
	const static double h0 = xRange / step;

	//double poten = 0;
	double dist;
	for (int i=0; i<= step; ++i)
	{
		dist = i*h0 - 0.5*xRange;
		//poten = m_pot -> potential(dist);
		//potential.push_back(poten);
		//potential.push_back(0.5*dist*dist);
		potential[i] = 0.5*dist*dist;

	}

	for(unsigned int i = 1; i <= m_nQuant; ++i)
	{
			Limits_are_defined = false;
			while (Limits_are_defined == false)
			{
				/* First, determine an upper limit for energy, so that the wave-
				function Psi[i] has one node more than physically needed.
				 */
				nodes = 0;
				trialEigenvalue = eigenValUpper;

				for (int i=2; i <= step; ++i)
				{
					Ksquare = 2.0*(trialEigenvalue - potential[i])/**factor*/;
					psi[i] = 2.0*psi[i-1]*(1.0 - (5.0*h0*h0*Ksquare / 12.0))
					/(1.0 + (h0*h0*Ksquare/12.0))-psi[i-2];

					if (psi[i]*psi[i-1] < 0)
						++nodes;
				}

				/* If one runs into the following condition, the modification
				of the upper limit was too aggressive. */
				if (eigenValUpper < eigenValLower)
					eigenValUpper = MAX(2*eigenValUpper, -2*eigenValUpper);
				if (nodes > m_nQuant) eigenValUpper *= 0.7;
					else if (nodes < m_nQuant) eigenValUpper *= 2.0;
				else Limits_are_defined = true; // At least one node should appear.
			}

			// Refine the energy by satisfying the right boundary condition.
			End_sign = -End_sign;
			while ((eigenValUpper - eigenValLower) > err)
			{
				trialEigenvalue = (eigenValUpper + eigenValLower) / 2.0;
				for (int i=2; i <= step; ++i)
				{
					Ksquare = 2.0*(trialEigenvalue - potential[i])/**factor*/;
					psi[i] = 2.0*psi[i-1]*(1.0 - (5.0*h0*h0*Ksquare / 12.0))
					/(1.0 + (h0*h0*Ksquare/12.0))-psi[i-2];
				}

				if (End_sign*psi[step] > psiR)
					eigenValLower = trialEigenvalue;
				else eigenValUpper = trialEigenvalue;
			}

			// Initialization for the next iteration in main loop
			trialEigenvalue = (eigenValUpper + eigenValLower)/2;
			//eigenEnergies.push_back(trialEigenvalue);
			eigenEnergies[m_nQuant] = trialEigenvalue;
			eigenValUpper = trialEigenvalue;
			eigenValLower = trialEigenvalue;

			// Now find the normalization coefficient
			double Integral = 0.0;
			for (int i=1; i <= step; ++i)
			{
				// Simple integration
				Integral += 0.5*h0*(psi[i-1]*psi[i-1]+psi[i]*psi[i]);
			}

			normalCoef = sqrt(1.0/Integral);
			for (int i=0; i <= step; ++i)
			{
				//normalPsi.push_back(normalCoef*psi[i]);
				normalPsi[i] = normalCoef*psi[i];
			}
		}

	    std::vector <double>::iterator walk1 = eigenEnergies.begin();
	    std::vector <double>::iterator walk2 = normalPsi.begin();
	    while (walk1 != eigenEnergies.end() && walk2 != normalPsi.end())
	    {
	    	//file << *walk1 << "\t\t" << *walk2 << std::endl;
	    	std::cout << *walk1*3.5 << "\t\t" << *walk2 << std::endl;
	    	walk1++;
	    	walk2++;
	    }return normalPsi;
}













