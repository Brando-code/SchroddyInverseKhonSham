/*
 * Potentials.cpp
 *
 *  Created on: 17 feb 2018
 *      Author: brando
 */
//Class for initial potentials
#include <cmath>
#include "Parameters.h"
#include "Potentials.h"

InitialPot::~InitialPot() {}

/*===================================================================
 * Woods-Saxon potential
 *=================================================================*/

WSaxPot::WSaxPot(double V0, double Rn, double a0): m_V0(V0), m_Rn(Rn), m_a0(a0) {};
double WSaxPot::potential(double x) const
{
    double wspot = -m_V0/(1+exp ((x-m_Rn)/m_a0));
    return wspot;
}
InitialPot* WSaxPot::clone() const
{
    return new WSaxPot(*this); //return a derived class object through a base class pointer
}


/*===================================================================
 * HO potential
 *=================================================================*/

HOPot::HOPot(double m, int anglmomentum): m_m(m), m_anglmomentum(anglmomentum){}
double HOPot::potential(double x) const
{
    double angularPart = 0;
    if (x != 0)
        angularPart = (Parameters::hbarc*Parameters::hbarc)*m_anglmomentum*(m_anglmomentum+1)/(2*m_m*x*x);

    const double c =(m_m*Parameters::hbar_omega*Parameters::hbar_omega)/(Parameters::hbarc*Parameters::hbarc);
    double hopot = angularPart+0.5*c*x*x;

    return hopot;
}
/*
double HOPot::getOmega() const
{
    return m_omega;
}*/

InitialPot* HOPot::clone() const
{
    return new HOPot(*this); //return a derived class object through a base class pointer
}

//=====================================================================
// Spin-orbit potential
//=====================================================================
/*
double SOPot(double k0, doubler0, double x, double hbar, double Rn, double a)
{
	double l, j;
	double LS=((hbar*hbar)/2)(j(j+1)-l(l+1)-3/4)
	double deriv
	sopot=*/

/*=====================================================================
 * Kohn-Sham potential
 *===================================================================*/

potOut::potOut(double potKS): m_potKS(potKS){}
double potOut::potential(double x) const
{
	//const double angularpart=((Parameters::hbar*Parameters::hbar)*Parameters::angularMomentum*(Parameters::angularMomentum+1))/(2*Parameters::mn*(x*x));
	return /*angularpart+*/m_potKS;

}
InitialPot* potOut::clone() const
{
    return new potOut(*this); //return a derived class object through a base class pointer
}












