/*
 * Potentials.h
 *
 *  Created on: 17 feb 2018
 *      Author: brando
 */
//class hierarchy interface goes in header
#pragma once

/*====================================================================
 * Pure virtual base class for in/out potentials
 *===================================================================*/

class InitialPot
{
public:
    virtual double potential(double x) const = 0;
    virtual ~InitialPot();

    virtual InitialPot* clone() const = 0; //virtual constructor

private:

};

/*=====================================================================
 * Woods-Saxon potential class
 *===================================================================*/

class WSaxPot: public InitialPot //derived class
{

public:
    WSaxPot(double V0, double Rn, double a0);
    double potential(double x) const override;

    InitialPot* clone() const override;

private:
    WSaxPot(); //make defaul ctor private as we must initialize parameters
    const double m_V0, m_Rn, m_a0;

};


/*=====================================================================
 * HO potential class
 *===================================================================*/

class HOPot: public InitialPot //derived class
{
public:

    HOPot(double m, int anglmomentum);
    double potential(double x) const override;

    //double getOmega() const;
    InitialPot* clone() const override;

private:
    HOPot(); //same as before
    const double m_m, m_anglmomentum;
};

/*=====================================================================
 * Spin-Orbit potential class
 *===================================================================*/
/*
class SOPot: public InitialPot //derived class
{
public:

    SOPot(double m, double omega);
    double potential(double x) const override;

    //double getOmega() const;
    InitialPot* clone() const override;

private:
    SOPot(); //same as before
    const double m_m, m_omega;
};
*/

/*======================================================================
 * Kohn-Sham potential class
 *====================================================================*/

class potOut: public InitialPot
{
public:
	potOut(double potKS);
	double potential (double x) const override;

	InitialPot* clone() const override;

private:
	double m_potKS;
};












