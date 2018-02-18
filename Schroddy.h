/*
 * Schroddy.h
 *
 *  Created on: 17 feb 2018
 *      Author: brando
 */
#pragma once
#include <vector>

class InitialPot;
class Schroddy;

/*==================================================================
 * Pure virtual class for eigenvalues
 *================================================================*/

class Eigenvalues
{
public:
    virtual double eigenvalue() const = 0;
    virtual ~Eigenvalues();

    virtual Eigenvalues* clone() const = 0;

private:

};

/*===================================================================
 * HO eigenvalues generator class
 *=================================================================*/

class HarmonicEigenvalues: public Eigenvalues
{
public:
    HarmonicEigenvalues(unsigned int n, int l);
    double eigenvalue() const override;
    //setEnergyLevel(unsigned int n);
    //setAngularMomentum(int l);

    Eigenvalues* clone() const override;

private:
    HarmonicEigenvalues();

    //if you want to be able to change quantum numbers within same object, implement set methods and remove const keyword
    const unsigned int m_n;
    const int m_l;
};

/*=======================================================================
 * Shroedinger class + solver
 *=====================================================================*/

class Schroddy
{
public:
    Schroddy(const InitialPot& pot, unsigned int nQuant);
    double solveSchroddyByRK(double x0, double x1, double psi0, double psiPrime0) const; //psiPrime0 is boundary condition on first derivative of eigenfunction

    ~Schroddy();

private:
    Schroddy();
    InitialPot* m_pot;
    unsigned int m_nQuant;
};



