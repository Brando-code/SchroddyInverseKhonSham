/*
 * main.cpp
 *
 *  Created on: 18 feb 2018
 *      Author: brando
 */
#include <iostream>
#include <fstream>
#include <vector>
#include "Includes.h"

int main(int argc, const char * argv[]) {

    unsigned int n=5;
    int l_mom=0;

	HOPot pot (Parameters::mn, l_mom);
	Schroddy sol (pot, n);
	//std::vector<double> solution = sol.solveSchroddyByRK(Parameters::x_in, Parameters::x_fin, Parameters::psi0, Parameters::psiPrime0);
	std::vector<double> solutiuon = sol.solveSchroddyByNV(Parameters::psiLeft, Parameters::psiRight);

	//std::cout << solution << std::endl;


	std::cout << "Program executed successfully." << std::endl;
    return 0;
}



