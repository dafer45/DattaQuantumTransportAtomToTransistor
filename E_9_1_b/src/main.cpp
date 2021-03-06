/* Copyright 2019 Kristofer Björnson
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/** @package TBTKQuantumTransportAtomToTransistor
 *  @file main.cpp
 *  @brief Exercise 9,1
 *
 *  Solution to exercise 9.1 in the book "Quantum Transport: Atom to
 *  Transistor, S. Datta (2005)".
 *
 *  @author Kristofer Björnson
 */

#include "TBTK/Array.h"
#include "TBTK/Model.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/PropertyExtractor/Greens.h"
#include "TBTK/Range.h"
#include "TBTK/Solver/Diagonalizer.h"
#include "TBTK/Solver/Greens.h"
#include "TBTK/Streams.h"
#include "TBTK/TBTK.h"
#include "TBTK/UnitHandler.h"
#include "TBTK/Visualization/MatPlotLib/Plotter.h"

#include <complex>

using namespace std;
using namespace TBTK;
using namespace Visualization::MatPlotLib;

complex<double> i(0, 1);

int main(int argc, char **argv){
	//Initialize TBTK.
	Initialize();

	//Set the natural units. Argument order: (angle, charge, count, energy,
	//length, temperature, time).
	UnitHandler::setScales(
		{"1 rad", "1 C", "1 pcs", "1 eV", "1 Ao", "1 K", "1 s"}
	);

	//Paramters.
	const double LOWER_BOUND = -0.2;		//eV
	const double UPPER_BOUND = 0.8;			//eV
	const int RESOLUTION = 1000;
	double U = 0.4;					//eV
	double a = 3;					//Angstrom
	double m_c = 0.25*UnitHandler::getConstantInNaturalUnits("m_e");
	double hbar = UnitHandler::getConstantInNaturalUnits("hbar");
	double t = hbar*hbar/(2*m_c*a*a);		//eV
	int sizeX = 46;
	Array<double> UB({(unsigned int)sizeX}, 0);

	//Setup the potential barrier.
	for(unsigned int x = 0; x < (unsigned int)sizeX; x++)
		if((x > 14 && x < 19) || (x > 26 && x < 31))
			UB[{x}] = U;

	//Setup the model.
	Model model;
	for(int x = 0; x < sizeX; x++){
		model << HoppingAmplitude(2*t, {x}, {x});
		model << HoppingAmplitude(UB[{(unsigned int)x}], {x}, {x});

		if(x + 1 < sizeX)
			model << HoppingAmplitude(-t, {x}, {x+1}) + HC;
	}
	model.construct();

	//Setup the solver that is used to calculate the bare Green's function.
	Solver::Diagonalizer solver;
	solver.setModel(model);

	//Setup the PropertyExtractor for calculating the Green's function.
	PropertyExtractor::Diagonalizer propertyExtractor(solver);
	propertyExtractor.setEnergyWindow(
		LOWER_BOUND,
		UPPER_BOUND,
		RESOLUTION
	);
	propertyExtractor.setEnergyInfinitesimal(1e-10);

	//Setup the solver that is used to calculate the full Green's function
	//and the transmission rate.
	Solver::Greens solver1;
	solver1.setModel(model);

	//Diagonalize the Hamiltonian.
	solver.run();

	//Calculate the Green's function.
	Property::GreensFunction greensFunction
		= propertyExtractor.calculateGreensFunction(
			{{Index({IDX_ALL}), Index({IDX_ALL})}}
		);

	//Setup the self-energies.
	IndexTree selfEnergyIndices;
	for(int x = 0; x < sizeX; x++)
		for(int xp = 0; xp < sizeX; xp++)
			selfEnergyIndices.add({Index({x}), Index({xp})});
	selfEnergyIndices.generateLinearMap();
	Property::SelfEnergy selfEnergy0(
		selfEnergyIndices,
		LOWER_BOUND,
		UPPER_BOUND,
		RESOLUTION
	);
	Property::SelfEnergy selfEnergy1(
		selfEnergyIndices,
		LOWER_BOUND,
		UPPER_BOUND,
		RESOLUTION
	);
	Property::SelfEnergy selfEnergy2(
		selfEnergyIndices,
		LOWER_BOUND,
		UPPER_BOUND,
		RESOLUTION
	);
	Range energies(LOWER_BOUND, UPPER_BOUND, RESOLUTION);
	for(unsigned int n = 0; n < RESOLUTION; n++){
		complex<double> ka = acos(
			complex<double>(
				1 - (energies[n] - UB[{0}])/(2*t),
				0
			)
		);
		selfEnergy0({Index({0}), Index({0})}, n) = -t*exp(i*ka);

		ka = acos(
			complex<double>(
				1 - (
					energies[n]
					- UB[{(unsigned int)(sizeX-1)}]
				)/(2*t),
				0
			)
		);
		selfEnergy1({Index({sizeX-1}), Index({sizeX-1})}, n) = -t*exp(i*ka);

		selfEnergy2({Index({22}), Index({22})}, n) = -0.25*i;
	}

	//Calculate the full Green's function including self-energies.
	solver1.setGreensFunction(greensFunction);
	Property::GreensFunction fullGreensFunction
		= solver1.calculateInteractingGreensFunction(
			selfEnergy0 + selfEnergy1 + selfEnergy2
		);

	//Calculate the transmission rates between the different terminals.
	solver1.setGreensFunction(fullGreensFunction);
	Property::TransmissionRate transmissionRate0to1
		= solver1.calculateTransmissionRate(
			selfEnergy0,
			selfEnergy1
		);
	Property::TransmissionRate transmissionRate0to2
		= solver1.calculateTransmissionRate(
			selfEnergy0,
			selfEnergy2
		);
	Property::TransmissionRate transmissionRate1to2
		= solver1.calculateTransmissionRate(
			selfEnergy1,
			selfEnergy2
		);

	//Calculate the total transmission rate using Eq. (9.5.7).
	Property::TransmissionRate transmissionRate = transmissionRate0to1
		+ transmissionRate0to2*transmissionRate1to2/(
			transmissionRate0to2 + transmissionRate1to2
		);

	//Plot the results.
	vector<double> energy;
	vector<double> transmissionRateData;
	for(
		unsigned int c = 0;
		c < transmissionRate.getResolution();
		c++
	){
		energy.push_back(energies[c]);
		transmissionRateData.push_back(transmissionRate(c));
	}
	Plotter plotter;
	plotter.setLabelX("Energy (eV)");
	plotter.setLabelY("Transmission rate");
	plotter.plot(energy, transmissionRateData);
	plotter.save("figures/TransmissionRate.png");

	return 0;
}
