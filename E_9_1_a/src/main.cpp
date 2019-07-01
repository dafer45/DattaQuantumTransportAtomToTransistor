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
#include "TBTK/Plotter2.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/PropertyExtractor/Greens.h"
#include "TBTK/Range.h"
#include "TBTK/Solver/Diagonalizer.h"
#include "TBTK/Solver/Greens.h"
#include "TBTK/Streams.h"
#include "TBTK/UnitHandler.h"

#include <complex>

using namespace std;
using namespace TBTK;
//using namespace Plot;

complex<double> i(0, 1);

class CallbackU : public HoppingAmplitude::AmplitudeCallback{
public:
	complex<double> getHoppingAmplitude(
		const Index &to,
		const Index &from
	) const{
		int x = from[0];

		switch(potentialType){
		case 0:
			return 0;
		case 1:
			if(x > 22 && x < 27)
				return U;
			else
				return 0;
		case 2:
			if((x > 14 && x < 19) || (x > 26 && x < 31))
				return U;
			else
				return 0;
		default:
			Streams::out << "Error: invalid potential type.\n";
			exit(1);
		}
	}

	void setPotentialType(unsigned int potentialType){
		this->potentialType = potentialType;
	}

	void setU(double U){
		this->U = U;
	}
private:
	unsigned int potentialType;
	double U;
} callbackU;

int main(int argc, char **argv){
	//Set the natural units. Argument order: (charge, count, energy,
	//length, temperature, time).
	UnitHandler::setScales({"1 C", "1 pcs", "1 eV", "1 Ao", "1 K", "1 s"});

	//Paramters.
	const double LOWER_BOUND = -0.2;		//eV
	const double UPPER_BOUND = 0.8;			//eV
	const int RESOLUTION = 1000;
	double U = 0.4;					//eV
	double a = 3;					//Angstrom
	double m_c = 0.25*UnitHandler::getM_eN();
	double hbar = UnitHandler::getHbarN();
	double t = hbar*hbar/(2*m_c*a*a);		//eV
	int sizeX = 50;
	Array<double> UB({(unsigned int)sizeX}, 0);
	callbackU.setU(U);

	//Setup the model.
	Model model;
	for(int x = 0; x < sizeX; x++){
		model << HoppingAmplitude(2*t, {x}, {x});
		model << HoppingAmplitude(callbackU, {x}, {x});

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

	for(unsigned int m = 0; m < 3; m++){
		callbackU.setPotentialType(m);

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
		}

		//Calculate the full Green's function including self-energies.
		solver1.setGreensFunction(greensFunction);
		Property::GreensFunction fullGreensFunction
			= solver1.calculateInteractingGreensFunction(
				selfEnergy0 + selfEnergy1
			);

		//Calculate the transmission rate.
		solver1.setGreensFunction(fullGreensFunction);
		Property::TransmissionRate transmissionRate
			= solver1.calculateTransmissionRate(
				selfEnergy0,
				selfEnergy1
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
		Plotter2 plotter;
		plotter.setLabelX("E (eV)");
		plotter.setLabelY("Transmission rate");
		plotter.plot(energy, transmissionRateData);
		plotter.save("figures/TransmissionRate_" + to_string(m) + ".png");
	}

	return 0;
}
