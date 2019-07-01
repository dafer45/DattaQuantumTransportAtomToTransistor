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
#include "TBTK/Functions.h"
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
	double U = 0.4;						//eV
	double a = 3;						//Angstrom
	double m_c = 0.25*UnitHandler::getM_eN();
	double hbar = UnitHandler::getHbarN();
	double q = UnitHandler::getEN();
	double t = hbar*hbar/(2*m_c*a*a);			//eV
	double E_f = 0.1;					//eV
	double temperature = 0.025/UnitHandler::getK_BN();	//K
	const double LOWER_BOUND = -0.2;			//eV
	const double UPPER_BOUND = 0.8;				//eV
	const int RESOLUTION = 101;
	int sizeX = 50;
	Array<double> U1({(unsigned int)sizeX}, 0);
	callbackU.setU(U);

	for(unsigned int m = 0; m < 3; m++){
		callbackU.setPotentialType(m);

		vector<double> voltages;
		vector<double> currents;
		Range V(0, 0.5, 26);
		for(unsigned int v = 0; v < 26; v++){
			voltages.push_back(V[v]);
			double mu0 = E_f + V[v]/2.;
			double mu1 = E_f - V[v]/2.;

			//Setup bias potential.
			switch(m){
			case 0:
			{
				Range slope(0.5, -0.5, 20);
				for(unsigned int x = 0; x < (unsigned int)sizeX; x++){
					if(x < 15)
						U1[{x}] = V[v]*0.5;
					else if(x < 35)
						U1[{x}] = V[v]*slope[x - 15];
					else
						U1[{x}] = -V[v]*0.5;
				}

				break;
			}
			case 1:
			{
				Range slope(0.5, -0.5, 4);
				for(unsigned int x = 0; x < (unsigned int)sizeX; x++){
					if(x < 23)
						U1[{x}] = V[v]*0.5;
					else if(x < 27)
						U1[{x}] = V[v]*slope[x - 23];
					else
						U1[{x}] = -V[v]*0.5;
				}

				break;
			}
			case 2:
			{
				Range slope(0.5, -0.5, 16);
				for(unsigned int x = 0; x < (unsigned int)sizeX; x++){
					if(x < 15)
						U1[{x}] = V[v]*0.5;
					else if(x < 31)
						U1[{x}] = V[v]*slope[x - 15];
					else
						U1[{x}] = -V[v]*0.5;
				}

				break;
			}
			default:
				Streams::out << "Error: invalid potential type.",
				exit(1);
			}

			//Setup the model.
			Model model;
			for(int x = 0; x < sizeX; x++){
				model << HoppingAmplitude(2*t, {x}, {x});
				model << HoppingAmplitude(callbackU, {x}, {x});
				model << HoppingAmplitude(
					U1[{(unsigned int)x}],
					{x},
					{x}
				);

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
						1 - (energies[n] - U1[{0}])/(2*t),
						0
					)
				);
				selfEnergy0({Index({0}), Index({0})}, n) = -t*exp(i*ka);

				ka = acos(
					complex<double>(
						1 - (
							energies[n]
							- U1[{(unsigned int)(sizeX-1)}]
						)/(2*t),
						0
					)
				);
				selfEnergy1({Index({sizeX-1}), Index({sizeX-1})}, n) = -t*exp(i*ka);

				selfEnergy2({Index({sizeX/2-1}), Index({sizeX/2-1})}, n) = -i*0.00025;
			}

			//Calculate the full Green's function including self-energies.
			solver1.setGreensFunction(greensFunction);
			Property::GreensFunction fullGreensFunction
				= solver1.calculateInteractingGreensFunction(
					selfEnergy0 + selfEnergy1 + selfEnergy2
				);

			//Calculate the transmission rate.
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

			//Calculate the total transmissionRate. The fraction
			//can lead to 0/0 division, which result in NaN.
			//Replace such Nan values by 0.
			Property::TransmissionRate transmissionRate
				= transmissionRate0to2*transmissionRate1to2/(
					transmissionRate0to2
					+ transmissionRate1to2
				);
			transmissionRate.replaceValues(NAN, 0);
			transmissionRate += transmissionRate0to1;

			//Calculate the current.
			double current = 0;
			double dE = transmissionRate.getEnergy(1)
				- transmissionRate.getEnergy(0);
			for(
				unsigned int n = 0;
				n < transmissionRate.getResolution();
				n++
			){
				current += dE*transmissionRate(n)*(
					Functions::fermiDiracDistribution(
						transmissionRate.getEnergy(n),
						mu0,
						temperature
					)
					- Functions::fermiDiracDistribution(
						transmissionRate.getEnergy(n),
						mu1,
						temperature
					)
				);
			}
			current *= q*q/(2*M_PI*hbar);
			currents.push_back(current);

			//Plot the results.
			Plotter2 plotter;
			plotter.setLabelX("Voltage (V)");
			plotter.setLabelY("Current (A)");
			plotter.plot(voltages, currents);
			plotter.save("figures/Current_" + to_string(m) + ".png");
		}
	}

	return 0;
}
