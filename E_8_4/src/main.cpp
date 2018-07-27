/* Copyright 2018 Kristofer Björnson
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
 *  @brief Exercise 8.4
 *
 *  Solution to exercise 8.4 in the book "Quantum Transport: Atom to
 *  Transistor, S. Datta (2005)".
 *
 *  @author Kristofer Björnson
 */

#include "TBTK/Model.h"
#include "TBTK/Plotter.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/PropertyExtractor/Greens.h"
#include "TBTK/Solver/Diagonalizer.h"
#include "TBTK/Solver/Greens.h"
#include "TBTK/Streams.h"
#include "TBTK/UnitHandler.h"

#include <complex>

using namespace std;
using namespace TBTK;
using namespace Plot;

complex<double> i(0, 1);

int main(int argc, char **argv){
	//Set the natural units. Argument order: (charge, count, energy,
	//length, temperature, time).
	UnitHandler::setScales({"1 C", "1 pcs", "1 eV", "1 Ao", "1 K", "1 s"});

	//Paramters.
	const double LOWER_BOUND = -1;	//eV
	const double UPPER_BOUND = 1;	//eV
	const int RESOLUTION = 1000;
	double epsilon = -0.25;		//eV
	double epsilon_1 = 0.25;	//eV
	double t = 0.5;			//eV

	////////////////////
	// Exercise 8.4.a //
	////////////////////

	//Setup the model.
	Model modelA;
	modelA << HoppingAmplitude(epsilon,	{0},	{0});
	modelA << HoppingAmplitude(epsilon_1,	{1},	{1});
	modelA << HoppingAmplitude(t,		{0},	{1}) + HC;
	modelA.construct();

	//Setup solver for calculating the Green's function.
	Solver::Diagonalizer solverA;
	solverA.setModel(modelA);
	solverA.run();

	//Setup the PropertyExtractor for calculating the Green's function.
	PropertyExtractor::Diagonalizer propertyExtractorA(solverA);
	propertyExtractorA.setEnergyWindow(
		LOWER_BOUND,
		UPPER_BOUND,
		RESOLUTION
	);
	propertyExtractorA.setEnergyInfinitesimal(0.01);

	//Calculate the Green's function.
	Property::GreensFunction greensFunctionA0
		= propertyExtractorA.calculateGreensFunction(
			{{Index({IDX_ALL}), Index({IDX_ALL})}}
		);

	//Setup the solver for calculating the LDOS.
	Solver::Greens solverA1;
	solverA1.setModel(modelA);
	solverA1.setGreensFunction(greensFunctionA0);

	//Setup PropertyExtractor and calculate the LDOS.
	PropertyExtractor::Greens propertyExtractorA1(solverA1);
	Property::LDOS ldosA = propertyExtractorA1.calculateLDOS({{IDX_ALL}});

	//Prepare the LDOS for plotting.
	Array<double> l({2, RESOLUTION});
	for(unsigned int n = 0; n < RESOLUTION; n++){
		l[{0, n}] = ldosA({0}, n);
		l[{1, n}] = ldosA({1}, n);
	}

	//Plot the LDOS.
	Plotter plotter;
	plotter.plot(l.getSlice({0, IDX_ALL}));
	plotter.save("figures/LDOS_A0.png");
	plotter.plot(l.getSlice({1, IDX_ALL}));
	plotter.save("figures/LDOS_A1.png");

	////////////////////
	// Exercise 8.4.b //
	////////////////////
	for(unsigned int n = 0; n < 2; n++){
		//Setup the model.
		Model modelB;
		if(n == 0)
			modelB << HoppingAmplitude(epsilon, {0}, {0});
		else
			modelB << HoppingAmplitude(epsilon_1, {0}, {0});
		modelB.construct();

		//Setup solver for calculating the non-interacting Green's
		//function.
		Solver::Diagonalizer solverB;
		solverB.setModel(modelB);
		solverB.run();

		//Setup the PropertyExtractor for calculating the
		//non-interacting Green's function.
		PropertyExtractor::Diagonalizer propertyExtractorB(solverB);
		propertyExtractorB.setEnergyWindow(
			LOWER_BOUND,
			UPPER_BOUND,
			RESOLUTION
		);
		propertyExtractorB.setEnergyInfinitesimal(0.01);

		//Calculate the non-interacting Green's function.
		Property::GreensFunction greensFunctionB0
			= propertyExtractorB.calculateGreensFunction(
				{{Index({0}), Index({0})}},
				Property::GreensFunction::Type::Retarded
			);

		//Setup the self-energy.
		IndexTree memoryLayoutSelfEnergy;
		memoryLayoutSelfEnergy.add({Index({0}), Index({0})});
		memoryLayoutSelfEnergy.generateLinearMap();
		Property::SelfEnergy selfEnergy(
			memoryLayoutSelfEnergy,
			LOWER_BOUND,
			UPPER_BOUND,
			RESOLUTION
		);
		for(unsigned int e = 0; e < RESOLUTION; e++){
			double E = LOWER_BOUND
				+ e*(UPPER_BOUND - LOWER_BOUND)/RESOLUTION;

			if(n == 0){
				selfEnergy({Index({0}), Index({0})}, e)
					= t*t/(E - epsilon_1 + i*0.01);
			}
			else{
				selfEnergy({Index({0}), Index({0})}, e)
					= t*t/(E - epsilon + i*0.01);
			}
		}

		//Setup the solver for calculating the interacting Green's
		//function and the LDOS.
		Solver::Greens solverC;
		solverC.setModel(modelB);
		solverC.setGreensFunction(greensFunctionB0);

		//Calculate the interacting Green's function and reconfigure
		//the solver to use it.
		Property::GreensFunction greensFunctionB
			= solverC.calculateInteractingGreensFunction(
				selfEnergy
			);
		solverC.setGreensFunction(greensFunctionB);

		//Setup the PropertyExtractor and calculate the LDOS.
		PropertyExtractor::Greens propertyExtractorC(solverC);
		Property::LDOS ldosB
			= propertyExtractorC.calculateLDOS({{IDX_ALL}});

		//Prepare the LDOS for plotting.
		for(unsigned int e = 0; e < RESOLUTION; e++)
			l[{n, e}] = ldosB({0}, e);
	}

	//Plot the LDOS.
	plotter.plot(l.getSlice({0, IDX_ALL}));
	plotter.save("figures/LDOS_B0.png");
	plotter.plot(l.getSlice({1, IDX_ALL}));
	plotter.save("figures/LDOS_B1.png");

	return 0;
}
