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
 *  @brief Exercise 8.3
 *
 *  Solution to exercise 8.3 in the book "Quantum Transport: Atom to
 *  Transistor, S. Datta (2005)".
 *
 *  @author Kristofer Björnson
 */

#include "TBTK/Model.h"
#include "TBTK/Plotter.h"
#include "TBTK/Property/Density.h"
#include "TBTK/Property/SelfEnergy.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/PropertyExtractor/Greens.h"
#include "TBTK/Range.h"
#include "TBTK/Solver/Diagonalizer.h"
#include "TBTK/Solver/Greens.h"
#include "TBTK/Solver/LinearEquationSolver.h"
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

	//Parameters.
	double hbar = UnitHandler::getHbarN();
	double m_e = UnitHandler::getM_eN();
	double m = 0.25*m_e;
	double a = 2;	//Ångström.
	double t = hbar*hbar/(2*m*a*a);
	double mu = 0.25;	//eV;
	const unsigned int SIZE = 50;

	////////////////////
	// Exercise 8.3.a //
	////////////////////

	//Setup the model with periodic boundary conditions.
	Model modelA;
	Range U(-0.05, 0.05, SIZE);
	for(unsigned int n = 0; n < SIZE; n++){
		modelA << HoppingAmplitude(U[{n}] + 2*t, {n}, {n});
		modelA << HoppingAmplitude(-t, {(n+1)%SIZE}, {n}) + HC;
	}
	modelA.construct();
	modelA.setChemicalPotential(mu);
	modelA.setTemperature(300);

	//Setup and run the solver.
	Solver::Diagonalizer solverA;
	solverA.setModel(modelA);
	solverA.run();

	//Extract the density.
	PropertyExtractor::Diagonalizer propertyExtractorA(solverA);
	Property::Density densityA
		= propertyExtractorA.calculateDensity({{IDX_ALL}});

	//Prepare the density for plotting.
	Array<double> d({SIZE});
	for(unsigned int n = 0; n < SIZE; n++)
		d[{n}] = densityA(n)/a;

	//Plot the density.
	Plotter plotter;
	plotter.setLabelX("x");
	plotter.setLabelY("Density (Ao^-1)");
	plotter.setHold(true);
	plotter.plot(d);

	////////////////////
	// Exercise 8.3.b //
	////////////////////

	//Setup the model without periodic boundary conditions.
	Model modelB;
	for(unsigned int n = 0; n < SIZE; n++){
		modelB << HoppingAmplitude(U[{n}] + 2*t, {n}, {n});
		if(n+1 < SIZE)
			modelB << HoppingAmplitude(-t, {n+1}, {n}) + HC;
	}
	modelB.construct();
	modelB.setChemicalPotential(mu);
	modelB.setTemperature(300);

	//Setup and run the solver.
	Solver::Diagonalizer solverB;
	solverB.setModel(modelB);
	solverB.run();

	//Parameters for the Green's function and self-ebergy.
	const double LOWER_BOUND = -0.1;
	const double UPPER_BOUND = 0.4;
	const int RESOLUTION = 250;

	//Setup the PropertyExtractor.
	PropertyExtractor::Diagonalizer propertyExtractorB(solverB);
	propertyExtractorB.setEnergyWindow(
		LOWER_BOUND,
		UPPER_BOUND,
		RESOLUTION
	);
	propertyExtractorB.setEnergyInfinitesimal(1e-12);

	//Extract the non-interacting Green's function.
	Property::GreensFunction greensFunction0
		= propertyExtractorB.calculateGreensFunction(
			{{Index({IDX_ALL}), Index({IDX_ALL})}},
			Property::GreensFunction::Type::Retarded
		);

	//Setup the self-energy.
	IndexTree memoryLayoutSelfEnergy;
	for(unsigned int x = 0; x < SIZE; x++)
		for(unsigned int xp = 0; xp < SIZE; xp++)
			memoryLayoutSelfEnergy.add({Index({x}), Index({xp})});
	memoryLayoutSelfEnergy.generateLinearMap();
	Property::SelfEnergy selfEnergy(
		memoryLayoutSelfEnergy,
		LOWER_BOUND,
		UPPER_BOUND,
		RESOLUTION
	);
	for(unsigned int x = 0; x < SIZE; x++){
		for(unsigned int xp = 0; xp < SIZE; xp++){
			for(unsigned int e = 0; e < RESOLUTION; e++){
				//Calculate the energy.
				double E = LOWER_BOUND + e*(
					UPPER_BOUND - LOWER_BOUND
				)/RESOLUTION;

				//Invert the relation
				//E = E_c + 2t_0(1 - cos(ka)).
				complex<double> k = acos(
					1. - (E - U[{x}] + i*1e-12)/(2*t)
				)/a;

				//Set the self-energy values.
				if(x == xp && (x == 0 || x == SIZE-1)){
					selfEnergy(
						{Index({x}), Index({xp})},
						e
					) = -t*exp(i*k*a);
				}
				else{
					selfEnergy(
						{Index({x}), Index({xp})},
						e
					) = 0;
				}
			}
		}
	}

	//Setup the Green's solver.
	Solver::Greens solverC;
	solverC.setModel(modelB);
	solverC.setGreensFunction(greensFunction0);

	//Extract the interacting Green's function and reset the Green's solver
	//to use it instead of the non-interacting Green's function.
	Property::GreensFunction greensFunction
		= solverC.calculateInteractingGreensFunction(selfEnergy);
	solverC.setGreensFunction(greensFunction);

	//Setup the PropertyExtractor and extract the density.
	PropertyExtractor::Greens propertyExtractorC(solverC);
	Property::Density densityB = propertyExtractorC.calculateDensity(
		{{IDX_ALL}}
	);

	//Prepare the density for plotting.
	for(unsigned int n = 0; n < SIZE; n++)
		d[{n}] = densityB(n)/a;

	//Plot the density.
	plotter.plot(
		d,
		Decoration({0, 0, 0}, Decoration::LineStyle::Point)
	);

	//Save the plot.
	plotter.save("figures/Density.png");

	//Calculate the LDOS at the two end points.
	Property::LDOS ldos = propertyExtractorC.calculateLDOS(
		{{0}, {SIZE-1}}
	);

	//Prepare the LDOS for plotting.
	Array<double> l({2, RESOLUTION});
	for(unsigned int e = 0; e < RESOLUTION; e++){
		l[{0, e}] = ldos({0}, e);
		l[{1, e}] = ldos({SIZE-1}, e);
	}

	//Plot the LDOS.
	plotter.setHold(false);
	plotter.setLabelX("Energy");
	plotter.setLabelY("LDOS");
	plotter.plot(l.getSlice({0, IDX_ALL}));
	plotter.save("figures/LDOSLeft.png");
	plotter.plot(l.getSlice({1, IDX_ALL}));
	plotter.save("figures/LDOSRight.png");

	return 0;
}
