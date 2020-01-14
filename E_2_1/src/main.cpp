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
 *  @brief Exercise 2.1
 *
 *  Solution to exercise 2.1 in the book "Quantum Transport: Atom to
 *  Transistor, S. Datta (2005)".
 *
 *  @author Kristofer Björnson
 */

#include "TBTK/Model.h"
#include "TBTK/Property/WaveFunctions.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/Solver/Diagonalizer.h"
#include "TBTK/Streams.h"
#include "TBTK/TBTK.h"
#include "TBTK/UnitHandler.h"
#include "TBTK/Visualization/MatPlotLib/Plotter.h"

#include <complex>

using namespace std;
using namespace TBTK;
using namespace Visualization::MatPlotLib;

int main(int argc, char **argv){
	//Initialize TBTK.
	Initialize();

	//Set the natural units. Argument order: (angle, charge, count, energy,
	//length, temperature, time).
	UnitHandler::setScales(
		{"1 rad", "1 C", "1 pcs", "1 eV", "1 Ao", "1 K", "1 s"}
	);

	//Parameters.
	double a = 1;
	double hbar = UnitHandler::getConstantInNaturalUnits("hbar");
	double m_e = UnitHandler::getConstantInNaturalUnits("m_e");
	double t_0 = pow(hbar/a, 2)/(2*m_e);
	const int SIZE_X = 100;

	//Set this flag to false for to run problem a and to true to run
	//problem b.
	const bool PERIODIC_BOUNDARY_CONDITIONS = false;

	//Setup the model.
	Model model;
	for(int x = 0; x < SIZE_X; x++){
		model << HoppingAmplitude(2*t_0, {x}, {x});

		if(x+1 < SIZE_X)
			model << HoppingAmplitude(-t_0, {x+1}, {x}) + HC;
		if(x+1 == SIZE_X && PERIODIC_BOUNDARY_CONDITIONS)
			model << HoppingAmplitude(-t_0, {0}, {x}) + HC;
	}
	model.construct();

	//Setup and run solver.
	Solver::Diagonalizer solver;
	solver.setModel(model);
	solver.run();

	//Extract eigenvalues.
	PropertyExtractor::Diagonalizer propertyExtractor(solver);
	Property::EigenValues eigenValues = propertyExtractor.getEigenValues();

	if(!PERIODIC_BOUNDARY_CONDITIONS){
		//This is for problem a.

		//Calculate the analytical eigenvalues.
		Array<double> analyticalEigenValues({SIZE_X});
		double L = (SIZE_X+1)*a;
		for(unsigned int alpha = 1; alpha <= SIZE_X; alpha++){
			analyticalEigenValues[{alpha-1}]
				= pow(hbar*M_PI*alpha/L, 2)/(2*m_e);
		}

		//Plot the results.
		Plotter plotter;
		plotter.setLabelX("Eigenvalue number");
		plotter.setLabelY("Energy");
		plotter.plot(eigenValues.getData());
		plotter.plot(analyticalEigenValues);
		plotter.save("figures/EigenValues.png");

		//Extract the wave functions for states 0 and 49. Note that zero based
		//indexing is used, in contrast to the one based indexing in the book.
		int states[2] = {0, 49};
		Property::WaveFunctions waveFunctions
			= propertyExtractor.calculateWaveFunctions(
				{{IDX_ALL}},
				{states[0], states[1]}
			);

		//Calculate the probability distribution.
		Array<double> probabilityDistribution({2, SIZE_X});
		for(unsigned int x = 0; x < SIZE_X; x++){
			probabilityDistribution[{0, x}]
				= pow(abs(waveFunctions({(int)x}, states[0])), 2);
			probabilityDistribution[{1, x}]
				= pow(abs(waveFunctions({(int)x}, states[1])), 2);
		}

		//Plot the probability distributions.
		plotter.clear();
		plotter.setLabelX("Lattice site number");
		plotter.setLabelY("Probability density");
		plotter.plot(probabilityDistribution.getSlice({0, IDX_ALL}));
		plotter.plot(probabilityDistribution.getSlice({1, IDX_ALL}));
		plotter.save("figures/ProbabilityDistribution.png");

	}
	else{
		//This is for problem b.

		//Plot the results.
		Plotter plotter;
		plotter.setLabelX("Eigenvalue number");
		plotter.setLabelY("Energy");
		plotter.plot(eigenValues);
		plotter.save("figures/EigenValuesPeriodicBoundaryConditions.png");
	}

	return 0;
}
