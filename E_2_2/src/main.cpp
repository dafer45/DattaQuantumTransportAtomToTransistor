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
 *  @brief Exercise 2.2
 *
 *  Solution to exercise 2.2 in the book "Quantum Transport: Atom to
 *  Transistor, S. Datta (2005)".
 *
 *  @author Kristofer Björnson
 */

#include "TBTK/Model.h"
#include "TBTK/Property/WaveFunctions.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/Range.h"
#include "TBTK/Solver/Diagonalizer.h"
#include "TBTK/Streams.h"
#include "TBTK/UnitHandler.h"
#include "TBTK/Visualization/MatPlotLib/Plotter.h"

#include <complex>

using namespace std;
using namespace TBTK;
using namespace Visualization::MatPlotLib;

int main(int argc, char **argv){
	//Initialize TBTK.
	Initialize();

	//Set the natural units. Argument order: (charge, count, energy,
	//length, temperature, time).
	UnitHandler::setScales(
		{"1 rad", "1 C", "1 pcs", "1 eV", "1 Ao", "1 K", "1 s"}
	);

	//Parameters.
	double a = 0.05;
	double hbar = UnitHandler::getConstantInNaturalUnits("hbar");
	double m_e = UnitHandler::getConstantInNaturalUnits("m_e");
	double e = UnitHandler::getConstantInNaturalUnits("e");
	double epsilon_0 = UnitHandler::getConstantInNaturalUnits("epsilon_0");
	double A = hbar*hbar/(2*m_e*a*a);
	double B = e*e/(4*M_PI*epsilon_0);
	double C = hbar*hbar/(2*m_e);
	double a_0 = 4*M_PI*epsilon_0*pow(hbar/e, 2)/m_e;
	const int SIZE_R = 100;
	double l = 0;

	//Setup the model.
	Model model;
	Range r(a, a*SIZE_R, SIZE_R);
	for(unsigned int n = 0; n < r.getResolution(); n++){
		model << HoppingAmplitude(
			2*A - B/r[n] + l*(l+1)*C/(r[n]*r[n]),
			{(int)n},
			{(int)n}
		);

		if(n+1 < r.getResolution()){
			model << HoppingAmplitude(
				-A,
				{(int)(n+1)},
				{(int)n}
			) + HC;
		}
	}
	model.construct();

	//Setup and run solver.
	Solver::Diagonalizer solver;
	solver.setModel(model);
	solver.run();

	//Extract the wave functions for states 1s and 2s. Note that zero based
	//indexing is used.
	int states[2] = {0, 1};
	PropertyExtractor::Diagonalizer propertyExtractor(solver);
	Property::WaveFunctions waveFunctions
		= propertyExtractor.calculateWaveFunctions(
			{{IDX_ALL}},
			{states[0], states[1]}
		);

	//Calculate the probability distribution.
	Array<double> probabilityDistribution({4, SIZE_R});
	for(unsigned int n = 0; n < SIZE_R; n++){
		//Numerical.
		probabilityDistribution[{0, n}]
			= pow(abs(waveFunctions({(int)n}, states[0])), 2);
		probabilityDistribution[{1, n}]
			= pow(abs(waveFunctions({(int)n}, states[1])), 2);

		//Analytical.
		probabilityDistribution[{2, n}]
			= (4*a*r[n]*r[n]/pow(a_0, 3))*exp(-2*r[n]/(a_0));
		probabilityDistribution[{3, n}]
			= (
				a*r[n]*r[n]/(8*pow(a_0, 3))
			)*pow(2 - r[n]/a_0, 2)*exp(-2*r[n]/(2*a_0));
	}

	//Plot the probability distributions.
	Plotter plotter;
	plotter.setLabelX("r");
	plotter.setLabelY("Probability density");

	plotter.plot(probabilityDistribution.getSlice({0, IDX_ALL}));
	plotter.plot(
		probabilityDistribution.getSlice({2, IDX_ALL}),
		{{"linestyle", " "}, {"marker", "o"}, {"markersize", "5"}}
	);
	plotter.save("figures/ProbabilityDistribution1s.png");

	plotter.clear();
	plotter.plot(probabilityDistribution.getSlice({1, IDX_ALL}));
	plotter.plot(
		probabilityDistribution.getSlice({3, IDX_ALL}),
		{{"linestyle", " "}, {"marker", "o"}, {"markersize", "5"}}
	);
	plotter.save("figures/ProbabilityDistribution2s.png");

	return 0;
}
