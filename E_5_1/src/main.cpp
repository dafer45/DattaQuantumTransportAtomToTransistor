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
 *  @brief Exercise 5.1
 *
 *  Solution to exercise 5.1 in the book "Quantum Transport: Atom to
 *  Transistor, S. Datta (2005)".
 *
 *  @author Kristofer Björnson
 */

#include "TBTK/Model.h"
#include "TBTK/PropertyExtractor/BlockDiagonalizer.h"
#include "TBTK/Range.h"
#include "TBTK/Solver/BlockDiagonalizer.h"
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

	//Parameters.
	double a = 1;
	double E_0 = 0;
	double E_ss = 2;
	double Ep_ss = 1;

	//Setup the model.
	Model model;
	Range k(-M_PI/a, M_PI/a, 100);
	for(int n = 0; n < (int)k.getResolution(); n++){
		//Diagonal entries.
		model << HoppingAmplitude(E_0, {n, 0}, {n, 0});
		model << HoppingAmplitude(E_0, {n, 1}, {n, 1});

		//Off-diagonal entries.
		model << HoppingAmplitude(
			E_ss + Ep_ss*exp(-i*k[n]*a),
			{n, 0},
			{n, 1}
		) + HC;
	}
	model.construct();

	//Setup and run the solver.
	Solver::BlockDiagonalizer solver;
	solver.setModel(model);
	solver.run();

	//Extract the eigenvalues.
	PropertyExtractor::BlockDiagonalizer propertyExtractor(solver);
	Array<double> eigenValues({2, k.getResolution()});
	for(unsigned int n = 0; n < k.getResolution(); n++){
		eigenValues[{0, n}] = propertyExtractor.getEigenValue({(int)n}, 0);
		eigenValues[{1, n}] = propertyExtractor.getEigenValue({(int)n}, 1);
	}

	//Plot the eigenvalues.
	Plotter plotter;
	plotter.setLabelX("k (in units of pi/a)");
	plotter.setLabelY("Energy (eV)");
	plotter.plot(eigenValues.getSlice({0, IDX_ALL}));
	plotter.plot(eigenValues.getSlice({1, IDX_ALL}));
	plotter.save("figures/Spectrum.png");

	return 0;
}
