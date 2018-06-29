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
 *  @brief Exercise 3.1
 *
 *  Solution to exercise 3.1 in the book "Quantum Transport: Atom to
 *  Transistor, S. Datta (2005)".
 *
 *  @author Kristofer Björnson
 */

#include "TBTK/Model.h"
#include "TBTK/Plotter.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/Range.h"
#include "TBTK/Solver/Diagonalizer.h"
#include "TBTK/Streams.h"
#include "TBTK/UnitHandler.h"
#include "TBTK/Property/WaveFunctions.h"

#include <complex>

using namespace std;
using namespace TBTK;
using namespace Plot;

const int SIZE_R = 100;
const double dr = 0.05;
Range r(dr, dr*SIZE_R, SIZE_R);
Array<double> U_SCF;

//Callback function that returns U_SCF.
complex<double> callbackU_SCF(const Index &to, const Index &from){
	return U_SCF[{(unsigned int)from[0]}];
}

//Callback function for calculating U_SCF.
bool selfConsistencyCallback(Solver::Diagonalizer &solver){
	PropertyExtractor::Diagonalizer propertyExtractor(solver);

	//Calculate the density. (Only the first state contributes for He,
	//multiplied by a factor 2 for spin.)
	Array<double> density({SIZE_R}, 0);
	for(unsigned int state = 0; state < 1; state++){
		for(int r = 0; r < SIZE_R; r++){
			complex<double> amplitude
				= propertyExtractor.getAmplitude(state, {r});
			density[{(unsigned int)r}] += 2*pow(abs(amplitude), 2)/dr;
		}
	}

	//Parameters
	double e = UnitHandler::getEN();
	double epsilon_0 = UnitHandler::getEpsilon_0N();

	//Calculate the new U_SCF
	Array<double> newU_SCF({SIZE_R}, 0);
	for(unsigned int n = 0; n < SIZE_R; n++){
		for(unsigned int np = 0; np < n; np++){
			newU_SCF[{n}] += (1/2.)*e*e/(4*M_PI*epsilon_0*r[n])*dr*density[{np}];
		}
		for(unsigned int np = n; np < SIZE_R; np++){
			newU_SCF[{n}] += (1/2.)*e*e/(4*M_PI*epsilon_0)*dr*density[{np}]/r[np];
		}
	}

	//Calculate the difference between the new and old solution.
	double difference = 0;
	for(unsigned int n = 0; n < SIZE_R; n++)
		difference += abs(newU_SCF[{n}] - U_SCF[{n}]);

	//Update U_SCF.
	U_SCF = newU_SCF;

	//Return true to indicate convergence if the difference between the new
	//and old solution is smaller than 10 meV.
	if(difference < 10e-3)
		return true;
	else
		return false;
}

int main(int argc, char **argv){
	//Set the natural units. Argument order: (charge, count, energy,
	//length, temperature, time).
	UnitHandler::setScales({"1 C", "1 pcs", "1 eV", "1 Ao", "1 K", "1 s"});

	//Parameters.
	double hbar = UnitHandler::getHbarN();
	double m_e = UnitHandler::getM_eN();
	double e = UnitHandler::getEN();
	double epsilon_0 = UnitHandler::getEpsilon_0N();
	double A = hbar*hbar/(2*m_e*dr*dr);
	double B = 2*e*e/(4*M_PI*epsilon_0);
	double C = hbar*hbar/(2*m_e);
	const int SIZE_R = 100;
	double l = 0;

	//Initial guess for U_SCF.
	U_SCF = Array<double>({SIZE_R}, 0);

	//Setup the model.
	Model model;
	for(unsigned int n = 0; n < r.getResolution(); n++){
		model << HoppingAmplitude(
			2*A - B/r[n] + l*(l+1)*C/(r[n]*r[n]),
			{(int)n},
			{(int)n}
		);

		model << HoppingAmplitude(callbackU_SCF, {(int)n}, {(int)n});

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
	solver.setSelfConsistencyCallback(selfConsistencyCallback);
	solver.setMaxIterations(100);
	solver.run();

	//Calculate the 1s energy.
	PropertyExtractor::Diagonalizer propertyExtractor(solver);
	Property::EigenValues eigenValues = propertyExtractor.getEigenValues();
	Streams::out << "1s energy:\t" << eigenValues(0) << "\n";

	//Plot the potentials U_N and U_SCF.
	Array<double> U_N({SIZE_R});
	for(unsigned int n = 0; n < SIZE_R; n++)
		U_N[{n}] = -B/r[n];
	Plotter plotter;
	plotter.setBoundsY(-100, 20);
	plotter.setLabelX("r");
	plotter.setLabelY("U (eV)");
	plotter.setHold(true);
	plotter.plot(U_N);
	plotter.plot(
		U_SCF,
		Decoration({0, 0, 0}, Decoration::LineStyle::Line, 2)
	);
	plotter.save("figures/Potentials.png");

	//Extract the wave functions for the 1s state. Note that zero based
	//indexing is used.
	Property::WaveFunctions waveFunctions
		= propertyExtractor.calculateWaveFunctions(
			{{IDX_ALL}},
			{0}
		);

	//Calculate the probability distribution.
	Array<double> probabilityDistribution({SIZE_R});
	for(unsigned int n = 0; n < SIZE_R; n++){
		probabilityDistribution[{n}]
			= pow(abs(waveFunctions({(int)n}, 0)), 2);
	}

	//Plot the probability distributions.
	plotter.clear();
	plotter.setAutoScaleY(true);
	plotter.setLabelX("r");
	plotter.setLabelY("Probability density");

	plotter.plot(probabilityDistribution);
	plotter.setHold(true);
	plotter.save("figures/ProbabilityDistribution1s.png");

	return 0;
}
