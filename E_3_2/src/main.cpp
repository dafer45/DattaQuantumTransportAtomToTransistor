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
 *  @brief Exercise 3.2
 *
 *  Solution to exercise 3.2 in the book "Quantum Transport: Atom to
 *  Transistor, S. Datta (2005)".
 *
 *  @author Kristofer Björnson
 */

#include "TBTK/Model.h"
#include "TBTK/Property/WaveFunctions.h"
#include "TBTK/PropertyExtractor/BlockDiagonalizer.h"
#include "TBTK/Range.h"
#include "TBTK/Solver/BlockDiagonalizer.h"
#include "TBTK/Streams.h"
#include "TBTK/UnitHandler.h"
#include "TBTK/Visualization/MatPlotLib/Plotter.h"

#include <complex>

using namespace std;
using namespace TBTK;
using namespace Visualization::MatPlotLib;

const int SIZE_R = 100;
const double dr = 0.05;
Range r(dr, dr*SIZE_R, SIZE_R);
Array<double> U_SCF;

//Callback that returns U_SCF.
class CallbackU_SCF : public HoppingAmplitude::AmplitudeCallback{
	complex<double> getHoppingAmplitude(
		const Index &to,
		const Index &from
	) const{
		return U_SCF[{(unsigned int)from[1]}];
	}
} callbackU_SCF;

//Callback for calculating U_SCF.
class SelfConsistencyCallback :
	public Solver::BlockDiagonalizer::SelfConsistencyCallback
{
	bool selfConsistencyCallback(Solver::BlockDiagonalizer &solver){
		PropertyExtractor::BlockDiagonalizer propertyExtractor(solver);

		//Calculate the density. The first and second loop adds
		//contributions from the s and p states, respectively.
		Array<double> density({SIZE_R}, 0);
		for(unsigned int state = 0; state < 3; state++){
			for(int r = 0; r < SIZE_R; r++){
				complex<double> amplitude
					= propertyExtractor.getAmplitude(
						{0},
						state,
						{r}
					);
				density[{(unsigned int)r}]
					+= 2*pow(abs(amplitude), 2)/dr;
			}
		}
		for(unsigned int state = 0; state < 2; state++){
			double numOccupiedStates;
			if(state == 0)
				numOccupiedStates = 6;
			else
				numOccupiedStates = 2;

			for(int r = 0; r < SIZE_R; r++){
				complex<double> amplitude
					= propertyExtractor.getAmplitude(
						{1},
						state,
						{r}
					);
				density[{(unsigned int)r}]
					+= numOccupiedStates*pow(
						abs(amplitude),
						2
					)/dr;
			}
		}

		//Parameters
		double e = UnitHandler::getConstantInNaturalUnits("e");
		double epsilon_0
			= UnitHandler::getConstantInNaturalUnits("epsilon_0");

		//Calculate the new U_SCF
		Array<double> newU_SCF({SIZE_R}, 0);
		for(unsigned int n = 0; n < SIZE_R; n++){
			for(unsigned int np = 0; np < n; np++){
				newU_SCF[{n}] += (13/14.)*e*e/(
					4*M_PI*epsilon_0*r[n]
				)*dr*density[{np}];
			}
			for(unsigned int np = n; np < SIZE_R; np++){
				newU_SCF[{n}] += (13/14.)*e*e/(
					4*M_PI*epsilon_0
				)*dr*density[{np}]/r[np];
			}
		}

		//Calculate the difference between the new and old solution.
		double difference = 0;
		for(unsigned int n = 0; n < SIZE_R; n++)
			difference += abs(newU_SCF[{n}] - U_SCF[{n}]);

		//Update U_SCF.
		U_SCF = newU_SCF;

		//Return true to indicate convergence if the difference between
		//the new and old solution is smaller than 10 meV.
		if(difference < 10e-3)
			return true;
		else
			return false;
	}
} selfConsistencyCallback;

int main(int argc, char **argv){
	//Initialize TBTK.
	Initialize();

	//Set the natural units. Argument order: (charge, count, energy,
	//length, temperature, time).
	UnitHandler::setScales(
		{"1 rad", "1 C", "1 pcs", "1 eV", "1 Ao", "1 K", "1 s"}
	);

	//Parameters.
	double hbar = UnitHandler::getConstantInNaturalUnits("hbar");
	double m_e = UnitHandler::getConstantInNaturalUnits("m_e");
	double e = UnitHandler::getConstantInNaturalUnits("e");
	double epsilon_0 = UnitHandler::getConstantInNaturalUnits("epsilon_0");
	double A = hbar*hbar/(2*m_e*dr*dr);
	double B = 14*e*e/(4*M_PI*epsilon_0);
	double C = hbar*hbar/(2*m_e);

	//Initial guess for U_SCF.
	U_SCF = Array<double>({SIZE_R}, 0);

	//Setup the model.
	Model model;
	for(int l = 0; l < 2; l++){
		for(unsigned int n = 0; n < r.getResolution(); n++){
			model << HoppingAmplitude(
				2*A - B/r[n] + l*(l+1)*C/(r[n]*r[n]),
				{l, (int)n},
				{l, (int)n}
			);

			model << HoppingAmplitude(callbackU_SCF, {l, (int)n}, {l, (int)n});

			if(n+1 < r.getResolution()){
				model << HoppingAmplitude(
					-A,
					{l, (int)(n+1)},
					{l, (int)n}
				) + HC;
			}
		}
	}
	model.construct();

	//Setup and run solver.
	Solver::BlockDiagonalizer solver;
	solver.setModel(model);
	solver.setSelfConsistencyCallback(selfConsistencyCallback);
	solver.setMaxIterations(1000);
	solver.run();

	//Calculate the 1s energy.
	PropertyExtractor::BlockDiagonalizer propertyExtractor(solver);
	Property::EigenValues eigenValues = propertyExtractor.getEigenValues();
	Streams::out << "1s energy:\t" << eigenValues(0) << "\n";

	//Plot the potentials U_N and U_SCF.
	Array<double> U_N({SIZE_R});
	for(unsigned int n = 0; n < SIZE_R; n++)
		U_N[{n}] = -B/r[n];

	//Calculate the probability distribution.
	Array<double> probabilityDistribution({2, SIZE_R});
	for(unsigned int n = 0; n < SIZE_R; n++){
		//1s state.
		complex<double> amplitude
			= propertyExtractor.getAmplitude({0}, 0, {(int)n});
		probabilityDistribution[{0, n}] = pow(abs(amplitude), 2);

		//3p state.
		amplitude = propertyExtractor.getAmplitude({1}, 1, {(int)n});
		probabilityDistribution[{1, n}] = pow(abs(amplitude), 2);
	}

	//Plot the probability distributions.
	Plotter plotter;
	plotter.setBoundsY(0, 0.08);
	plotter.setLabelX("r");
	plotter.setLabelY("Probability density");
	plotter.plot(probabilityDistribution.getSlice({0, IDX_ALL}));
	plotter.plot(
		probabilityDistribution.getSlice({1, IDX_ALL}),
		{{"linestyle", " "}, {"marker", "o"}, {"markersize", "5"}}
	);
	plotter.save("figures/ProbabilityDistribution.png");

	return 0;
}
