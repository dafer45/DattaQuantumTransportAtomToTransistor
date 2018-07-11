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
 *  @brief Exercise 7.2.a
 *
 *  Solution to exercise 7.2.a in the book "Quantum Transport: Atom to
 *  Transistor, S. Datta (2005)".
 *
 *  @author Kristofer Björnson
 */

#include "TBTK/Model.h"
#include "TBTK/Plotter.h"
#include "TBTK/Property/Density.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/Range.h"
#include "TBTK/Solver/Diagonalizer.h"
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
	double epsilon_0 = UnitHandler::getEpsilon_0N();
	double k_B = UnitHandler::getK_BN();
	double e = UnitHandler::getEN();
	double a = 3;	//Ångström
	double mu = 0;	//eV
	double temperature = 300;	//Kelvin

	//Widths.
	unsigned int width_oxide = 21/a;
	unsigned int width_silicon = 30/a;
	unsigned int width_total = 2*width_oxide + width_silicon;

	//Energies.
	double E_silicon = 0;	//eV
	double E_oxide = 3;	//eV

	//Dielectric constant.
	double epsilon = 4*epsilon_0;

	//Mass.
	double m_c = 0.25*m_e;

	//Gate voltage.
	double V_g = UnitHandler::convertVoltageDtN(
		0.25,
		UnitHandler::VoltageUnit::V
	);

	//Constant in Eq. 7.2.8.
	double N_0 = m_c*k_B*temperature/(2*M_PI*hbar*hbar);

	//Setup the electric potential.
	Array<double> U({width_total}, -e*V_g);
	Array<double> density;

	//Self-consistency loop.
	double maxDifference = 1;
	while(maxDifference > 1e-3){
		/////////////////
		// Schrödinger //
		/////////////////

		//Create the model.
		Model schrodingerModel;
		schrodingerModel.setVerbose(false);

		//Prefactor for the discretized second derivative.
		double t_silicon = hbar*hbar/(2*m_c*a*a);
		double t_oxide = hbar*hbar/(2*m_c*a*a);

		//Add diagonal elements.
		for(unsigned int n = 0; n < width_total; n++){
			if(n < width_oxide || n >= width_oxide + width_silicon){
				//Oxide.
				schrodingerModel << HoppingAmplitude(
					E_oxide + U[{n}] + 2*t_oxide,
					{n},
					{n}
				);
			}
			else if(
				n == width_oxide
				|| n == width_oxide + width_silicon - 1
			){
				//Boundary sites.
				schrodingerModel << HoppingAmplitude(
					(E_oxide + E_silicon)/2.
					+ U[{n}] + t_oxide + t_silicon,
					{n},
					{n}
				);
			}
			else{
				//Silicon.
				schrodingerModel << HoppingAmplitude(
					E_silicon + U[{n}]
					+ 2*t_silicon,
					{n},
					{n}
				);
			}
		}

		//Add off-diagonal elements.
		for(unsigned int n = 0; n < width_total - 1; n++){
			if(
				n < width_oxide
				|| n >= width_oxide + width_silicon - 1
			){
				//Oxide.
				schrodingerModel << HoppingAmplitude(
					-t_oxide,
					{n+1},
					{n}
				) + HC;
			}
			else{
				//Silicon.
				schrodingerModel << HoppingAmplitude(
					-t_silicon,
					{n+1},
					{n}
				) + HC;
			}
		}

		//Construct the Hilbert space basis.
		schrodingerModel.construct();

		//Setup and run the solver.
		Solver::Diagonalizer schrodingerSolver;
		schrodingerSolver.setVerbose(false);
		schrodingerSolver.setModel(schrodingerModel);
		schrodingerSolver.run();

		//Calculate the density.
		PropertyExtractor::Diagonalizer propertyExtractor(
			schrodingerSolver
		);
		Property::EigenValues eigenValues
			= propertyExtractor.getEigenValues();
		Property::WaveFunctions waveFunctions
			= propertyExtractor.calculateWaveFunctions(
				{{IDX_ALL}},
				{IDX_ALL}
			);
		density = Array<double>({width_total}, 0);
		for(unsigned int n = 0; n < width_total; n++){
			for(
				unsigned int state = 0;
				state < width_total;
				state++
			){
				complex<double> amplitude
					= waveFunctions({n}, state);
				double energy = eigenValues(state);

				//Eq 7.2.13.
				density[{n}] += pow(
					abs(amplitude), 2
				)*2*N_0*log(
					1 + exp(
						(mu - energy)/(k_B*temperature)
					)
				)/a;
			}
		}

		/////////////
		// Poisson //
		/////////////

		//Create the model.
		Model poissonModel;
		poissonModel.setVerbose(false);
		for(unsigned int n = 0; n < width_total; n++){
			//Add diagonal elements.
			poissonModel << HoppingAmplitude(
				2*epsilon/(a*a),
				{n},
				{n}
			);

			//Add off-diagonal elements.
			if(n+1 < width_total){
				poissonModel << HoppingAmplitude(
					-epsilon/(a*a),
					{n+1},
					{n}
				) + HC;
			}

			//Add the right hand side of the equation (except for
			//boundary terms).
			poissonModel << SourceAmplitude(
				e*e*density[{n}],
				{n}
			);
		}

		//Add boundary terms.
		poissonModel << SourceAmplitude(-e*V_g*epsilon/(a*a), {0});
		poissonModel << SourceAmplitude(-e*V_g*epsilon/(a*a), {width_total-1});

		//Create the Hilbert space basis.
		poissonModel.construct();

		//Setup and run the solver.
		Solver::LinearEquationSolver poissonSolver;
		poissonSolver.setVerbose(false);
		poissonSolver.setModel(poissonModel);
		poissonSolver.run();

		//Extract the potential.
		Array<double> newU({width_total});
		for(unsigned int n = 0; n < width_total; n++)
			newU[{n}] = real(poissonSolver.getAmplitude({n}));

		//Calculate the maximum difference between the old an new
		//solution.
		maxDifference = 0;
		Array<double> difference = U - newU;
		for(unsigned int n = 0; n < width_total; n++)
			if(abs(difference[{n}]) > maxDifference)
				maxDifference = abs(difference[{n}]);

		//Replace U with a mix of the old and new U.
		U = 0.1*newU + 0.9*U;
	}

	//Plot the density.
	Plotter plotter;
	plotter.setBoundsY(0, 1.5e-5);
	plotter.setLabelX("z");
	plotter.setLabelY("Density (Ao^-3)");
	plotter.plot(density);
	plotter.save("figures/Density.png");

	//Calculate the conduction band profile.
	Array<double> conductionBandProfile = U;
	for(unsigned int n = 0; n < width_total; n++){
		if(n < width_oxide || n >= width_oxide + width_silicon)
			conductionBandProfile[{n}] += E_oxide;
		else
			conductionBandProfile[{n}] += E_silicon;
	}

	//Plot the conduction band profile.
	plotter.setBoundsY(-0.5, 3);
	plotter.setLabelX("z");
	plotter.setLabelY("Energy (eV)");
	plotter.plot(conductionBandProfile);
	plotter.save("figures/ConductionBandProfile.png");

	return 0;
}
