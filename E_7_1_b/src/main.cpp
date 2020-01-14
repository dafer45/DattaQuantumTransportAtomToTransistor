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
 *  @brief Exercise 7.1.b
 *
 *  Solution to exercise 7.1.b in the book "Quantum Transport: Atom to
 *  Transistor, S. Datta (2005)".
 *
 *  @author Kristofer Björnson
 */

#include "TBTK/Model.h"
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
	double hbar = UnitHandler::getConstantInNaturalUnits("hbar");
	double m_e = UnitHandler::getConstantInNaturalUnits("m_e");
	double a = 3;	//Ångström
	int width_AlGaAs = 100;

	//Masses.
	double m_GaAs = 0.07*m_e;
	double m_AlAs = 0.15*m_e;
	double m_AlGaAs = 0.3*m_AlAs + 0.7*m_GaAs;

	//Energies.
	double E_GaAs = 0;	//eV
	double E_AlAs = 1.25;	//eV
	double E_AlGaAs = 0.3*E_AlAs + 0.7*E_GaAs;

	//Calculate the two lowest eigenvalues as a function of the GaAs width.
	Array<double> lowestEigenValues({2, 24});
	for(int width = 10; width < 34; width++){
		Model model;

		//Prefactor for the discretized second derivative.
		double t_AlGaAs = hbar*hbar/(2*m_AlGaAs*a*a);
		double t_GaAs = hbar*hbar/(2*m_GaAs*a*a);

		//Add diagonal elements.
		for(int n = 0; n < 2*width_AlGaAs + width; n++){
			if(n < width_AlGaAs || n >= width_AlGaAs + width){
				//AlGaAs.
				model << HoppingAmplitude(
					E_AlGaAs + 2*t_AlGaAs,
					{n},
					{n}
				);
			}
			else if(
				n == width_AlGaAs
				|| n == width_AlGaAs + width - 1
			){
				//Boundary sites.
				model << HoppingAmplitude(
					(E_AlGaAs + E_GaAs)/2.
					+ t_AlGaAs + t_GaAs,
					{n},
					{n}
				);
			}
			else{
				//GaAs.
				model << HoppingAmplitude(
					E_GaAs + 2*t_GaAs,
					{n},
					{n}
				);
			}
		}

		//Add off-diagonal elements.
		for(int n = 0; n < 2*width_AlGaAs + width - 1; n++){
			if(n < width_AlGaAs || n >= width_AlGaAs + width - 1){
				//AlGaAs.
				model << HoppingAmplitude(
					-t_AlGaAs,
					{n+1},
					{n}
				) + HC;
			}
			else{
				//GaAs.
				model << HoppingAmplitude(
					-t_GaAs,
					{n+1},
					{n}
				) + HC;
			}
		}

		//Construct the Hilbert space basis.
		model.construct();

		//Setup and run solver.
		Solver::Diagonalizer solver;
		solver.setModel(model);
		solver.run();

		//Extract eigenvalues.
		PropertyExtractor::Diagonalizer propertyExtractor(solver);
		Property::EigenValues eigenValues
			= propertyExtractor.getEigenValues();

		//Store the two lowest eigenvalues.
		unsigned int n = width - 10;
		lowestEigenValues[{0, n}] = eigenValues(0);
		lowestEigenValues[{1, n}] = eigenValues(1);
	}

	//Plot the two lowest eigenvalues as a function of GaAs width.
	Plotter plotter;
	plotter.setBoundsY(0, 0.4);
	plotter.setLabelX("Width");
	plotter.setLabelY("Energy (eV)");
	plotter.plot(lowestEigenValues.getSlice({0, IDX_ALL}));
	plotter.plot(lowestEigenValues.getSlice({1, IDX_ALL}));
	plotter.save("figures/EigenValues.png");

	return 0;
}
