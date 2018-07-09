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
#include "TBTK/Plotter.h"
#include "TBTK/PropertyExtractor/BlockDiagonalizer.h"
#include "TBTK/Range.h"
#include "TBTK/Solver/BlockDiagonalizer.h"
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

	//Parameter.
	double t = 3;

	//Calculate and plot Eq 6.1.12 for m=65 and m=66.
	for(unsigned int m = 65; m < 67; m++){
		//calculate.
		Array<double> E({2, 2, 101}, 0);
		for(unsigned int mode = 0; mode < 2; mode++){
			int nu = round(2*m/3) + mode;

			//kxa means k_x*a.
			Range kxa(-0.05*M_PI, 0.05*M_PI, 101);
			for(unsigned int n = 0; n < kxa.getResolution(); n++){
				//Eq 6.1.12 rewritten by using Eq. 5.2.10 to
				//eliminate a_0 and b in terms of a.
				E[{0, mode, n}] = t*sqrt(
					pow(kxa[n], 2)
					+ pow(M_PI*nu/m - 2*M_PI/3, 2)/3
				);
				E[{1, mode, n}] = -t*sqrt(
					pow(kxa[n], 2)
					+ pow(M_PI*nu/m - 2*M_PI/3, 2)/3
				);
			}
		}

		//Plot.
		Plotter plotter;
		plotter.setLabelX("k_x*a/pi");
		plotter.setLabelY("Energy (eV)");
		plotter.setHold(true);
		plotter.plot(E.getSlice({0, 0, IDX_ALL}));
		plotter.plot(E.getSlice({1, 0, IDX_ALL}));
		plotter.plot(E.getSlice({0, 1, IDX_ALL}));
		plotter.plot(E.getSlice({1, 1, IDX_ALL}));
		plotter.save("figures/Spectrum_m_" + to_string(m) + ".png");
	}

	return 0;
}
