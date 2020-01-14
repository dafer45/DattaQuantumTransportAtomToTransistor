/* Copyright 2016 Kristofer Björnson
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

/** @package TBTKtemp
 *  @file main.cpp
 *  @brief New project
 *
 *  Empty template project.
 *
 *  @author Kristofer Björnson
 */

#include "SelfConsistentField.h"
#include "TBTK/Streams.h"
#include "TBTK/TBTK.h"
#include "TBTK/UnitHandler.h"

#include <complex>

using namespace std;
using namespace TBTK;

int main(int argc, char **argv){
	//Initialize TBTK.
	Initialize();

	//Set the natural units. Argument order: (charge, count, energy,
	//length, temperature, time).
	UnitHandler::setScales(
		{"1 rad", "1 C", "1 pcs", "1 eV", "1 m", "1 K", "1 s"}
	);

	//Get the fundamental charge in natural units.
	double q = UnitHandler::getConstantInNaturalUnits("e");

	//Convert voltage from V to natural units.
	double V_D = UnitHandler::convertArbitraryToNatural<Quantity::Voltage>(
		0.2,
		Quantity::Voltage::Unit::V
	);

	//Calculate the temperature T = 0.025/k_B.
	double temperature
		= 0.025/UnitHandler::getConstantInNaturalUnits("k_B");

	QTAT::Solver::SelfConsistentField solver;

	solver.setMu1(0);
	solver.setMu2(0 - q*V_D);
	solver.setTemperature(temperature);
	solver.setGamma1(0.05);
	solver.run();

	return 0;
}
