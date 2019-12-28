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

/** @file SelfConsistentField.cpp
 *
 *  @author Kristofer Björnson
 */

#include "SelfConsistentField.h"
#include "TBTK/Functions.h"
#include "TBTK/Plotter.h"

namespace TBTK{
namespace QTAT{
namespace Solver{

SelfConsistentField::SelfConsistentField(){
	gamma1 = 0.005;
	gamma2 = 0.005;
	temperature = 0.025/UnitHandler::getK_BN();
	mu1 = 0;
	mu2 = 2;
}

void SelfConsistentField::run(){
	auto f = Functions::fermiDiracDistribution;
	auto q = UnitHandler::getEN();

	Array<double> N({100});
	Array<double> I({100});
	Array<double> U0({100});
	for(unsigned int n = 0; n < 100; n++){
		const double V_G = UnitHandler::convertVoltageDtN(
			0,
			UnitHandler::VoltageUnit::V
		);
		const double V_D = UnitHandler::convertVoltageDtN(
			(int)n/100.,
			UnitHandler::VoltageUnit::V
		);

		Streams::out << V_D << " " << UnitHandler::getVoltageUnitString() << "\n";

		Property::DOS dos = calculateDOS();
		double U_L = -q*V_D/2.;
		double U_0 = 0.1;
		double U = U_L;
		mu2 = mu1 - q*V_D;
		Streams::out << "U_L:\t" << U_L << "\n";
		Streams::out << "mu2:\t" << mu2 << "\n";

		double dE = (
			dos.getUpperBound() - dos.getLowerBound()
		)/dos.getResolution();

		double N_0 = 0;
		for(int n = 0; n < dos.getResolution(); n++){
			double E = dos.getLowerBound() + n*dE;

			N_0 += dos(n)*f(E, 0, temperature)*dE;
		}

		double dU = 1;
		double numParticles;
		double current;
		while(dU > 1e-6){
			numParticles = 0;
			current = 0;
			for(int n = 0; n < dos.getResolution(); n++){
				double E = dos.getLowerBound() + n*dE;

				numParticles += dos(n)*(
					gamma1*f(E + U, mu1, temperature)
					+ gamma2*f(E + U, mu2, temperature)
				)/(gamma1 + gamma2)*dE;

				current += dos(n)*gamma1*gamma2/(gamma1 + gamma2)*(
					f(E + U, mu1, temperature)
					- f(E + U, mu2, temperature)
				)*dE;
			}
			current *= UnitHandler::getEN()/UnitHandler::getHbarN();

			double dN = numParticles - N_0;
			double newU = U_L + U_0*dN;
			dU = fabs(newU - U);
			U = newU;
			Streams::out << dU << "\n";
		}

		Plot::Plotter plotter;
		plotter.plot(dos);
		plotter.save("figures/DOS.png");

		Streams::out << "Number of particles:\t" << numParticles << "\n";
		Streams::out << "Current:\t" << current << "\n";

		N[{n}] = numParticles;
		I[{n}] = current;
		U0[{n}] = U;
	}

	Plot::Plotter plotter;
	plotter.plot(N);
	plotter.save("figures/NumParticles.png");
	plotter.plot(I);
	plotter.save("figures/Current.png");
	plotter.plot(U0);
	plotter.save("figures/U.png");
}

Property::DOS SelfConsistentField::calculateDOS(){
	const double gamma = gamma1 + gamma2;
	const double RESOLUTION = 1000;

	Property::DOS dos(-10, 10, RESOLUTION);
	double dE = (
			dos.getUpperBound() - dos.getLowerBound()
		)/dos.getResolution();
	for(unsigned int n = 0; n < RESOLUTION; n++){
		double energy = dos.getLowerBound() + n*dE;
		dos(n) = (gamma/(2.*M_PI))/(
			pow(energy - 0.2, 2) + pow(gamma/2, 2)
		);
	}

	return dos;
}

};	//End of namespace Solver
};	//End of namespace QTAT
};	//End of namespace TBTK
