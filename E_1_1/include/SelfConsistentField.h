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

/** @package TBTKcalc
 *  @file SelfConsistentField.h
 *  @brief Base class for solving the self-consistent field equation outlined
 *  in Quantum Transport - Atom to Transistor, S. Datta (2005).
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_TBTK_QTAT_SELF_CONSISTENT_FIELD
#define COM_DAFER45_TBTK_QTAT_SELF_CONSISTENT_FIELD

#include "TBTK/Property/DOS.h"
#include "TBTK/Solver/Solver.h"

namespace TBTK{
namespace QTAT{
namespace Solver{

class SelfConsistentField : public TBTK::Solver::Solver{
public:
	/** Constructor. */
	SelfConsistentField();

	/***/
	void run();

	/** */
	virtual Property::DOS calculateDOS();

	/** Set gamma1. */
	void setGamma1(double gamma1);

	/** Set gamma2. */
	void setGamma2(double gamma1);

	/** Set temperature. */
	void setTemperature(double temperature);

	/** Set mu1. */
	void setMu1(double mu1);

	/** Set mu2. */
	void setMu2(double mu2);
private:
	double gamma1;
	double gamma2;
	double temperature;
	double mu1;
	double mu2;
};

inline void SelfConsistentField::setGamma1(double gamma1){
	this->gamma1 = gamma1;
}

inline void SelfConsistentField::setGamma2(double gamma2){
	this->gamma2 = gamma2;
}

inline void SelfConsistentField::setTemperature(double temperature){
	this->temperature = temperature;
}

inline void SelfConsistentField::setMu1(double mu1){
	this->mu1 = mu1;
}

inline void SelfConsistentField::setMu2(double mu2){
	this->mu2 = mu2;
}

};	//End namespace Solver
};	//End of namespace QTAT
};	//End namespace TBTK

#endif
