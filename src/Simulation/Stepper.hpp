/*
 * Copyright (c) <2016>, <Robert Reinecke>
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef GLOBAL_FLOW_STEPPER_HPP
#define GLOBAL_FLOW_STEPPER_HPP

#include "../Solver/Equation.hpp"

namespace GlobalFlow {
    namespace Simulation {

        /** @class Enum for step-sizes
         *  @bug cannot use double value e.g. for week: 7.5 should be a struct instead
         */
        enum TimeFrame {
            DAY = 1,
            TWO_DAYS = 2,
            FOUR_DAYS = 4,
            WEEK = 7,
            TEN_DAYS = 10,
            FORTNIGHT = 15,
	        MONTH = 30,
	        YEAR = 365,
            TWO_YEARS = YEAR * 2,
            TEN_YEARS = YEAR * 10,
            HUNDRED_YEARS = YEAR * 100,
            THOUSAND_YEARS = YEAR * 1000
        };

        typedef std::pair<Solver::Equation *, double> step;

        /**
         * @class AbstractStepper An iterator in order to iterate simply over simulation steps
         * Holds a pointer to the equation and the chosen step-size
         */
        class AbstractStepper {
        public:
            virtual Solver::Equation *get(int col) const = 0;
        };

        /**
         * @class Iterator The internal iterator holding the current simulation step
         */
        class Iterator {
        public:
            //Iterator(const AbstractStepper *stepper, TimeFrame time, int steps, double pos, bool dynStep = false)
            Iterator(const AbstractStepper *stepper, int step_size, int steps, double delta_t, bool dynStep = false)
                    : _delta_t(delta_t), _stepper(stepper), _step_size(step_size), _dynStep(dynStep), _totalSteps(steps) {
                assert( ((dynStep) ? steps>1 : true) && "Dynamic steps not valid for 1 step");
                if(_dynStep){calcInit();}
            }

            /**
             * @note This should only be used in a for-each loop
             * This allows for a temporal error of 1e-2
             */
            bool operator!=(const Iterator &other) const {
                if (_dynStep) {
                    return not((_delta_t - 1e-2) >= _totalSteps and other._delta_t >= _totalSteps);
                } else {
                    return static_cast<int>(_delta_t) != static_cast<int>(other._delta_t);
                }
            }

            step operator*() const {
                return {_stepper->get(0), _delta_t};
            };

            const Iterator &operator++() {
                if (_dynStep) {
                    double __delta{0};
                    __delta = _delta_t0 * _p;
                    _delta_t0 = __delta;

                    LOG(debug) << "Stepsize delta " << __delta;
                    LOG(debug) << "Stepsize: " << double(_step_size) *  __delta;
                    _stepper->get(0)->updateStepSize(double(_step_size) * __delta);
                    _delta_t = _delta_t + __delta;
                    LOG(debug) << "Current position " << _delta_t;

                } else {
                    ++_delta_t;
                }
                return *this;
            }

            void calcInit() {
                double __delta{0};
                __delta = _totalSteps * ((_p - 1) / (std::pow(_p, _totalSteps) - 1));
                _delta_t0 = __delta;
                _stepper->get(0)->updateStepSize(double(_step_size) * __delta);
		            LOG(debug) << "Stepsize: " << double(_step_size) * __delta;
                _delta_t = _delta_t + __delta;
            }

        private:
            double _delta_t{0}; //delta t
            double _delta_t0{0}; //last delta t
            const double _p{1.2}; //The step multiplier
            const bool _dynStep{false};
            const int _totalSteps;
            const AbstractStepper *_stepper;
            const int _step_size;
        };

        /**
         * @class Stepper An instance holding the simulation iterator
         */
        class Stepper : public AbstractStepper {
        public:
            //Stepper(Solver::Equation *eq, const TimeFrame time, const size_t steps, bool dynStep = false)
            Stepper(Solver::Equation *eq, const int stepSize, bool isSteadyState, bool isDensityVariable,
                    const size_t steps, bool dynStep = false)
                    : _equation(eq), _stepSize(stepSize), _isSteadyState(isSteadyState),
                      _isDensityVariable(isDensityVariable), _steps(steps), _dyn(dynStep) {
                _equation->updateStepSize(_stepSize);
                _equation->updateIsSteadyState(_isSteadyState);
                _equation->updateIsDensityVariable(_isDensityVariable);
            }

            Stepper(Solver::Equation *eq, const std::string& stepSize, bool isSteadyState,
                    bool isDensityVariable, const size_t steps, bool dynStep = false)
                    : _equation(eq), _stepSize(getStepSizeWithString(stepSize)), _isSteadyState(isSteadyState),
                      _isDensityVariable(isDensityVariable), _steps(steps), _dyn(dynStep) {
                _equation->updateStepSize(_stepSize);
                _equation->updateIsSteadyState(_isSteadyState);
                _equation->updateIsDensityVariable(_isDensityVariable);
            }

            virtual Solver::Equation *
            get(int col) const {
                return _equation;
            }

            Iterator
            begin() const {
                return Iterator(this, _stepSize, this->_steps, 0, _dyn);
            }

            Iterator
            end() const {
                return Iterator(this, _stepSize, this->_steps, this->_steps, _dyn);
            }

            const int getStepSize() {
                return _stepSize;
            };

            const int getStepSizeWithString(const std::string& stepSizeString){
                if (stepSizeString == "DAY") { return DAY; }
                else if (stepSizeString == "TWO_DAYS") { return TWO_DAYS; }
                else if (stepSizeString == "FOUR_DAYS") { return FOUR_DAYS; }
                else if (stepSizeString == "WEEK") { return WEEK; }
                else if (stepSizeString == "TEN_DAYS") { return TEN_DAYS; }
                else if (stepSizeString == "FORTNIGHT") { return FORTNIGHT; }
                else if (stepSizeString == "MONTH") { return MONTH; }
                else if (stepSizeString == "YEAR") { return YEAR; }
                else if (stepSizeString == "TWO_YEARS") { return TWO_YEARS; }
                else if (stepSizeString == "TEN_YEARS") { return TEN_YEARS; }
                else if (stepSizeString == "HUNDRED_YEARS") { return HUNDRED_YEARS; }
                else if (stepSizeString == "THOUSAND_YEARS") { return THOUSAND_YEARS; }
                else {throw "Provided step size " + stepSizeString + " not availabe";}
            }

        private:
            Solver::Equation *_equation;
            int _stepSize;
            const bool _isSteadyState;
            const bool _isDensityVariable;
            const size_t _steps;
            const bool _dyn;
        };

    }
}//ns
#endif
