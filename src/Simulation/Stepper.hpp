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

        /** @class Enum for stepsizes
         *  @bug cannot use double value e.g. for week: 7.5 should be a struct instead
         */
        enum TimeFrame {
            DAY = 1,
            TWO_DAYS = 2,
            WEEK = 7,
            TEN_DAYS = 10,
            FORTNIGHT = 15,
	        MONTH = 30,
	        YEAR = 360
        };

        typedef std::pair<Solver::Equation *, double> step;

        /**
         * @class AbstractStepper An iterator in order to iterate simply over simulation steps
         * Holds a pointer to the equation and the choosen stepsize
         */
        class AbstractStepper {
        public:
            virtual Solver::Equation *
            get(int col) const = 0;
        };

        /**
         * @class Iterator The internal iterator holding the current simulation step
         */
        class Iterator {
        public:
            Iterator(const AbstractStepper *stepper, TimeFrame time, int steps, double pos, bool dynStep = false)
                    : _pos(pos), _stepper(stepper), _time(time), _dynStep(dynStep), _totalSteps(steps) {
                assert( ((dynStep) ? steps>1 : true) && "Dynamic steps not valid for 1 step");
                if(_dynStep){calcInit();}
            }

            /**
             * @note This should only be used in a for-each loop
             * This allows for a temporal error of 1e-2
             */
            bool operator!=(const Iterator &other) const {
                if (_dynStep) {
                    return not((_pos - 1e-2) >= _totalSteps and other._pos >= _totalSteps);
                } else {
                    return static_cast<int>(_pos) != static_cast<int>(other._pos);
                }
            }

            step operator*() const {
                return {_stepper->get(0), _pos};
            };

            const Iterator &operator++() {
                if (_dynStep) {
                    double __delta{0};
                    __delta = _delta_t_n * _p;
                    _delta_t_n = __delta;

                    LOG(debug) << "Stepsize delta " << __delta;
                    LOG(debug) << "Stepsize: " << _time *  __delta;
                    _stepper->get(0)->updateStepSize(_time * __delta);
                    _pos = _pos + __delta;
                    LOG(debug) << "Current position " << _pos;

                } else {
                    ++_pos;
                }
                return *this;
            }

            const void calcInit() {
                double __delta{0};
                __delta = _totalSteps * ((_p - 1) / (std::pow(_p, _totalSteps) - 1));
                _delta_t_n = __delta;
                _stepper->get(0)->updateStepSize(_time * __delta);
		            LOG(debug) << "Stepsize: " << _time * __delta;
                _pos = _pos + __delta;
            }

        private:
            double _pos; //delta t
            double _delta_t_n{0}; //last delta
            const double _p{1.2}; //The step multiplier
            const bool _dynStep{false};
            const int _totalSteps;
            const AbstractStepper *_stepper;
            const TimeFrame _time;
        };

        /**
         * @class Stepper An instance holding the simulation iterator
         */
        class Stepper : public AbstractStepper {
        public:

            Stepper(Solver::Equation *eq, const TimeFrame time, const size_t steps, bool dynStep = false)
                    : _equation(eq), _timeFrame(time), _steps(steps), _dyn(dynStep) {
                _equation->updateStepSize(_timeFrame);
            }

            virtual Solver::Equation *
            get(int col) const {
                return _equation;
            }

            Iterator
            begin() const {
                return Iterator(this, _timeFrame, this->_steps, 0, _dyn);
            }

            Iterator
            end() const {
                return Iterator(this, _timeFrame, this->_steps, this->_steps, _dyn);
            }

            const TimeFrame
            getStepSize() {
                return _timeFrame;
            };

        private:
            Solver::Equation *_equation;
            const TimeFrame _timeFrame;
            const size_t _steps;
            const bool _dyn;
        };

    }
}//ns
#endif