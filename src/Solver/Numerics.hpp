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

#ifndef GLOBAL_FLOW_NUMERICS_HPP
#define GLOBAL_FLOW_NUMERICS_HPP
namespace GlobalFlow {
namespace Solver {


using pr_t = double;
using vector = Eigen::Matrix<pr_t, Eigen::Dynamic, 1>;

/**
 * @class AdaptiveDamping
 * Can be used to dampen the change per iteration
 */
class AdaptiveDamping {

    public:

        AdaptiveDamping() : min_damping(0), max_damping(0), max_allowed_change(0.01) {};

        AdaptiveDamping(pr_t min_damping, pr_t max_damping, pr_t max_allowed_change, vector x_t0)
                : min_damping(min_damping), max_damping(max_damping), max_allowed_change(max_allowed_change), x_t0(x_t0) {
            damping_t0 = sqrt(min_damping * max_damping);
        }


        vector getChanges(vector &residuals, vector &x, bool isAdaptiveDamping) {
            assert(x.rows() == x_t0.rows() && "Damping hasn't been properly initialized");
            vector changes = x - x_t0;
            pr_t damping{1.0};
            if (isAdaptiveDamping) {
                damping = getDamping(changes.maxCoeff(), getDampNorm(residuals, x));
                LOG(numerics) << "Damping applied: " << damping;
            }
            x_t0 = x;
            return changes * damping;
        }

        pr_t getNorm() { return norm_t0; }

    private:
        bool firstOuterIteration{true};

        vector x_t0;
        pr_t norm_t0;
        pr_t max_change_t0;
        pr_t damping_t0;

        pr_t min_damping;
        pr_t max_damping;
        pr_t max_allowed_change;


        int cnt{0};

        /**
         * @brief Get the norm of the current residuals for adaptive damping
         * @note Use only once per Iteration! Modifies x_t0
         */
        pr_t getDampNorm(vector &residuals, vector &x) {
            vector changes = x_t0 - x;
            pr_t res = residuals.transpose() * residuals;
            pr_t ch = changes.transpose() * changes;
            pr_t norm = sqrt(res * ch);

            LOG(numerics) << "square root norm of residuals: " << norm;
            return norm;
        }


        pr_t getDamping(pr_t max_change, pr_t norm) {
            pr_t PHI{0.01};
            pr_t sigma{0};

            pr_t p_n;
            pr_t p_h;

            if (firstOuterIteration) {
                p_n = norm;
                p_h = max_change;
                firstOuterIteration = false;
            } else {
                p_n = norm / norm_t0;
                p_h = max_change / max_change_t0;
            }

            if (p_n < 1 and p_h < 1) {
                auto lambda = log10(p_n) / log10(PHI);
                if (lambda < 1) {
                    sigma = damping_t0 + lambda * (max_damping - damping_t0);
                } else {
                    sigma = max_damping;
                }
                cnt = 0;
            }

            if (p_n > 1) {
                sigma = damping_t0 / p_n;
            }

            if (p_h > 1) {
                sigma = damping_t0 / p_h;
            }

            pr_t damping_t{sqrt(sigma * damping_t0)};

            if (std::abs(max_change) > max_allowed_change and damping_t > (max_allowed_change / std::abs(max_change))) {
                damping_t = max_allowed_change / std::abs(max_change);
            }
            if (damping_t < min_damping) {
                damping_t = min_damping;
                cnt++;
                if (cnt > 10) {
                    damping_t = pow(pow(min_damping, 2) * max_damping, 1.0 / 3.0);
                }
            }

            damping_t0 = damping_t;
            norm_t0 = norm;
            max_change_t0 = max_change;
            return damping_t;
        }

};

}
}//ns
#endif //GLOBAL_FLOW_NUMERICS_HPP
