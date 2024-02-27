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

        AdaptiveDamping() : Sigma_MIN(0), Sigma_MAX(0), Change_MAX(0.01) {};

        AdaptiveDamping(pr_t Sigma_Min, pr_t Sigma_Max, pr_t Change_Max, vector x_t0)
                : Sigma_MIN(Sigma_Min), Sigma_MAX(Sigma_Max), Change_MAX(Change_Max), x_t0(x_t0) {
            Sigma_t0 = sqrt(Sigma_Min * Sigma_Max);
        }


        vector getChanges(vector &residuals, vector &x, bool isAdaptiveDamping) {
            assert(x.rows() == x_t0.rows() && "Damping hasn't been properly initialized");
            vector changes = x - x_t0;
            pr_t sigma{1.0};
            if (isAdaptiveDamping) {
                sigma = applyAdaptiveDamping(changes.maxCoeff(), getDampNorm(residuals, x));
                LOG(numerics) << "Damping applied: " << sigma;
            }
            x_t0 = x;
            return changes * sigma;
        }

        pr_t getNorm() { return norm_t0; }

    private:
        bool first{true};

        vector x_t0;
        pr_t norm_t0;
        pr_t max_change_t0;
        pr_t Sigma_t0;

        pr_t Sigma_MIN;
        pr_t Sigma_MAX;
        pr_t Change_MAX;


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


        pr_t applyAdaptiveDamping(pr_t max_change, pr_t norm) {
            pr_t PHI{0.01};
            pr_t sigma{0};

            pr_t p_n;
            pr_t p_h;

            if (first) {
                p_n = norm;
                p_h = max_change;
		first = false;
            } else {
                p_n = norm / norm_t0;
                p_h = max_change / max_change_t0;
            }

            if (p_n < 1 and p_h < 1) {
                auto lambda = log10(p_n) / log10(PHI);
                if (lambda < 1) {
                    sigma = Sigma_t0 + lambda * (Sigma_MAX - Sigma_t0);
                } else {
                    sigma = Sigma_MAX;
                }
                cnt = 0;
            }

            if (p_n > 1) {
                sigma = Sigma_t0 / p_n;
            }

            if (p_h > 1) {
                sigma = Sigma_t0 / p_h;
            }

            pr_t sig_t{sqrt(sigma * Sigma_t0)};

            if (std::abs(max_change) > Change_MAX and sig_t > (Change_MAX / std::abs(max_change))) {
                sig_t = Change_MAX / std::abs(max_change);
            }
            if (sig_t < Sigma_MIN) {
                sig_t = Sigma_MIN;
                cnt++;
                if (cnt > 10) {
                    sig_t = pow(pow(Sigma_MIN, 2) * Sigma_MAX, 1 / 3);
                }
            }

            Sigma_t0 = sig_t;
            norm_t0 = norm;
            max_change_t0 = max_change;
            return sig_t;
        }

};

}
}//ns
#endif //GLOBAL_FLOW_NUMERICS_HPP
