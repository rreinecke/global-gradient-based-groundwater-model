#include "VariableDensityFlow.hpp"

namespace GlobalFlow {
    namespace Model {

        std::vector<t_meter> VariableDensityFlow::calculateZoneThicknesses(
                std::unordered_map<t_dim, t_meter> zetas,
                std::unordered_map<t_dim, t_meter> zetas_neig,
                t_meter edgeLength_neig,
                t_meter edgeLength_self) noexcept {
            // Question: how to adapt this or the definition of zetas so we can track across multiple layers?
            // calculate zone thicknesses
            std::vector<t_meter> out;
            for (int p; p <= zetas.size() - 1; p++) {
                out.push_back((edgeLength_neig * (zetas[p] - zetas[p + 1])) / (edgeLength_neig + edgeLength_self)) +
                ((edgeLength_self * (zetas_neig[p] - zetas_neig[p + 1])) / (edgeLength_neig + edgeLength_self));
            } // todo: check if this works correctly (e.g by summing up all zone thicknesses. Expected result: aquifer thickness
            return out;
        }

        std::vector<t_s_meter_t> VariableDensityFlow::calculateDensityZoneConductances(
                std::vector<t_meter> zoneThicknesses, t_s_meter_t conductance)noexcept {
            // calculate the sum of the zone thicknesses
            t_meter sumOfZoneThicknesses = 0 * si::meter;
            std::for_each(zoneThicknesses.begin(), zoneThicknesses.end(), [&](t_meter zoneThickness) {
                sumOfZoneThicknesses += zoneThickness;
            }); // todo: check if this works correctly
            // calculate the density zone conductances
            std::vector<t_s_meter> out;
            for (int n; n <= zoneThicknesses.size(); n++) {
                out.push_back(conductance * (zoneThicknesses[n] / sumOfZoneThicknesses));
            } // todo: check if this works correctly
            return out;
        }

        std::vector<t_s_meter_t> VariableDensityFlow::calculateCumulativeDensityZoneConductances(
                std::vector<t_s_meter_t> densityZoneCond)noexcept {
            // calculate the sum of density zone conductances below a zeta surface n and add to vector out
            std::vector<t_s_meter_t> out;
            t_s_meter_t sumDensityZoneCondBelow = 0 * si::square_meter / day;
            for (int n; n <= densityZoneCond.size(); n++) {
                std::for_each(densityZoneCond.begin() + n, densityZoneCond.end(),
                              [&](t_s_meter_t densityZoneCondBelow) {
                    sumDensityZoneCondBelow += densityZoneCondBelow;
                });
                out.push_back(sumDensityZoneCondBelow);
                sumDensityZoneCondBelow = 0;
            } // todo: check if this works correctly
            return out;
        }

    }
}
