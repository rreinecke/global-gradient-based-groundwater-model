#include "VariableDensityFlow.hpp"
#include "../Misc/Helpers.hpp"

namespace GlobalFlow {
    namespace Model {

        std::vector<t_meter>
        VariableDensityFlow::calculateZoneThicknesses(
                std::vector<t_meter> zetas, // std::unordered_map<t_dim, t_meter>
                std::vector<t_meter> zetas_neig, // std::unordered_map<t_dim, t_meter>
                t_meter edgeLength_neig,
                t_meter edgeLength_self) noexcept {
            // Question: how to adapt this or the definition of zetas so we can track across multiple layers?
            // calculate zone thicknesses
            std::vector<t_meter> out;
            t_meter zoneThickness;
            int numOfZones = options.getNumberOfDensityZones();
            for (int p = 0; p < numOfZones; p++) {
                zoneThickness = ((edgeLength_neig * (zetas[p] - zetas[p + 1])) / (edgeLength_neig + edgeLength_self)) + ((edgeLength_self * (zetas_neig[p] - zetas_neig[p + 1])) / (edgeLength_neig + edgeLength_self));
                out.push_back(zoneThickness);
            }
            return out;
        }

        std::vector<t_s_meter_t>
        VariableDensityFlow::calculateDensityZoneConductances(
                std::vector<t_meter> zoneThicknesses, t_s_meter_t conductance)noexcept {
            // calculate the sum of the zone thicknesses
            t_meter sumOfZoneThicknesses = 0 * si::meter;
            std::for_each(zoneThicknesses.begin(), zoneThicknesses.end(), [&](t_meter zoneThickness) {
                sumOfZoneThicknesses += zoneThickness;
            });
            // calculate the density zone conductances
            std::vector<t_s_meter_t> out;
            t_s_meter_t densityZoneConductance;
            int numOfZones = options.getNumberOfDensityZones();
            for (int n = 0; n < numOfZones; n++) {
                densityZoneConductance = conductance * (zoneThicknesses[n] / sumOfZoneThicknesses);
                out.push_back(densityZoneConductance);
            }
            return out;
        }

        std::vector<t_s_meter_t>
        VariableDensityFlow::calculateCumulativeDensityZoneConductances(
                std::vector<t_s_meter_t> densityZoneCond)noexcept {
            // calculate the sum of density zone conductances below a zeta surface n and add to vector out
            std::vector<t_s_meter_t> out;
            t_s_meter_t sumDensityZoneCondBelow = 0 * si::square_meter / day;
            int numOfZones = options.getNumberOfDensityZones();
            for (int n = 0; n < numOfZones; n++) {
                std::for_each(densityZoneCond.begin() + n, densityZoneCond.end(), [&](t_s_meter_t densityZoneCondBelow) {
                    sumDensityZoneCondBelow += densityZoneCondBelow;
                });
                out.push_back(sumDensityZoneCondBelow);
                sumDensityZoneCondBelow = 0;
            }
            return out;
        }

        /*std::vector<t_s_meter_t>
        VariableDensityFlow::calculateZetaMovementConductances(
                std::vector<t_s_meter_t> densityZoneCond, std::vector<t_s_meter_t> densityZoneCondCum,
                std::vector<t_dim> delnus, std::vector<t_dim> eps) noexcept {

            std::vector<t_s_meter_t> out;
            t_s_meter_t tmp;
            int numOfZones = options.getNumberOfDensityZones();
            for (int n = 0; n < numOfZones; n++) {
                tmp = (delnus[n] * densityZoneCondCum[n]) - (eps[n] * densityZoneCond[n]);
                out.push_back(tmp);
                NANChecker(out[n].value(), "zetaMovementConductance");
            }

            return out;
        }*/


    }
}