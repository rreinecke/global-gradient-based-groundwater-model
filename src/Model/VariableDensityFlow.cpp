/*
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
            t_meter deltaZeta;
            t_meter deltaZeta_neig;
            int numOfZones = options.getNumberOfDensityZones();

            for (int localZetaID = 0; localZetaID <= numOfZones; localZetaID++) {
                deltaZeta = zetas[localZetaID] - zetas[localZetaID + 1];
                deltaZeta_neig = zetas_neig[localZetaID] - zetas_neig[localZetaID + 1];
                if (deltaZeta <= (0 * si::meter) or deltaZeta_neig <= (0 * si::meter)){ // adapted from SWI2 code line 1149
                    zoneThickness = 0 * si::meter;
                } else {
                    zoneThickness = ((edgeLength_neig * deltaZeta) / (edgeLength_neig + edgeLength_self)) +
                                    ((edgeLength_self * deltaZeta_neig) / (edgeLength_neig + edgeLength_self));
                }
                //LOG(debug) << "zoneThickness (in calculateZoneThicknesses):" << zoneThickness.value() << std::endl;
                NANChecker(zoneThickness.value(), "zoneThickness");
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
            //LOG(debug) << "sumOfZoneThicknesses:" + std::to_string(sumOfZoneThicknesses.value());

            // calculate the density zone conductances
            std::vector<t_s_meter_t> out;
            t_s_meter_t densityZoneConductance;
            std::for_each(zoneThicknesses.begin(), zoneThicknesses.end(), [&](t_meter zoneThickness) {
                if (sumOfZoneThicknesses == (0 * si::meter)) { // adapted from SWI2 code line 1159
                    densityZoneConductance = 0 * si::square_meter / day;
                } else {
                    densityZoneConductance = conductance * (zoneThickness / sumOfZoneThicknesses);
                }
                //LOG(debug) << "densityZoneConductance (in calculateDensityZoneConductances):" << densityZoneConductance.value() << std::endl;
                out.push_back(densityZoneConductance);

                NANChecker(densityZoneConductance.value(), "densityZoneConductance");
            });
            return out;
        }

        t_s_meter_t
        VariableDensityFlow::calculateZoneConductanceCum(int n, std::vector<t_s_meter_t> densityZoneCond)noexcept {
            // calculate the sum of density zone conductances below a zeta surface n and add to vector out
            t_s_meter_t out = 0 * si::square_meter / day;

            std::for_each(densityZoneCond.begin() + n, densityZoneCond.end(), [&](t_s_meter_t densityZoneCondBelow) {
                out += densityZoneCondBelow;
            });
            //LOG(debug) << "zoneConductancesCum (in calculateZoneConductanceCum):" << out.value() << std::endl;
            NANChecker(out.value(), "calculateZoneConductanceCum");
            return out;
        }

        */
/*t_s_meter_t
        VariableDensityFlow::calculateZetaMovementConductance(int n,
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
        }*//*



    }
}*/