//
// Created by dk on 23.02.22.
//

#ifndef VARIABLE_DENSITY_FLOW_HPP
#define VARIABLE_DENSITY_FLOW_HPP

#include "Units.hpp"
#include <unordered_map>
#include "../Simulation/Options.hpp"

using namespace std;

namespace GlobalFlow {
    namespace Model {

        /**
         * @class VariableDensityFlow
         * Provides helper functions for variable density calulcations
         */
        class VariableDensityFlow {
        private:
            Simulation::Options options;
        public:
            VariableDensityFlow(){}

            /**
             * @brief calculates the thicknesses of the density zones between current and neighbouring node by linear interpolation (in SWI2: THICKRF, THICKFF)
             * @param zetas zeta surface height (zeta[p+1]: surface below)
             * @param zetas_neig zeta surface height in neighbouring cell (zeta_neig[p+1]: surface below)
             * @param edgeLength_neig edge width/length of neighbouring column/row
             * @param edgeLength_self edge width/length of this column/row
             * @return zone thickness
             */
            std::vector<t_meter> calculateZoneThicknesses(std::vector<t_meter> zetas,
                                                          std::vector<t_meter> zetas_neig,
                                                          t_meter edgeLength_neig,
                                                          t_meter edgeLength_self) noexcept;

            /**
             * @brief calculates the conductance in the column/row direction for a density zone (in SWI2: SWICR & SWICC)
             * @param zoneThicknesses thicknesses of zones
             * @param conductance conductance of cell towards respective neighbouring cell
             * @return vector with conductances of density zones
             */
            std::vector<t_s_meter_t> calculateDensityZoneConductances(
                    std::vector<t_meter> zoneThicknesses, t_s_meter_t conductance) noexcept;

            /**
             * @brief calculates the cumulative conductance in column/row direction below a density surface n (in SWI2: SWICUMCR & SWICUMCC)
             * @param densityZoneCond
             * @param conductance conductance of cell towards respective neighbouring cell
             * @return vector with conductances of density zones
             */
            std::vector<t_s_meter_t> calculateCumulativeDensityZoneConductances(
                    std::vector<t_s_meter_t> densityZoneCond) noexcept;


            /**
             * @brief calculates the horizontal conductances used to solve ZETA surface heights (in SWI2: SWISOLCC/R)
             * @param densityZoneCond vector of density zone condutctances
             * @param densityZoneCondCum vector of density zone condutctances below each surface
             * @param delnus difference in dimensionless density between successive zeta surfaces
             * @param eps variation of dimensionless density over a density zone (for continuous option)
             * @return vector with horizontal conductances used to solve ZETA surface heights
             */
            std::vector<t_s_meter_t> calculateZetaMovementConductances(
                    std::vector<t_s_meter_t> densityZoneCond,
                    std::vector<t_s_meter_t> densityZoneCondCum,
                    std::vector<t_dim> delnus,
                    std::vector<t_dim> eps) noexcept;

        };

        /**
         * @class DensityProperties
         * Class to store the dimensionless density of zones and surfaces, & info about density change between surfaces
         */
        class DensityProperties{
        protected:
            bool densityVariable;
            vector<t_dim> nusZetas; // dimensionless density on the zeta surfaces
            vector<t_dim> nusZones; // dimensionless density in the density zones between successive zeta surfaces
            vector<t_dim> delnus; // difference in dimensionless density between successive zeta surfaces
            vector<t_dim> eps; // variation of dimensionless density over a density zone (for continuous option)
            t_dim maxToeSlope;
            t_dim maxTipSlope;
            int numOfZones;

        public:


            static DensityProperties setDensityProperties(bool densityVariable, bool densityStratified, double densityFresh,
                                                          vector<double> densityZones, int numberOfDensityZones,
                                                          double maxToeSlope, double maxTipSlope){
                DensityProperties densityProps;
                densityProps.densityVariable = densityVariable;
                vector<t_dim> nusZetaVec;
                vector<t_dim> nusZoneVec;
                vector<t_dim> delnusVec;
                vector<t_dim> epsVec;
                double densityZeta;

                for (int id = 0; id <= numberOfDensityZones; id++) {
                    if (id == 0) {
                        densityZeta = densityZones[id];
                    } else if (id == numberOfDensityZones) {
                        densityZeta = densityZones[id-1];
                    } else {
                        densityZeta = (densityZones[id-1] + densityZones[id]) * 0.5;
                    }
                    nusZetaVec.push_back((densityZeta - densityFresh ) / densityFresh * Model::si::si_dimensionless);
                }

                for (int id = 0; id <= numberOfDensityZones-1; id++) {

                    if (densityStratified) {
                        nusZoneVec.push_back((densityZones[id] - densityFresh) / densityFresh * Model::si::si_dimensionless); // nus of zones is equal to nus of zeta surface below
                        epsVec.push_back(0 * Model::si::si_dimensionless);
                    } else { // if continuous todo: test for continuous case whether the results are correct
                        nusZoneVec.push_back(((densityZones[id] - densityFresh) * 0.5 +
                                              (densityZones[id+1] - densityFresh) * 0.5)
                                             * 0.5 * Model::si::si_dimensionless); // nus of zones is mean of zeta surfaces above and below
                        epsVec.push_back((( nusZetaVec[id+1] - nusZetaVec[id] ) / 6));
                    }

                    if (id == 0) {
                        delnusVec.push_back(nusZoneVec[id]); // density difference in top zone
                    } else {
                        delnusVec.push_back((nusZoneVec[id] - nusZoneVec[id-1]));
                    }
                }
                // todo sort nusZetas and nusZones ascending or throw error if not ascending
                densityProps.nusZetas = nusZetaVec;
                densityProps.nusZones = nusZoneVec;
                densityProps.delnus = delnusVec;
                densityProps.eps = epsVec;
                densityProps.maxToeSlope = maxToeSlope * Model::si::si_dimensionless;
                densityProps.maxTipSlope = maxTipSlope * Model::si::si_dimensionless;
                densityProps.numOfZones = numberOfDensityZones;
                return densityProps;
            }

            vector<t_dim> getNusZones(){return nusZones;}

            vector<t_dim> getDelnus(){return delnus;}

            vector<t_dim> getEps(){return eps;}

            bool isDensityVariable(){return densityVariable;}

            t_dim getMaxToeSlope(){return maxToeSlope;}

            t_dim getMaxTipSlope(){return maxTipSlope;}
        };
    }
}
#endif //VARIABLE_DENSITY_FLOW_HPP