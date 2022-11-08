//
// Created by dk on 23.02.22.
//

#ifndef VARIABLE_DENSITY_FLOW_HPP
#define VARIABLE_DENSITY_FLOW_HPP

#include "Units.hpp"
#include <unordered_map>
#include "../Simulation/Options.hpp"
#include "../Misc/Helpers.hpp"

using namespace std;

namespace GlobalFlow {
    namespace Model {
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

                for (int id = 0; id < numberOfDensityZones; id++) {

                    if (densityStratified) {
                        nusZoneVec.push_back((densityZones[id] - densityFresh) / densityFresh * Model::si::si_dimensionless); // nus of zones is equal to nus of zeta surface below
                    } else { // if continuous todo: test for continuous case whether the results are correct
                        nusZoneVec.push_back((nusZetaVec[id] + nusZetaVec[id+1]) * 0.5 * Model::si::si_dimensionless); // nus of zones is mean of zeta surfaces at top and bottom
                    }
                    if (id == 0) {
                        delnusVec.push_back(nusZoneVec[id]); // density difference in top zone
                    } else {
                        delnusVec.push_back((nusZoneVec[id] - nusZoneVec[id-1]));
                    }
                    //LOG(debug) << "delnus (for zone " << id << "): " << delnusVec[id].value() << std::endl;
                }
                // todo sort nusZetas and nusZones ascending or throw error if not ascending
                densityProps.nusZetas = nusZetaVec;
                densityProps.nusZones = nusZoneVec;
                densityProps.delnus = delnusVec;
                densityProps.maxToeSlope = maxToeSlope * Model::si::si_dimensionless;
                densityProps.maxTipSlope = maxTipSlope * Model::si::si_dimensionless;
                densityProps.numOfZones = numberOfDensityZones;
                return densityProps;
            }

            vector<t_dim> getNusZones(){return nusZones;}

            vector<t_dim> getDelnus(){return delnus;}

            bool isDensityVariable(){return densityVariable;}

            t_dim getMaxToeSlope(){return maxToeSlope;}

            t_dim getMaxTipSlope(){return maxTipSlope;}

            int getNumOfZones(){return numOfZones;}
        };
    }
}
#endif //VARIABLE_DENSITY_FLOW_HPP