//
// Created by dk on 23.02.22.
//

#ifndef VARIABLEDENSITYFLOW_HPP
#define VARIABLEDENSITYFLOW_HPP

#include "Units.hpp"

using namespace std;

namespace GlobalFlow {
    namespace Model {
        class DensityProperties{
        protected:
            vector<t_dim> nusZeta; // dimensionless density on the zeta surfaces
            vector<t_dim> nusZone; // dimensionless density in the density zones between successive zeta surfaces
            vector<t_dim> delnus; // difference in dimensionless density between successive zeta surfaces
            vector<t_dim> eps; // variation of dimensionless density over a density zone

        private:


            /*// Question: what is this good for?
            friend class boost::serialization::access;
            template<class Archive>
            void serialize(Archive & ar, const unsigned int version) {
                LOG(debug) << "Serializing node-independent density properties";
                ar & delnus;
                ar & eps;
                ar & nusZeta;
                ar & nusZone;
                }*/
        public:


            static DensityProperties setDensityProperties(bool densityStratified, double densityFresh,
                                                          vector<double> densityZones,
                                                          int numberOfDensityZones){
                DensityProperties densityProps;

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
                densityProps.nusZeta = nusZetaVec;
                densityProps.nusZone = nusZoneVec;
                densityProps.delnus = delnusVec;
                densityProps.eps = epsVec;
                return densityProps;
            }

            /*void addNusZone(int id, t_dim value){
                std::pair<int,t_dim> newEntry (id, value);
                // todo: what to do if there is already a value at that id?
                nusZone.insert (newEntry); // Question: or zetas[id] = zetaHeight; ?
            }

            void addDelnus(int id, t_dim value){
                std::pair<int,t_dim> newEntry (id, value);
                // todo: what to do if there is already a value at that id?
                delnus.insert (newEntry); // Question: or zetas[id] = zetaHeight; ?
            }

            void addEps(int id, t_dim value){
                std::pair<int,t_dim> newEntry (id, value);
                // todo: what to do if there is already a value at that id?
                eps.insert (newEntry); // Question: or zetas[id] = zetaHeight; ?
            }*/

            /* todo remove density properties
            void removeNusZeta(int id){
                set<vector<t_dim>, NusZeta>(nusZeta);
            }

            void removeNusZone(int id){
                set<vector<t_dim>, NusZeta>(nusZone);
            }

            void removeDelnus(int id){
                set<vector<t_dim>, Delnus>(delnus);
            }

            void removeEps(int id){
                set<vector<t_dim>, Eps>(eps);
            }*/
        };
    }
}
#endif //VARIABLEDENSITYFLOW_HPP
