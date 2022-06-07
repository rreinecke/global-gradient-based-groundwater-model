//
// Created by dk on 23.02.22.
//

#ifndef VARIABLEDENSITYFLOW_HPP
#define VARIABLEDENSITYFLOW_HPP

#include "Units.hpp"
#include <unordered_map>

// #include "../Misc/Helpers.hpp"
using namespace std;

namespace GlobalFlow {
    namespace Model {
        class VariableDensityFlow{
        protected:
            unordered_map<int, t_dim> nusZeta; // dimensionless density on the zeta surfaces
            vector<t_dim> nusZone; // dimensionless density in the density zones between successive zeta surfaces
            vector<t_dim> delnus; // difference in dimensionless density between successive zeta surfaces
            vector<t_dim> eps; // variation of dimensionless density over a density zone
        public:
            friend class boost::serialization::access;
            template<class Archive>
            void serialize(Archive & ar, const unsigned int version) {
                LOG(debug) << "Serializing node-independent density properties";
                ar & delnus;
                ar & eps;
                ar & nusZeta;
                ar & nusZone;
            }

            // todo add density properties
            void addZoneProperties(bool densityStratified, double densityFresh, vector<double> densityZetas,
                                   int numberOfDensityZones){
                pair<int,t_dim> tmp;
                vector<double> densities = densityZetas;
                densities.insert(densities.begin(),densityFresh);
                densities.push_back(densityFresh);

                for (int id = 0; id <= numberOfDensityZones; id++) {
                    pair<int,t_dim> tmp (id, (densities[id] - densityFresh ) / densityFresh * Model::si::si_dimensionless);
                    nusZeta.insert(tmp);
                }

                for (int id = 0; id <= numberOfDensityZones-1; id++) {
                    if (densityStratified) {
                        pair<int,t_dim> tmp (id,  nusZeta[id+1] * Model::si::si_dimensionless);
                        nusZone.insert(tmp); // nus of zones is equal to nus of zeta surface below

                        pair<int,t_dim> tmp (id,  0 * Model::si::si_dimensionless);
                        eps.insert(tmp);
                    } else { // if continuous
                        pair<int,t_dim> tmp (id, 0.5 * ( nusZeta[id] + nusZeta[id+1] ) * Model::si::si_dimensionless);
                        nusZone.insert(tmp); // nus of zones is mean of zeta surfaces above and below

                        pair<int,t_dim> tmp (id, (( nusZeta[id+1] - nusZeta[id] ) / 6) * Model::si::si_dimensionless)
                        eps.insert(id, tmp);
                    }

                    if (id == 0) {
                        pair<int,t_dim> tmp (id, nusZone[id] * Model::si::si_dimensionless)
                        delnus.insert(tmp); // by density difference in top zone
                    } else {
                        pair<int,t_dim> tmp (id, (nusZone[id] - nusZone[id-1]) * Model::si::si_dimensionless);
                        delnus.insert(tmp);
                    }
                }
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
