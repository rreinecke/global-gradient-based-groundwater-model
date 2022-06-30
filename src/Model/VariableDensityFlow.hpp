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
            vector<t_dim> nusZetas; // dimensionless density on the zeta surfaces
            vector<t_dim> nusZones; // dimensionless density in the density zones between successive zeta surfaces
            vector<t_dim> delnus; // difference in dimensionless density between successive zeta surfaces
            vector<t_dim> eps; // variation of dimensionless density over a density zone

        /*private:
            // Question: what is this good for?
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
                // todo sort nusZetas and nusZones ascending or throw error if not ascending
                densityProps.nusZetas = nusZetaVec;
                densityProps.nusZones = nusZoneVec;
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

        // todo: add THICKRF, THICKFF (lineraly interpolated zone thickness at the interface between columns/rows)
        /**
         * Calculates the density zone thicknesses between the density surfaces
         * @param zetas zeta surface height (zeta[p+1]: surface below)
         * @param zetas_neig zeta surface height in neighbouring cell (zeta_neig[p+1]: surface below)
         * @param edgeLength_neig edge width/length of neighbouring column/row
         * @param edgeLength_self edge width/length of this column/row
         * @return zone thickness
         */
        //
        std::vector<quantity<Meter>> VariableDensityFlow::calculateZoneThicknesses(
                std::unordered_map<t_dim, t_meter> zetas,
                std::unordered_map<t_dim, t_meter> zetas_neig,
                t_meter edgeLength_neig,
                t_meter edgeLength_self)noexcept {
            // Question: how to adapt this or the definition of zetas so we can track across multiple layers?
            // calculate zone thicknesses
            std::vector<quantity<Meter>> out;
            for (int p; p <= zetas.size() - 1; p++){
                out.push_back((edgeLength_neig*(zetas[p] - zetas[p+1]))/(edgeLength_neig+ edgeLength_self)) +
                ((edgeLength_self*(zetas_neig[p] - zetas_neig[p+1]))/(edgeLength_neig + edgeLength_self));
            } // todo: check if this works correctly (e.g by summing up all zone thicknesses. Expected result: aquifer thickness
            return out;
        }

        // todo: add SWICR, SWICC (conductance in the column/row direction for a zone)
        std::vector<quantity<MeterSquaredPerTime>> VariableDensityFlow::calculateDensityZoneConductances(std::vector<t_meter> zoneThickness, t_s_meter_t conductance)noexcept {
        // calculate the sum of the zone thicknesses
        quantity<Meter> sumOfZoneThicknesses = 0 * si::meter;
        std::for_each(zoneThicknesses.begin(), zoneThicknesses.end(),[&] (t_meter zoneThickness) {
        sumOfZoneThicknesses += zoneThickness;
        }); // todo: check if this works correctly
        // calculate the density zone conductances
        std::vector<quantity<MeterSquaredPerTime>> out;
        for (int n; n <= zoneThicknesses.size(); n++){
        out.push_back(conductance * (zoneThicknesses[n] / sumOfZoneThicknesses));
        } // todo: check if this works correctly
        return out;
        }

        // todo: add SWICUMCR, SWICUMCC (cumulative conductance in column/row direction below a density surface n)
        std::vector<quantity<MeterSquaredPerTime>> VariableDensityFlow::calculateCumulativeDensityZoneConductances(std::vector<t_s_meter_t> densityZoneCond, conductance)noexcept {
        // calculate the sum of density zone conductances below a zeta surface n and add to vector out
        std::vector<quantity<MeterSquaredPerTime>> out;
        quantity<MeterSquaredPerTime> sumDensityZoneCondBelow = 0 * si::square_meter / day;
        for (int n; n <= densityZoneCond.size(); n++){
        std::for_each(densityZoneCond.begin()+n, densityZoneCond.end(), [&] (t_s_meter_t densityZoneCondBelow) {
        sumDensityZoneCondBelow += densityZoneCondBelow;
        });
        out.push_back(sumDensityZoneCondBelow);
        sumDensityZoneCondBelow = 0;
        } // todo: check if this works correctly
        return out;
        }
    }
}
#endif //VARIABLEDENSITYFLOW_HPP
