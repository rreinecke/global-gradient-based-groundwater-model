//
// Created by dk on 23.02.22.
//

#ifndef TESTING_VARIABLEDENSITYFLOW_HPP
#define TESTING_VARIABLEDENSITYFLOW_HPP

#include "Units.hpp"
#include "../Misc/Helpers.hpp"

namespace GlobalFlow {
    namespace Model {



        // Density zone properties:
        t_dim getNusZone(int nZone, bool stratified){ // returns dimensionless density of the zone n (which reaches from surface n up to surface n-1, or if n=1 up to top of aquifer)

        }

        t_dim getDelnus(int zetaN, bool stratified) { // returns the difference in dimensionless density between surface zetaIdTop and zetaIdTop-1
            if (zetaIdTop == 1) {
                return getNusZone(zetaN, stratified);
            } else {
                return getNusZone(zetaN, stratified) -getNusZone(zetaN - 1, stratified);
            }
        }

        t_dim getEps(int zetaN, bool stratified){
            if (stratified) {
                return 0;
            } else { // if continuous
                return (getNus(zetaN + 1, true) - getNus(zetaN, true)) / 6;
            }
        }

        // TODO: t_dim getSourceType() {} // returns the type of source (in SWI2:
        // sourceTypeZone = >0: sources and sinks are of the same type as water in this sourceTypeZone
        //                    =0: sources and sinks are of the same type of water as the water at the top of the aquifer
        //                    <0: (for SGD!) sources are of the same type as the water in this sourceTypeZone, sinks are the same type as water at the top of the aquifer


        // Density surface properties:
        // TODO: t_meter getZetaPosition(int zetaId) {} // returns position relative to top and bottom of the cell (in SWI2: zetaPosition = 0: between, 1: at top, 2: at bottom, 3: cell inactive)
        // TODO: t_dim getToeAngle(int zetaId)
        // TODO: t_dim getTipAngle(int zetaId)


        /**
         * @param n : density zone identifier n (n=1 at top of aquifer, n=numberDensityZones+1 at bottom of aquifer)
         * @param elevation : elevation of density surface n (located at the top of the density zone n)
         * @param nus : difference in the dimensionless density between density surface n (above) & n+1 (below)
         */






    }
}

#endif //TESTING_VARIABLEDENSITYFLOW_HPP
