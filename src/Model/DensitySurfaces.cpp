#include "DensitySurfaces.hpp"

namespace GlobalFlow {
    namespace Model {


        // Density zone properties:
        t_dim getNusZone(int zetaN, bool stratified){ // returns dimensionless density of the zone n (which reaches from surface n up to surface n-1, or if n=1 up to top of aquifer)
            if (stratified) {
                return this->nus;
            } else { // if continuous
                return 0.5 * (getNusZone(zetaN, true) + getNusZone(zetaN + 1, true));
            }
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

