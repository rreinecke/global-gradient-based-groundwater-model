//
// Created by dk on 23.02.22.
//

#ifndef TESTING_DENSITYSURFACES_HPP
#define TESTING_DENSITYSURFACES_HPP

#include "Units.hpp"
#include "../Misc/Helpers.hpp"

namespace GlobalFlow {
    namespace Model {

        // int maxNumberZetas = 100; // max number of Zetas (Density Surfaces)

        std::list<int> ZetaN{}; // Empty list of Zeta IDs

        struct ZetaHash {
            template<typename T>
            std::size_t operator()(T t) const {
                return static_cast<std::size_t>(t);
            }
        };

        /**
         * @class Zeta
         *
         */
        class Zeta{
        public:
            /**
             * @brief Constructor for Density Surface "Zeta"
             * @param n : density zone identifier n (n=1 at top of aquifer, n=numberDensityZones+1 at bottom of aquifer)
             * @param elevation : elevation of density surface n (located at the top of the density zone n)
             * @param nus : difference in the dimensionless density between density surface n (above) & n+1 (below)
             */
            Zeta(int n, t_meter elevation, t_dim nus) : n(n), elevation(elevation), nus(nus) {};

            // Density surface properties:
            // TODO: t_meter getZetaPosition(int zetaId) {} // returns position relative to top and bottom of the cell (in SWI2: zetaPosition = 0: between, 1: at top, 2: at bottom, 3: cell inactive)
            // TODO: t_dim getToeAngle(int zetaId)
            // TODO: t_dim getTipAngle(int zetaId)

            int getZetaN() const noexcept { return n; }

            t_meter getElevation() const noexcept { return elevation; }

            t_dim getNus() const noexcept { return nus; }


        private:
            const int n;
            const t_meter elevation;
            const t_dim nus;
        };


    }
}

#endif //TESTING_DENSITYSURFACES_HPP
