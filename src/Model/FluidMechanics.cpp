#include "FluidMechanics.hpp"
#include "../Misc/Helpers.hpp"

namespace GlobalFlow {
    namespace Model {

        t_meter FluidMechanics::calcDeltaV(t_meter head, t_meter elevation, t_meter verticalSize) noexcept {
            t_meter epsilon{1e-4 * si::meter};
            t_meter bottom = elevation - verticalSize;
            if (head >= elevation) { return verticalSize; }
            if (elevation > head and head + epsilon > bottom) { return head - bottom; }
            if (head + epsilon <= bottom) {
                return 0.0 * si::meter;
            }
            return 0.0 * si::meter;
        };

        quantity<MeterSquaredPerTime> calcEfoldingTrans(t_vel
                                                        k,
                                                        t_meter f, t_meter
                                                        z,
                                                        t_meter h
        ) {
            quantity<MeterSquaredPerTime> t; //transmissivity
            const quantity<Meter> d{100 * si::meter};
            if (f == 0 * si::meter) { f = 1 * si::meter; }
            if (h >= (z - d)) {
                //head is above d0
                t = k * f;
            } else {
                //head is below d0
                t_dim fold = exp(-((z - h - d) / f));
                if (fold.value() < 1e-50) { fold = 1e-50; }
                NANChecker(fold.value(), "E-folding problem");
                t = k * f * fold;
            }
            NANChecker(t.value(), "E-folding based transmissivity");
            return t;
        }

        quantity<MeterSquaredPerTime> FluidMechanics::calculateEFoldingConductance(FlowInputHor const& flow,
                                                                                   t_meter folding_self,
                                                                                   t_meter folding_neig) {
            t_vel k_neig;
            t_vel k_self;
            t_meter nodeLength_neig; // node length in flow direction (of neighbour)
            t_meter nodeLength_self; // node length in flow direction (of this node)
            t_meter nodeWidth_self; // node length perpendicular to flow direction (of this node)
            t_meter head_neig;
            t_meter head_self;
            t_meter ele_neig;
            t_meter ele_self;
            t_meter deltaV_neig;
            t_meter deltaV_self;
            bool confined;
            std::tie(k_neig, k_self, nodeLength_neig, nodeLength_self, nodeWidth_self, head_neig,
                     head_self, ele_neig, ele_self,deltaV_neig, deltaV_self, confined) = flow;
            quantity<MeterSquaredPerTime> out = 0 * si::square_meter / day;

            quantity<MeterSquaredPerTime> transmissivity_self =
                    calcEfoldingTrans(k_self, folding_self, ele_self, head_self);
            quantity<MeterSquaredPerTime> transmissivity_neig =
                    calcEfoldingTrans(k_neig, folding_neig, ele_neig, head_neig);

            if (transmissivity_neig != 0 * si::square_meter / day and
                transmissivity_self != 0 * si::square_meter / day) {
                out = (2.0 * nodeWidth_self) * ((transmissivity_self * transmissivity_neig)
                                                / (transmissivity_self * nodeLength_neig +
                                                   transmissivity_neig * nodeLength_self));
            }
            NANChecker(out.value(), "E-folding based Harmonic Mean Conductance");
            return out;
        }

        quantity<MeterSquaredPerTime> FluidMechanics::calculateHarmonicMeanConductance(FlowInputHor const& flow) noexcept {
            t_vel k_neig;
            t_vel k_self;
            t_meter nodeLength_neig; // node length in flow direction (of neighbour)
            t_meter nodeLength_self; // node length in flow direction (of this node)
            t_meter nodeWidth; // node length perpendicular to flow direction (of smaller node of the two interacting)
            t_meter head_neig;
            t_meter head_self;
            t_meter ele_neig;
            t_meter ele_self;
            t_meter deltaV_neig;
            t_meter deltaV_self;
            bool confined;
            std::tie(k_neig, k_self, nodeLength_neig, nodeLength_self, nodeWidth,head_neig,
                     head_self, ele_neig, ele_self,deltaV_neig, deltaV_self, confined) = flow;

            quantity<MeterSquaredPerTime> out = 0 * si::square_meter / day;
            //LOG(debug) << "nodeLength_self = " << nodeLength_self.value() << ", nodeLength_neig = " << nodeLength_neig.value();
            /*
            //Used if non dry-out approach is used
            // FIXME need to be checked here
            t_meter threshold_saturated_thickness{1e-5 * si::meter};

            //Cell is not confined layer -> we need to calculate transmissivity dependent on head
            enum e_dry {
                NONE, SELF, NEIG
            };
            e_dry dry{NONE};
            if (not confined) {
                deltaV_self = calcDeltaV(head_self, ele_self, deltaV_self);
                deltaV_neig = calcDeltaV(head_neig, ele_neig, deltaV_neig);
                if (deltaV_self == 0 * si::meter) {
                    deltaV_self = threshold_saturated_thickness;
                    dry = SELF;
                }
                if (deltaV_neig == 0 * si::meter) {
                    deltaV_neig = threshold_saturated_thickness;
                    dry = NEIG;
                }
                //TODO One of the cells is "dry" - Upstream weighting instead of harmonic mean
            }*/

            //Transmissivity = K * thickness of prism
            quantity<MeterSquaredPerTime> transmissivity_self = deltaV_self * k_self;
            quantity<MeterSquaredPerTime> transmissivity_neig = deltaV_neig * k_neig;

            if (transmissivity_neig != 0 * si::square_meter / day and
                transmissivity_self != 0 * si::square_meter / day) {
                out = nodeWidth * ((transmissivity_self * transmissivity_neig)
                        / (transmissivity_self * nodeLength_neig * 0.5 + // half of neighbour node's length
                        transmissivity_neig * nodeLength_self * 0.5)); // half of this node's length
            }
            NANChecker(out.value(), "Harmonic Mean Conductance");
            return out;
        };


        quantity<MeterSquaredPerTime>
        FluidMechanics::calculateVerticalConductance(FlowInputVert const& flowInputVer) noexcept {
            t_vel k_vert_neig;
            t_vel k_vert_self;
            t_meter verticalSize_self;
            t_meter verticalSize_neig;
            t_meter head_self;
            t_meter head_neig;
            t_meter elevation_self;
            t_meter elevation_neig;
            t_s_meter area_self;
            bool confined;

            std::tie(k_vert_neig,
                     k_vert_self,
                     verticalSize_self,
                     verticalSize_neig,
                     head_self,
                     head_neig,
                     elevation_self,
                     elevation_neig,
                     area_self,
                     confined) = flowInputVer;

            quantity<MeterSquaredPerTime> out = 0.0 * si::square_meter / day;
            t_meter deltaV_self = verticalSize_self;
            t_meter deltaV_neig = verticalSize_neig;

            //Used if non dry-out approach is used
            // FIXME need to be checked here
            t_meter threshhold_saturated_thickness{1e-5 * si::meter};

            enum e_dry {
                NONE, SELF, NEIG
            };
            e_dry dry{NONE};
            if (not confined) {
                deltaV_self = calcDeltaV(head_self, elevation_self, verticalSize_self);
                deltaV_neig = calcDeltaV(head_neig, elevation_neig, verticalSize_neig);
                if (deltaV_self == 0 * si::meter) {
                    deltaV_self = threshhold_saturated_thickness;
                    dry = SELF;
                }
                if (deltaV_neig == 0 * si::meter) {
                    deltaV_neig = threshhold_saturated_thickness;
                    dry = NEIG;
                }
                //TODO One of the cells is "dry" - Upstream weighting instead of harmonic mean
            }

            if (deltaV_self != 0.0 * si::meter and deltaV_neig != 0.0 * si::meter) {
                out = area_self / (((deltaV_self * 0.5) / k_vert_self) + ((deltaV_neig * 0.5) / k_vert_neig));
            }
            //LOG(debug) << "area_self: " << area_self.value() << ". deltaV_self: " << deltaV_self.value() << ". k_vert_self:" << k_vert_self.value();
            //LOG(debug) << "deltaV_neig: " << deltaV_neig.value() << ". k_vert_neig:" << k_vert_neig.value();

            NANChecker(out.value(), "Vertical Conductance");
            return out;
        };

        quantity<MeterSquaredPerTime> FluidMechanics::getHCOF(bool steadyState, quantity<Dimensionless> stepModifier,
                                                              t_s_meter storageCapacity,
                                                              quantity<MeterSquaredPerTime> P) noexcept {
            if (steadyState) {
                // LOG(debug) << "HCOF for steady state sim. (= P) (in getHCOF): " << P.value() << std::endl;
                return P;
            } else {
                quantity<MeterSquaredPerTime> out = P - (storageCapacity / (day * stepModifier) );
                //LOG(debug) << "P = " << P.value() << std::endl;
                //LOG(debug) << "storageCapacity = " << storageCapacity.value() << std::endl;
                //LOG(debug) << "HCOF for transient sim. (= P - storage capacity/time step): " << out.value() << std::endl;
                NANChecker(out.value(), "HCOF");
                return out;}
        }
    }
}//ns
