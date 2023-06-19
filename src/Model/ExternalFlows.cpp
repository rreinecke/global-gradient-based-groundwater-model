#include "ExternalFlows.hpp"

namespace GlobalFlow {
namespace Model {

t_s_meter_t ExternalFlow::getP(t_meter eq_head, t_meter head,
                               t_vol_t recharge,
                               t_vol_t eqFlow) const noexcept {
    t_s_meter_t out = 0.0 * (si::square_meter / day);
    switch (type) {
        case RECHARGE:
            return out;
        case FAST_SURFACE_RUNOFF:
            return out;
        case NET_ABSTRACTION:
            return out;
        case EVAPOTRANSPIRATION:
            //flowHead = surface, bottom = extinction depth
            if ((head < flowHead - bottom) xor (head > flowHead)) { return out;
            } else { return -special_flow / bottom; }
        case FLOODPLAIN_DRAIN:
            return out;
        case RIVER:
            return -conductance;
        case RIVER_MM:
            //Can happen in transient coupling
            if (flowHead <= bottom){
                //stil allow gaining conditions!
                if(head >= bottom){return -calcERC(recharge, eq_head, head, eqFlow);}
                return out;}
            return -calcERC(recharge, eq_head, head, eqFlow);

            /*if (flowHead <= bottom){
		    //stil allow gaining conditions!
                if(head >= bottom){ return -calcERC(recharge, eq_head, head, eqFlow);
                } else { return out; }
            } else { return -calcERC(recharge, eq_head, head, eqFlow); }*/
        case WETLAND:
            //Can happen in transient coupling
            if (flowHead <= bottom){
                if(head >= bottom){return -conductance;}
                return out; }
            return -conductance;

            /*if (flowHead <= bottom) {
		        if (head >= bottom) { return -conductance;
                } else { return out; }
            } else { return -conductance; }*/
        case GLOBAL_WETLAND:
            //Can happen in transient coupling
            if (flowHead <= bottom){
                if(head >= bottom){return -conductance;}
                return out; }
            return -conductance;

            /*if (flowHead <= bottom) {
		        if (head >= bottom) { return -conductance; }
                else { return out; }
            } else { return -conductance; }*/
        case LAKE:
            //Can happen in transient coupling
            if (flowHead <= bottom){
                if(head >= bottom){return -conductance;}
                return out; }
            return -conductance;

            /*if (flowHead <= bottom) {
		        if (head >= bottom) { return -conductance;
                } else { return out; }
            } else { return -conductance; }*/
        /*case GLOBAL_LAKE:
            //Can happen in transient coupling
            if (flowHead <= bottom) {
                if (head >= bottom) { return -conductance;
                } else { return out; }
            } else { return -conductance; }*/
        case DRAIN:
            if (head > flowHead) { return -calcERC(recharge, eq_head, head, eqFlow);
            } else { return out; }
        case GENERAL_HEAD_BOUNDARY:
            return -conductance;
    }
    return out;
}

t_vol_t ExternalFlow::getQ(t_meter eq_head, t_meter head,
                           t_vol_t recharge,
                           t_vol_t eqFlow) const noexcept {
    quantity<VolumePerTime, double> out = 0.0 * (si::cubic_meter / day);
    switch (type) {
        case RECHARGE:
            return this->special_flow;
        case FAST_SURFACE_RUNOFF:
            return this->special_flow;
        case NET_ABSTRACTION:
            return this->special_flow;
        case EVAPOTRANSPIRATION:
            if (head < flowHead - bottom) {
                return out;
            } else if (flowHead - bottom <= head and head <= flowHead) {
                return -this->special_flow + (this->special_flow * flowHead) / bottom;
            } else {
                return -this->special_flow;
            }
        case FLOODPLAIN_DRAIN:
            return -calculateFloodplainDrainage(head);
        case RIVER:
            return conductance * flowHead;
        case RIVER_MM:
            //Can happen in transient coupling
            if (flowHead <= bottom){
                if(head >= bottom){return calcERC(recharge, eq_head, head, eqFlow) * flowHead;}
                return out; }
            return calcERC(recharge, eq_head, head, eqFlow) * flowHead;

            /*if (flowHead <= bottom) {
                if (head >= bottom) { return calcERC(recharge, eq_head, head, eqFlow) * flowHead;
                } else { return out; }
            } else { return calcERC(recharge, eq_head, head, eqFlow) * flowHead; }*/
        case WETLAND:
            //Can happen in transient coupling
            if (flowHead <= bottom){
                if(head >= bottom){return conductance * flowHead;}
                return out; }
            return conductance * flowHead;

            /*if (flowHead <= bottom) {
		        if (head >= bottom) { return conductance * flowHead;
                } else { return out; }
            } else { return conductance * flowHead; }*/
        case GLOBAL_WETLAND:
            //Can happen in transient coupling
            if (flowHead <= bottom){
                if(head >= bottom){return conductance * flowHead;}
                return out; }
            return conductance * flowHead;

            /*if (flowHead <= bottom) {
		        if (head >= bottom) { return conductance * flowHead;
                } else { return out; }
            } else { return conductance * flowHead; }*/
        case LAKE:
            //Can happen in transient coupling
            if (flowHead <= bottom){
                if(head >= bottom){return conductance * flowHead;}
                return out; }
            return conductance * flowHead;

            /*if (flowHead <= bottom) {
		        if (head >= bottom) { return conductance * flowHead;
                } else { return out; }
            } else { return conductance * flowHead; }*/
        /*case GLOBAL_LAKE:
            //Can happen in transient coupling
            if (flowHead <= bottom) {
                if (head >= bottom) { return conductance * flowHead;
                } else { return out; }
            } else { return conductance * flowHead; }*/
        case DRAIN:
            if (head > flowHead) { return calcERC(recharge, eq_head, head, eqFlow) * flowHead;
            } else { return out; }
        case GENERAL_HEAD_BOUNDARY:
            return conductance * flowHead;
    }
    return out;
}

t_vol_t ExternalFlow::calculateFloodplainDrainage(t_meter head) const noexcept {
    quantity<VolumePerTime, double> out = 0.0 * (si::cubic_meter / day);
    t_meter headAboveFloodplain = head - flowHead;
    if (headAboveFloodplain > 0 * si::meter) {
        const double PI = std::atan(1.0) * 4;
        double J = (PI * conductance.value()) / (4 * 0.15 * (500 * 500));
        return bottom * bottom * headAboveFloodplain * (J * 1 / day);
    }
    return out;
}

t_meter ExternalFlow::getRiverDiff(t_meter eqHead) const noexcept {
    return eqHead - flowHead;
}

/*
 * From ATI pixel shaders
 */
double clamp(double x, double lowerlimit, double upperlimit) {
    if (x < lowerlimit) x = lowerlimit;
    if (x > upperlimit) x = upperlimit;
    return x;
}

double smoothstep(double edge0, double edge1, double x) {
    // Scale, bias and saturate x to 0..1 range
    x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    // Evaluate polynomial
    return x * x * (3 - 2 * x);
}

/**
 * @brief Equilibrium river conductance (ERC) following Miguez-Macho et al. (2007):
 * the idea is to set river conductance to let rivers take up the cells "drainage demand"
 * (recharge and lateral flow at equilibrium groundwater head)
 * @param current_recharge Groundwater recharge
 * @param eq_head Equilibrium groundwater head
 * @param current_head Current groundwater head
 * @param eq_flow Equilibrium lateral groundwater flow
 * @return
 */
t_s_meter_t ExternalFlow::calcERC(t_vol_t current_recharge,
                                  t_meter eq_head,
                                  t_meter current_head,
                                  t_vol_t eq_flow) const noexcept {
    //possibility to lock conductance equation with former recharge e.g. from steady-state model
    if (lock_recharge) {
        //current_recharge = locked_recharge;
        return locked_conductance * mult;
    }

    //LOG(debug) << "recharge:" << current_recharge.value() << "head:" << eq_head.value() << "StreamStage: " << flowHead.value() << "Bottom" << bottom.value() << "EQFlow" << eq_flow.value() << "AltConduct" << conductance.value();
    t_s_meter_t out = 0 * si::square_meter / day;

    // River conductance of gaining rivers in steady state following Miguez-Macho et al. (2007)
    // sets river conductance to let rivers take up the cells "drainage demand" (recharge and lateral flow at equilibrium groundwater head)
    t_meter stage = eq_head - flowHead;
    NANChecker(stage.value(), "ERC stage problem");
    if (stage.value() <= 0) { stage = .1 * si::meter; }
    // set scale parameter p (not in use)
    t_dim p = 1 * si::si_dimensionless;
    // calculate conductance
    out = (current_recharge * p + eq_flow) / stage;
    NANChecker(out.value(), "ERC Recharge Problem");

    if (out < conductance) {out = conductance;} //Only happens if cell was loosing in eq and is now gaining; stage

    if (current_head < flowHead - 1 * si::meter) { // for losing rivers: use conductance from input data
        out = conductance;
        if (out.value() > 1e+10) { out = 1e+10 * si::square_meter / day; }
        NANChecker(out.value(), "ERC Problem low flow head");
        if (out.value() <= 0) { LOG(critical) << "conductance <= 0"; }
        return out * mult;
    } else if (current_head > flowHead + 1 * si::meter) { // for gaining rivers: use approach by Miguez-Macho et al. (2007)
        if (out.value() > 1e+10) { out = 1e+10 * si::square_meter / day; }
        NANChecker(out.value(), "ERC Problem high flow head");
        if (out.value() <= 0) { LOG(critical) << "conductance <= 0"; }
        return out * mult; // mult is only used for sensitivity analysis
    } else { // when current head and flow head are less than 1 meter apart
        double delta = smoothstep(flowHead.value() - 1, flowHead.value() + 1, current_head.value());
        double range = std::abs(out.value() - conductance.value());
        //double lower_bound = out.value() > conductance.value() ? out.value() : conductance.value();
        out = (out.value() + range * delta) * si::square_meter / day;

        NANChecker(out.value(), "ERC Problem");

        if (out.value() > 1e+10) { out = 1e+10 * si::square_meter / day; }
        if (out.value() <= 0) { LOG(critical) << "conductance <= 0"; }
        //LOG(debug) << "calcERC out = " << out.value() << ", recharge = " << current_recharge.value() <<
        //              ", equilibrium flow = " << eq_flow.value() << ", stage = " << stage.value();
        return out * mult;
    }
}
}
}//ns
