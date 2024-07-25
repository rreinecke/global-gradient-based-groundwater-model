#include "Node.hpp"

BOOST_CLASS_IMPLEMENTATION(GlobalFlow::Model::ExternalFlow, boost::serialization::object_serializable);
BOOST_CLASS_IMPLEMENTATION(GlobalFlow::Model::NodeInterface, boost::serialization::object_serializable);

namespace GlobalFlow {
namespace Model {

/**
 * Initialize the physical properties with default values
 * @return
 */
PhysicalProperties initProperties() {
    PhysicalProperties fields;
    fields.set < large_num, ID > (0);
    fields.set < large_num, SpatID > (0);
    fields.set<double, Lat>(0);
    fields.set<double, Lon>(0);
    fields.set<int, Layer>(0);
    fields.set<quantity<Dimensionless>, StepSize>(1 * si::si_dimensionless);
    fields.set<quantity<SquareMeter>, Area>(0 * si::square_meter);
    fields.set<quantity<Meter>, Elevation>(5 * si::meter);
    fields.set<quantity<Meter>, TopElevation>(5 * si::meter);
    fields.set<quantity<Meter>, VerticalSize>(10 * si::meter);
    fields.set<quantity<Meter>, EFolding>(1 * si::meter);
    fields.set<bool, Confinement>(true);
    fields.set<quantity<Velocity>, K>(0.03 * (si::meter / day));
    fields.set<quantity<Dimensionless>, Anisotropy>(10 * si::si_dimensionless);
    fields.set<quantity<VolumePerTime>, OUT>(0.0 * si::cubic_meter/day);
    fields.set<quantity<VolumePerTime>, IN>(0.0 * si::cubic_meter/day);
    fields.set<quantity<Meter>, Head>(1 * si::meter);
    fields.set<quantity<Meter>, EQHead>(1 * si::meter);
    fields.set<quantity<Meter>, HeadChange>(0 * si::meter);
    fields.set<quantity<Meter>, Head_TZero>(0 * si::meter);
    fields.set<quantity<Meter>, HeadChange_TZero>(0 * si::meter);
    return fields;
}

NodeInterface::NodeInterface(NodeVector nodes,
                             double lat,
                             double lon,
                             quantity<SquareMeter> area,
                             quantity<Meter> edgeLengthLeftRight,
                             quantity<Meter> edgeLengthFrontBack,
                             unsigned long int spatID,
                             unsigned long int identifier,
                             quantity<Velocity> conduct,
                             quantity<Meter> head,
                             double aquiferDepth,
                             double anisotropy,
                             double specificYield,
                             double specificStorage,
                             bool useEfolding,
                             bool confined,
                             bool isSteadyState,
                             bool isDensityVariable,
                             std::vector<quantity<Dimensionless>> delnus,
                             std::vector<quantity<Dimensionless>> nusInZones,
                             double effPorosity,
                             double maxTipSlope,
                             double maxToeSlope,
                             double minDepthFactor,
                             double slopeAdjFactor,
                             quantity<Meter> vdfLock,
                             int sourceZoneGHB,
                             int sourceZoneRecharge): nodes(nodes) {
    fields = initProperties();
    fields.set<double, Lat>(lat);
    fields.set<double, Lon>(lon);
    fields.set<quantity<SquareMeter>, Area>(area);
    fields.set < unsigned long int, SpatID > (spatID);
    fields.set < unsigned long int, ID > (identifier);
    fields.set<quantity<Velocity>, K>(conduct);
    fields.set<quantity<Meter>, Head>(head);
    fields.set<bool, UseEfolding>(useEfolding);
    fields.set<bool, Confinement>(confined);
    fields.set<quantity<Dimensionless>, SpecificYield>(specificYield * si::si_dimensionless);
    fields.set<quantity<perUnit>, SpecificStorage>(specificStorage * perMeter);
    fields.set<quantity<Meter>, VerticalSize>(aquiferDepth * si::meter);
    fields.set<quantity<Dimensionless>, Anisotropy>(anisotropy * si::si_dimensionless);
    fields.set<quantity<Meter>, EdgeLengthLeftRight>(edgeLengthLeftRight);
    fields.set<quantity<Meter>, EdgeLengthFrontBack>(edgeLengthFrontBack);
    fields.set<quantity<SquareMeter>, SurfaceLeftRight>(
            fields.get<quantity<Meter>, EdgeLengthLeftRight>() * fields.get<quantity<Meter>, VerticalSize>());
    fields.set<quantity<SquareMeter>, SurfaceFrontBack>(
            fields.get<quantity<Meter>, EdgeLengthFrontBack>() * fields.get<quantity<Meter>, VerticalSize>());
    fields.set<quantity<CubicMeter>, VolumeOfCell>(
            fields.get<quantity<SquareMeter>, Area>() * fields.get<quantity<Meter>, VerticalSize>());
    fields.set<bool, IsSteadyState> (isSteadyState);
    fields.set<bool, IsDensityVariable> (isDensityVariable);
    fields.set<std::vector<quantity<Dimensionless>>, Delnus> (delnus);
    fields.set<std::vector<quantity<Dimensionless>>, NusInZones> (nusInZones);
    fields.set<quantity<Dimensionless>, EffectivePorosity> (effPorosity * si::si_dimensionless);
    fields.set<quantity<Dimensionless>, MaxTipSlope> (maxTipSlope * si::si_dimensionless);
    fields.set<quantity<Dimensionless>, MaxToeSlope> (maxToeSlope * si::si_dimensionless);
    fields.set<quantity<Dimensionless>, MinDepthFactor> (minDepthFactor * si::si_dimensionless);
    fields.set<quantity<Dimensionless>, SlopeAdjFactor> (slopeAdjFactor * si::si_dimensionless);
    fields.set<quantity<Meter>, VDFLock> (vdfLock);
    fields.set<int, SourceZoneGHB> (sourceZoneGHB);
    fields.set<int, SourceZoneRecharge> (sourceZoneRecharge);
}
}
}//ns
