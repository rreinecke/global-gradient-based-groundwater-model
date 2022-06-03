#include "Options.hpp"
#include "../Logging/Logging.hpp"

namespace GlobalFlow {
namespace Simulation {

namespace pt = boost::property_tree;
using namespace std;

inline vector<string> getArray(string child_name, pt::ptree subtree) {
    vector<string> out;
    pt::ptree array = subtree.get_child(child_name);
    BOOST_FOREACH(const pt::ptree::value_type &child, array) {
                    assert(child.first.empty());
                    string t = child.second.get<string>("");
                    out.push_back(child_name + "/" + t);
                }
    return out;
}

template<class T>
inline vector<T> getTypeArray(string child_name, pt::ptree subtree) {
    vector<T> out;
    pt::ptree array = subtree.get_child(child_name);
    BOOST_FOREACH(const pt::ptree::value_type &child, array) {
                    assert(child.first.empty());
                    T t = child.second.get<T>("");
                    out.push_back(t);
                }
    return out;
}


inline string getOptional(string child_name, pt::ptree subtree) {
    boost::optional<string> data = subtree.get_optional<string>(child_name);
    if (data) {
        return data.value();
    }
    return string("");
}

void
Options::load(const std::string &filename) {
    pt::ptree tree;
    pt::read_json(filename, tree);
    tree = tree.get_child("config");

    pt::ptree config = tree.get_child("model_config");
    NODES = config.get<string>("nodes");
    ROW_COLS = config.get<bool>("row_cols");
    NUMBER_OF_NODES = config.get<long>("number_of_nodes");
    NUMBER_OF_ROWS = config.get<long>("number_of_rows");
    NUMBER_OF_COLS = config.get<long>("number_of_cols");
    EDGE_LENGTH_ROWS = config.get<double>("edge_length_rows");
    EDGE_LENGTH_COLS = config.get<double>("edge_length_cols");
    THREADS = config.get<int>("threads");
    LAYERS = config.get<int>("layers");
    ONE_LAYER = config.get<bool>("one_layer_approach");
    CONFINED = getTypeArray<bool>("confinement", config);
    if (LAYERS != CONFINED.size()) {
        LOG(critical) << "mismatching layers";
        exit(3);
    }
    //if (ONE_LAYER and LAYERS > 1) {
    //    LOG(critical) << "Approach only viable with one layer";
    //    exit(3);
    //}
    CACHE = config.get<bool>("cache");
    ADAPTIVE_STEPSIZE = config.get<bool>("adaptivestepsize");
    BOUNDARY_CONDITION = config.get<string>("boundarycondition");
    SENSITIVITY = config.get<bool>("sensitivity");

    pt::ptree numerics = tree.get_child("numerics");
    SOLVER = numerics.get<string>("solver");
    IITER = numerics.get<int>("iterations");
    I_ITTER = numerics.get<int>("inner_itter");
    RCLOSE = numerics.get<double>("closingcrit");
    MAX_HEAD_CHANGE = numerics.get<double>("headchange");
    DAMPING = numerics.get<bool>("damping");
    MIN_DAMP = numerics.get<double>("min_damp");
    MAX_DAMP = numerics.get<double>("max_damp");
    string tmp = numerics.get<string>("stepsize");
    if (tmp == "DAILY") {
        stepsize = DAILY;
    }
    if (tmp == "MONTHLY") {
        stepsize = MONTHLY;
    }
    WETTING_APPROACH = numerics.get<string>("wetting_approach");

    pt::ptree input = tree.get_child("input");

//BASE_PATH = input.get<string>("base_path");

    pt::ptree data_config = input.get_child("data_config");
    k_from_lith = data_config.get<bool>("k_from_lith");
    k_ghb_from_file = data_config.get<bool>("k_ghb_from_file");
    specificstorage_from_file = data_config.get<bool>("specificstorage_from_file");
    specificyield_from_file = data_config.get<bool>("specificyield_from_file");
    k_river_from_file = data_config.get<bool>("k_river_from_file");
    aquifer_depth_from_file = data_config.get<bool>("aquifer_depth_from_file");
    //heads_from_file = data_config.get<bool>("initial_head_from_file");

    bool asArray = data_config.get<bool>("data_as_array");

    pt::ptree default_data = input.get_child("default_data");
    K = default_data.get<double>("K");
    INITAL_HEAD = default_data.get<double>("initial_head");
    GHB_K = default_data.get<double>("ghb_K");
    AQUIFER_DEPTH = getTypeArray<int>("aquifer_thickness", default_data);
    ANISOTROPY = default_data.get<double>("anisotropy");
    SPECIFIC_YIELD = default_data.get<double>("specificyield");
    SPECIFIC_STORAGE = default_data.get<double>("specificstorage");

    DENSITY_VARIABLE = config.get<bool>("density_variable");
    DENSITY_STRATIFIED = config.get<bool>("density_stratified");
    NUMBER_OF_DENSITY_ZONES = config.get<int>("number_of_density_zones");
    DENSITIES = getTypeArray<double>("densities", config);

    pt::ptree data = input.get_child("data");

    if (asArray) {
        //ELEVATION_a = getArray("Elevation", data.get_child("elevation"));
        SLOPE_a = getArray("Slope", data.get_child("slope"));
        //EQ_WTD_a = getArray("WTD", data.get_child("eq_wtd"));
        EFOLDING_a = getArray("E-Folding", data.get_child("e-folding"));
        //BLUE_ELEVATION_a = getArray("Blue", data.get_child("blue_cells"));
    } else {
        EFOLDING = getOptional("e-folding", data);
        SLOPE = getOptional("slope", data);
    }

    ELEVATION = getOptional("elevation", data);
    EQ_WTD = getOptional("eq_wtd", data);
    BLUE_ELEVATION = getOptional("blue_cells", data);


    LITHOLOGY = getOptional("lithology", data);
    RECHARGE = getOptional("recharge", data);
    PSEUDO_SOURCE_FLOW = getOptional("pseudo_source_flow", data);
    RIVERS = getOptional("rivers", data);
    GLOBAL_WETLANDS = getOptional("globalwetlands", data);
    GLOBAL_LAKES = getOptional("globallakes", data);
    LOCAL_LAKES = getOptional("lokallakes", data);
    LOCAL_WETLANDS = getOptional("lokalwetlands", data);

    //Optional
    K_DIR = getOptional("conductance", data);
    RIVER_K_DIR = getOptional("river_conductance", data);
    GHB_K_DIR = getOptional("ghb_conductance", data);
    SS_FILE = getOptional("specificstorage", data);
    SY_FILE = getOptional("specificyield", data);
    AQ_DEPTH = getOptional("aquiferdepth", data);

    INITIAL_HEADS = getOptional("initial_head", data);


    boost::optional<pt::ptree &> mappings = input.get_child_optional("mapping");
    if (mappings) {
        NODEID_SPATID = mappings.value().get<string>("nodeID-spatID");
    }
}

}
}//ns
