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
            ADAPTIVE_STEP_SIZE = config.get<bool>("adaptive_step_size");
            BOUNDARY_CONDITION = config.get<string>("boundary_condition");
            SENSITIVITY = config.get<bool>("sensitivity");

            pt::ptree numerics = tree.get_child("numerics");
            SOLVER = numerics.get<string>("solver");
            IITER = numerics.get<int>("iterations");
            I_ITTER = numerics.get<int>("inner_itter");
            RCLOSE = numerics.get<double>("closing_crit");
            MAX_HEAD_CHANGE = numerics.get<double>("head_change");
            MAX_ZETA_CHANGE = numerics.get<double>("zeta_change");
            DAMPING = numerics.get<bool>("damping");
            MIN_DAMP = numerics.get<double>("min_damp");
            MAX_DAMP = numerics.get<double>("max_damp");
            string tmp = numerics.get<string>("step_size");
            if (tmp == "DAILY") {
                step_size = DAILY;
            }
            if (tmp == "MONTHLY") {
                step_size = MONTHLY;
            }
            WETTING_APPROACH = numerics.get<string>("wetting_approach");

            pt::ptree input = tree.get_child("input");

//BASE_PATH = input.get<string>("base_path");

            pt::ptree data_config = input.get_child("data_config");
            k_from_lith = data_config.get<bool>("k_from_lith");
            k_ghb_from_file = data_config.get<bool>("k_ghb_from_file");
            specificstorage_from_file = data_config.get<bool>("specific_storage_from_file");
            specificyield_from_file = data_config.get<bool>("specific_yield_from_file");
            k_river_from_file = data_config.get<bool>("k_river_from_file");
            aquifer_depth_from_file = data_config.get<bool>("aquifer_depth_from_file");
            eq_wtd_from_file = data_config.get<bool>("eq_wtd_from_file");
            initial_head_from_file = data_config.get<bool>("initial_head_from_file");
            initial_zetas_from_file = data_config.get<bool>("initial_zetas_from_file");
            effective_porosity_from_file = data_config.get<bool>("effective_porosity_from_file");
            zones_sources_sinks_from_file = data_config.get<bool>("zones_sources_sinks_from_file");

            pt::ptree default_data = input.get_child("default_data");
            K = default_data.get<double>("K");
            INITIAL_HEAD = default_data.get<double>("initial_head");
            GHB_K = default_data.get<double>("ghb_K");
            AQUIFER_DEPTH = getTypeArray<int>("aquifer_thickness", default_data);
            INITIAL_ZETAS = getTypeArray<double>("initial_zetas", default_data);

            ANISOTROPY = default_data.get<double>("anisotropy");
            SPECIFIC_YIELD = default_data.get<double>("specific_yield");
            SPECIFIC_STORAGE = default_data.get<double>("specific_storage");

            DENSITY_VARIABLE = config.get<bool>("density_variable");
            DENSITY_ZONES = getTypeArray<double>("density_zones", config);
            MAX_TIP_TOE_SLOPE = config.get<double>("max_tip_toe_slope");
            MIN_DEPTH_FACTOR = config.get<double>("min_depth_factor");
            SLOPE_ADJ_FACTOR = config.get<double>("slope_adj_factor");
            VDF_LOCK = config.get<double>("vdf_lock");

            EFFECTIVE_POROSITY = default_data.get<double>("effective_porosity");
            ZONES_SOURCES_SINKS = getTypeArray<int>("zones_sources_sinks", default_data);


            bool slopeAsArray = data_config.get<bool>("slope_as_array");
            bool efoldAsArray = data_config.get<bool>("efold_as_array");
            pt::ptree data = input.get_child("data");

            if (slopeAsArray) {
                SLOPE_a = getArray("Slope", data.get_child("slope"));
                //ELEVATION_a = getArray("Elevation", data.get_child("elevation"));
                //EQUAL_WATER_TABLE_DEPTH_a = getArray("WTD", data.get_child("equal_water_table_depth"));
                //SURFACE_WATER_ELEVATION_a = getArray("SurfaceWaterElevation", data.get_child("surface_water_elevation"));
            }
            SLOPE = getOptional("slope", data);


            if (efoldAsArray){
                EFOLDING_a = getArray("E-Folding", data.get_child("e-folding"));
            }
            EFOLDING = getOptional("e-folding", data);

            ELEVATION = getOptional("elevation", data);
            EQUAL_WATER_TABLE_DEPTH = getOptional("equal_water_table_depth", data);
            SURFACE_WATER_ELEVATION = getOptional("surface_water_elevation", data);


            LITHOLOGY = getOptional("lithology", data);
            RECHARGE = getOptional("recharge", data);
            ZONES_SOURCES_SINKS_FILE = getOptional("zones_sources_sinks", data);
            PSEUDO_SOURCE_FLOW = getOptional("pseudo_source_flow", data);
            RIVER = getOptional("river", data);
            GLOBAL_WETLANDS = getOptional("global_wetlands", data);
            GLOBAL_LAKES = getOptional("global_lakes", data);
            LOCAL_LAKES = getOptional("local_lakes", data);
            LOCAL_WETLANDS = getOptional("local_wetlands", data);

            //Optional
            K_DIR = getOptional("conductance", data);
            RIVER_K = getOptional("river_conductance", data);
            GHB_K_DIR = getOptional("ghb_conductance", data);
            SS_FILE = getOptional("specific_storage", data);
            SY_FILE = getOptional("specific_yield", data);
            AQ_DEPTH = getOptional("aquifer_depth", data);

            INITIAL_HEAD_FILE = getOptional("initial_head", data);

            INITIAL_ZETAS_FILE = getOptional("initial_zetas", data);

            INITIAL_ZONES = getOptional("initial_zones", data);

            EFFECTIVE_POROSITY_FILE = getOptional("effective_porosity", data);

            SPATID_ARCID = getOptional("spatID-arcID", data);
        }

    }
}//ns
