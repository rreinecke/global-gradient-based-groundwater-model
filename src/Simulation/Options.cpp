#include "Options.hpp"
#include "../Logging/Logging.hpp"

namespace GlobalFlow {
    namespace Simulation {

        namespace pt = boost::property_tree;

        inline std::vector<std::string> getArray(std::string child_name, pt::ptree subtree) {
            std::vector<std::string> out;
            pt::ptree array = subtree.get_child(child_name);
            BOOST_FOREACH(const pt::ptree::value_type &child, array) {
                assert(child.first.empty());
                std::string t = child.second.get<std::string>("");
                out.push_back(child_name + "/" + t);
            }
            return out;
        }

        template<class T>
        inline std::vector<T> getTypeArray(std::string child_name, pt::ptree subtree) {
            std::vector<T> out;
            pt::ptree array = subtree.get_child(child_name);
            BOOST_FOREACH(const pt::ptree::value_type &child, array) {
                assert(child.first.empty());
                T t = child.second.get<T>("");
                out.push_back(t);
            }
            return out;
        }


        inline std::string getOptional(std::string child_name, pt::ptree subtree) {
            boost::optional<std::string> data = subtree.get_optional<std::string>(child_name);
            if (data) {
                return data.value();
            }
            return std::string("");
        }

        void
        Options::load(const std::string &filename) {
            pt::ptree tree;
            pt::read_json(filename, tree);
            tree = tree.get_child("config");

            pt::ptree config = tree.get_child("model_config");
            STRESS_PERIOD_STEADY_STATE = getTypeArray<bool>("stress_period_steady_state", config);
            STRESS_PERIOD_STEPS = getTypeArray<int>("stress_period_time_steps", config);
            STRESS_PERIOD_STEP_SIZES = getTypeArray<std::string>("stress_period_time_step_sizes", config);
            STRESS_PERIOD_VARIABLE_DENSITY = getTypeArray<bool>("stress_period_variable_density", config);
            NODES = config.get<std::string>("nodes");
            NUMBER_OF_NODES_PER_LAYER = config.get<unsigned long int>("number_of_nodes_per_layer");
            Y_RANGE = config.get<long>("y_range");
            X_RANGE = config.get<long>("x_range");
            IS_GLOBAL = config.get<bool>("is_global");
            RESOLUTION_IN_DEGREE = config.get<double>("resolution_in_degree");
            EDGE_LENGTH_LEFT_RIGHT = config.get<double>("edge_length_left_right");
            EDGE_LENGTH_FRONT_BACK = config.get<double>("edge_length_front_back");
            LAYERS = config.get<int>("layers");
            USE_EFOLDING = config.get<bool>("use_efolding");
            CONFINED = getTypeArray<bool>("confinement", config);
            if (LAYERS != CONFINED.size()) {
                LOG(critical) << "mismatching layers";
                exit(3);
            }

            DEFAULT_BOUNDARY_CONDITION = config.get<std::string>("default_boundary_condition");
            SENSITIVITY = config.get<bool>("sensitivity");

            pt::ptree numerics = tree.get_child("numerics");
            THREADS = numerics.get<int>("threads");
            SOLVER = numerics.get<std::string>("solver");
            MAX_INNER_ITERATIONS = numerics.get<int>("max_inner_iterations");
            MAX_OUTER_ITERATIONS_HEAD = numerics.get<int>("max_outer_iterations_head");
            MAX_OUTER_ITERATIONS_ZETA  = numerics.get<int>("max_outer_iterations_zeta");
            RCLOSE_HEAD = numerics.get<double>("closing_crit_head");
            RCLOSE_ZETA = numerics.get<double>("closing_crit_zeta");
            MAX_HEAD_CHANGE = numerics.get<double>("head_change");
            MAX_ZETA_CHANGE = numerics.get<double>("zeta_change");
            DAMPING = numerics.get<bool>("damping");
            MIN_DAMP = numerics.get<double>("min_damp");
            MAX_DAMP = numerics.get<double>("max_damp");

            pt::ptree input = tree.get_child("input");

            pt::ptree data_config = input.get_child("data_config");
            k_from_file = data_config.get<bool>("k_from_file");
            k_ghb_from_file = data_config.get<bool>("k_ghb_from_file");
            specificstorage_from_file = data_config.get<bool>("specific_storage_from_file");
            specificyield_from_file = data_config.get<bool>("specific_yield_from_file");
            k_river_from_file = data_config.get<bool>("k_river_from_file");
            aquifer_depth_from_file = data_config.get<bool>("aquifer_depth_from_file");
            eq_wtd_from_file = data_config.get<bool>("eq_wtd_from_file");
            initial_head_from_file = data_config.get<bool>("initial_head_from_file");
            effective_porosity_from_file = data_config.get<bool>("effective_porosity_from_file");
            zones_sources_from_file = data_config.get<bool>("zones_sources_from_file");

            pt::ptree default_data = input.get_child("default_data");
            K = getTypeArray<double>("K", default_data);
            INITIAL_HEAD = default_data.get<double>("initial_head");
            GHB_K = default_data.get<double>("ghb_K");
            RIVER_CONDUCTIVITY = default_data.get<double>("river_conductivity");
            SWB_ELEVATION_FACTOR = default_data.get<double>("swb_elevation_factor");
            AQUIFER_DEPTH = getTypeArray<int>("aquifer_thickness", default_data);

            ANISOTROPY = getTypeArray<double>("anisotropy", default_data);
            SPECIFIC_YIELD = default_data.get<double>("specific_yield");
            SPECIFIC_STORAGE = default_data.get<double>("specific_storage");

            EFFECTIVE_POROSITY = default_data.get<double>("effective_porosity");
            SOURCE_ZONE_GHB = default_data.get<int>("source_zone_ghb");
            SOURCE_ZONE_RECHARGE = default_data.get<int>("source_zone_recharge");

            pt::ptree vdf = tree.get_child("vdf_config");
            DENSITY_ZONES = getTypeArray<double>("density_zones", vdf);
            MAX_TIP_SLOPE = vdf.get<double>("max_tip_slope");
            MAX_TOE_SLOPE = vdf.get<double>("max_toe_slope");
            MIN_DEPTH_FACTOR = vdf.get<double>("min_depth_factor");
            SLOPE_ADJ_FACTOR = vdf.get<double>("slope_adj_factor");
            VDF_LOCK = vdf.get<double>("vdf_lock");

            pt::ptree data = input.get_child("data");

            bool efoldAsArray = data_config.get<bool>("efold_as_array");
            if (efoldAsArray){
                EFOLDING_a = getArray("E-Folding", data.get_child("e-folding"));
            }
            EFOLDING = getOptional("e-folding", data);

            INITIAL_ZETAS_AS_ARRAY = data_config.get<bool>("initial_zetas_as_array");
            if (INITIAL_ZETAS_AS_ARRAY){
                INITIAL_ZETAS_a = getArray("Zetas", data.get_child("zetas"));
            }
            INITIAL_ZETAS = getOptional("zetas", data);

            ELEVATION = getOptional("elevation", data);
            EQUAL_WATER_TABLE_DEPTH = getOptional("equal_water_table_depth", data);
            RIVER_ELEVATION = getOptional("river_elevation", data);

            LITHOLOGY = getOptional("lithology", data);
            RECHARGE = getOptional("recharge", data);
            ZONES_SOURCES_FILE = getOptional("zones_sources", data);
            PSEUDO_SOURCE_FLOW = getOptional("pseudo_source_flow", data);
            RIVER = getOptional("river_extent", data);
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

            INITIAL_ZONES = getOptional("initial_zones", data);

            EFFECTIVE_POROSITY_FILE = getOptional("effective_porosity", data);

            SPATID_ARCID = getOptional("spatID-arcID", data);
        }

    }
}//ns
