#include <gtest/gtest.h>
#include <random>
#include <algorithm>
#include "../../../src/Simulation/Stepper.hpp"
#include "../../../src/Simulation/Options.hpp"
#include "../Solver/Equation__MOCK.hpp"

using namespace GlobalFlow::Simulation;

class StepperFixture : public ::testing::Test {
public:
    NodeVector nodes;
    Options op;
    void SetUp(){
        NodeVector ptr(new std::vector<std::unique_ptr<GlobalFlow::Model::NodeInterface>>);
        nodes = std::move(ptr);
        nodes->emplace_back(new GlobalFlow::Model::StandardNode(
                nodes, 0, 0, 1 * si::square_meter, 1 * si::meter, 1 * si::meter, 0, 0, 0.1 * si::meter / day,
                1 * si::meter, 10, 1,
                0.2, 0.1, true, true, true, false, {0, 0.25}, {0, 0.025}, 0.2, 0.1, 0.1,
                0.001, 0.4, 0.001 * si::meter, 0, 0
                ));
        nodes->emplace_back(new GlobalFlow::Model::StandardNode(
                nodes, 1, 0, 1 * si::square_meter, 1 * si::meter, 1 * si::meter, 1, 1, 0.2 * si::meter / day,
                1 * si::meter, 10, 1,
                0.2, 0.1, true, true, true, false, {0, 0.25}, {0, 0.025}, 0.2, 0.1, 0.1,
                0.001, 0.4, 0.001 * si::meter, 0, 0
                ));
        nodes->emplace_back(new GlobalFlow::Model::StandardNode(
                nodes, 0, 1, 1 * si::square_meter, 1 * si::meter, 1 * si::meter, 2, 2, 0.1 * si::meter / day,
                1 * si::meter, 10, 1,
                0.2, 0.1, true, true, true, false, {0, 0.25}, {0, 0.025}, 0.2, 0.1, 0.1,
                0.001, 0.4, 0.001 * si::meter, 0, 0
                ));
        nodes->emplace_back(new GlobalFlow::Model::StandardNode(
                nodes, 1, 1, 1 * si::square_meter, 1 * si::meter, 1 * si::meter, 3, 3, 0.1 * si::meter / day,
                1 * si::meter, 10, 1,
                0.2, 0.1, true, true, true, false, {0, 0.25}, {0, 0.025}, 0.2, 0.1, 0.1,
                0.001, 0.4, 0.001 * si::meter, 0, 0
                ));

        nodes->at(0)->setNeighbour(1,RIGHT);
        nodes->at(1)->setNeighbour(0,LEFT);

        nodes->at(0)->setNeighbour(2,DOWN);
        nodes->at(2)->setNeighbour(0, TOP);

        nodes->at(2)->setNeighbour(3,DOWN);
        nodes->at(3)->setNeighbour(2, TOP);
    }

    std::vector<int> make_rnd(){
        std::random_device rnd_device;
        std::mt19937 mersenne_engine {rnd_device()};
        std::uniform_int_distribution<int> dist {2, 30};

        auto gen = [&dist, &mersenne_engine](){
            return dist(mersenne_engine);
        };
        std::vector<int> vec(10);
        std::generate(begin(vec), end(vec), gen);
        return vec;
    }

    volatile int writeMe{0};
    void IwillCrash(){
        MockEquation equation(nodes, op);

        Stepper stepper = Stepper(reinterpret_cast<GlobalFlow::Solver::Equation *>(&equation), "MONTH", true, false, 1);
        for (step step : stepper) {writeMe = static_cast<int>(10);}
    }

};

//FIXME currently boost log causes a free(): invalid pointer
TEST_F(StepperFixture,DayLoop){
    MockEquation equation(nodes, op);
    Stepper stepper = Stepper(reinterpret_cast<GlobalFlow::Solver::Equation *>(&equation), "DAY", true, false, 2);
    int p{0};
    double a{0};
    for (step step : stepper) {
        a = step.second;
        ++p;
    }
    ASSERT_EQ(p,2);
    ASSERT_EQ(a,1);
    //FIXME currently not possible to test as stepper holds an equation pointer; intro of abstract EQ would solve this
}

TEST_F(StepperFixture,MonthLoop){
    MockEquation equation(nodes, op);
    Stepper stepper = Stepper(reinterpret_cast<GlobalFlow::Solver::Equation *>(&equation), "MONTH", true, false, 10);
    ASSERT_EQ(stepper.getStepSize(),30.0);
    int p{0};
    double a{0};
    for (step step : stepper) {
        a = step.second;
        ++p;
    }
    ASSERT_EQ(p,10);
    ASSERT_EQ(a,9);
}

using DeathStepperFixture = StepperFixture;

TEST_F(DeathStepperFixture,DynmicSteps){
    testing::FLAGS_gtest_death_test_style="threadsafe";
    ASSERT_DEATH(IwillCrash(), "\"Dynamic steps not valid for 1 step\"' failed");
}

TEST_F(StepperFixture,DynmicStepsRand){
    MockEquation equation(nodes, op);

    for(int p : make_rnd()){
        Stepper stepper = Stepper(reinterpret_cast<GlobalFlow::Solver::Equation *>(&equation), "MONTH", true, false, p);
        ASSERT_EQ(stepper.getStepSize(),30.0);
        int i{0};
        double a{0};
        for (step step : stepper) {
            a = step.second;
            ++i;
        }
        ASSERT_EQ(i,p); //it will take p-1 iterations to run through complete time frame
        ASSERT_DOUBLE_EQ(p,a);
    }
}