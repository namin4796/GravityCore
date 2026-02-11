#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "../src/gravity.hpp"
#include <vector>

// test 1: does the universe init correctly
TEST_CASE("Universe Initialization", "[core]") {
    int N = 100;
    Universe sim(N, 100.0, 0.001, 1.0, 0.01);

    // check if memory was allocated
    REQUIRE(sim.get_x().size() == N);
    REQUIRE(sim.get_mass().size() == N);
}

// test 2: does gravity force calculation work for
// simple 2-body case
TEST_CASE("Two Body Attraction", "[physics]") {
    Universe sim(2, 100.0, 0.0, 1.0, 1.0); // No DM, rho=0

    std::vector<double> x = {0.0, 10.0};
    std::vector<double> y = {0.0, 0.0};
    std::vector<double> vx = {0.0, 0.0};
    std::vector<double> vy = {0.0, 0.0};
    std::vector<double> mass = {1000.0, 1000.0};

    sim.set_state(x, y, vx, vy, mass);

    // run 1 step
    sim.step();

    std::vector<double> new_vx = sim.get_vx();

    //particle 0 should not move; central BH
    //therefore new_vx > 0.0 (expectation)
    REQUIRE(new_vx[0] == 0.0);
    //particle 1 should move LEFT (-), attracted to BH
    //therefore new_vx < 0.0 (expectation)
    REQUIRE(new_vx[1] < 0.0);
}

// test 3: checking quadtree integrity by checking if 
// creating a sim doesn't crash the memory pool
TEST_CASE("Memory Pool Stress Test", "[stress]") {
    int N = 5000;
    Universe sim(N, 100.0, 0.001, 1.0, 0.01);

    std::vector<double> x(N, 0.0);
    std::vector<double> y(N, 0.0);
    std::vector<double> vx(N, 0.0);
    std::vector<double> vy(N, 0.0);
    std::vector<double> mass(N, 1.0);

    //random positions
    for (int i=0; i<N; i++) {
        x[i] = (i % 100) * 1.0;
        y[i] = (i % 100) * 1.0;
    }

    sim.set_state(x, y, vx, vy, mass);

    // if the mem. pool is broken this will give a segfault!
    REQUIRE_NOTHROW(sim.step());
}
