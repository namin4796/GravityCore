#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "../src/gravity.hpp" // Adjust path if necessary

TEST_CASE("Galaxy Initialization and Mass Allocation", "[universe]") {
    // Universe(n, rs, rho, g, dt)
    Universe engine(100, 100.0, 0.01, 2.0, 0.01);
    
    // init_galaxy(radius, central_mass)
    engine.init_galaxy(500.0, 50000.0);

    SECTION("Vectors are sized correctly") {
        REQUIRE(engine.get_x().size() == 100);
        REQUIRE(engine.get_y().size() == 100);
    }

    SECTION("Central Black Hole is initialized correctly") {
        auto masses = engine.get_mass();
        auto vx = engine.get_vx();
        auto vy = engine.get_vy();
        
        REQUIRE(masses[0] == 50000.0); // Check central mass
        REQUIRE(vx[0] == 0.0);         // Center should be stationary
        REQUIRE(vy[0] == 0.0);
    }
}

TEST_CASE("Barnes-Hut Quadtree Generation", "[quadtree]") {
    Universe engine(500, 100.0, 0.01, 2.0, 0.01);
    engine.init_galaxy(500.0, 50000.0);
    
    // Taking a step forces the compute_forces() function to run,
    // which builds the Quadtree and populates the debug rectangles.
    engine.step(); 
    
    auto rects = engine.get_tree_rects();
    
    SECTION("Tree successfully partitioned space") {
        REQUIRE(rects.size() > 0);
        // Rectangles are stored as flat triplets: [x, y, size, x, y, size...]
        REQUIRE(rects.size() % 3 == 0); 
    }
}
