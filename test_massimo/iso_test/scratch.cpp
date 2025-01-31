TEST(mesh_test, simple_surface){

    int order = 2;  // Quadratic B-spline

    // Define knot vectors with padding
    std::array<std::vector<double>,2> knots;
    knots[0] = {0.0, 1.0};
    knots[1] = {0.0, 1.0};

    // Define weights
    MdArray<double,full_dynamic_extent_t<2>> weights(3,3);
    weights(0,0) = 1.0;  weights(0,1) = 2.0;  weights(0,2) = 1.0;
    weights(1,0) = 2.0;  weights(1,1) = 4.0;  weights(1,2) = 2.0;
    weights(2,0) = 1.0;  weights(2,1) = 2.0;  weights(2,2) = 4.0;

    // Define control points (3x3 grid)
    MdArray<double,full_dynamic_extent_t<3>> control_points(3,3,3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            control_points(i, j, 0) = i - 1; // X-coordinates (-1, 0, 1)
            control_points(i, j, 1) = j - 1; // Y-coordinates (-1, 0, 1)
            control_points(i, j, 2) = (i == 1 && j == 1) ? 2.0 : 0.0; // Peak at (0,0)
        }
    }

    // Create the NURBS surface
    IsoMesh<2,3> mesh(knots, weights, control_points, order);

    // Evaluate at midpoint (u = 0.5, v = 0.5)
    std::array<double, 2> uv = {0.5, 0.5};
    auto result1 = mesh.eval_param(uv);
    auto result2 = mesh.eval_param2(uv);

    auto param_deriv = mesh.eval_param_derivative(uv);
    auto param_deriv2 = mesh.eval_param_derivative2(uv);

    for(int i=0;i<3;i++){
        for(int j=0;j<2;j++){
            EXPECT_TRUE(almost_equal(param_deriv(i,j),param_deriv2(i,j)));
        }
    }

}