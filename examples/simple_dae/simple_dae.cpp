/*
 * Solves a very simple system of differential algebraic equation as a test:
 *
 * x' = y
 * 0  = x*x + y*y - 1
 *
 * Initial conditions are: x = 0, y = 1 for t = 0.
 *
 * The solution of this system is
 *
 * x = sin(t), y = cos(t), 0 <= t <= pi/2;
 * x = 1, y = 0, t > pi/2.
 *
 * Each time step we will check that
 * (1) x*x + y*y = 1 for any t, and
 * (2) x(t) = sin(t) for t <= pi/2, x(t) = 1 for t > pi/2
 *
 * with the absolute tolerance at least 1e-6.
 */

#include <iostream>
#include <cmath>
#include <algorithm>  // for std::max_element

#include "../../src/solver.h"  // the main header of dae-cpp library solver

using namespace daecpp;

// python3 + numpy + matplotlib should be installed in order to enable plotting
// #define PLOTTING

#ifdef PLOTTING
#include "../../src/external/matplotlib-cpp/matplotlibcpp.h"
namespace plt = matplotlibcpp;
#endif

/*
 * Singular mass matrix in simplified 3-array sparse format
 * =============================================================================
 * The matrix has the following form:
 * M = |1 0|
 *     |0 0|
 */
class MyMassMatrix : public MassMatrix
{
public:
    void operator()(sparse_matrix_holder &M)
    {
        M.A.resize(2);   // Number of non-zero and diagonal elements
        M.ia.resize(2);  // Number of non-zero and diagonal elements
        M.ja.resize(2);  // Number of non-zero and diagonal elements

        // Non-zero and diagonal elements
        M.A[0] = 1;
        M.A[1] = 0;

        // Column index of each element given above
        M.ja[0] = 0;
        M.ja[1] = 1;

        // Row index of each element in M.A:
        M.ia[0] = 0;
        M.ia[1] = 1;
    }
};

/*
 * RHS of the problem
 * =============================================================================
 */
class MyRHS : public RHS
{
public:
    /*
     * Receives current solution vector x and the current time t.
     * Defines the RHS f.
     */
    void operator()(const state_type &x, state_type &f, const double t)
    {
        f[0] = x[1];
        f[1] = x[0] * x[0] + x[1] * x[1] - 1.0;
    }
};

/*
 * (Optional) Observer
 * =============================================================================
 * Every time step checks that
 * (1) x*x + y*y = 1, and
 * (2) x(t) - sin(t) = 0 for t <= pi/2, x(t) = 1 for t > pi/2
 * and prints solution and errors to console.
 */
class MySolver : public Solver
{
public:
    MySolver(RHS &rhs, Jacobian &jac, MassMatrix &mass, SolverOptions &opt)
        : Solver(rhs, jac, mass, opt)
    {
    }

#ifdef PLOTTING
    state_type x_axis, x0, x1;  // For plotting
#endif
    state_type err1, err2;  // To check errors

    /*
     * Overloaded observer.
     * Receives current solution vector and the current time every time step.
     */
    void observer(state_type &x, const double t)
    {
        double e1 = std::abs(x[0] * x[0] + x[1] * x[1] - 1.0);
        double e2 = 0;

        if(t <= 1.5707963)
            e2 = std::abs(std::sin(t) - x[0]);
        else
            e2 = std::abs(x[0] - 1.0);

        std::cout << t << '\t' << x[0] << '\t' << x[1] << '\t' << e1 << '\t'
                  << e2 << '\n';

        err1.push_back(e1);
        err2.push_back(e2);

#ifdef PLOTTING
        // Save data for plotting
        x_axis.push_back(t);
        x0.push_back(x[0]);
        x1.push_back(x[1]);
#endif
    }
};

/*
 * (Optional) Analytical Jacobian in simplified 3-array sparse format
 * =============================================================================
 */
class MyJacobian : public Jacobian
{
public:
    explicit MyJacobian(RHS &rhs) : Jacobian(rhs) {}

    /*
     * Receives current solution vector x and the current time t. Defines the
     * analytical Jacobian matrix J.
     */
    void operator()(sparse_matrix_holder &J, const state_type &x,
                    const double t)
    {
        // Initialize Jacobian in simplified sparse format
        J.A.resize(4);
        J.ia.resize(4);
        J.ja.resize(4);

        // Non-zero and diagonal elements
        J.A[0] = 0.0;
        J.A[1] = 1.0;
        J.A[2] = 2.0 * x[0];
        J.A[3] = 2.0 * x[1];

        // Column index of each element given above
        J.ja[0] = 0;
        J.ja[1] = 1;
        J.ja[2] = 0;
        J.ja[3] = 1;

        // Row index of each non-zero or diagonal element of A
        J.ia[0] = 0;
        J.ia[1] = 0;
        J.ia[2] = 1;
        J.ia[3] = 1;
    }
};

/*
 * MAIN FUNCTION
 * =============================================================================
 * Returns '0' if solution comparison is OK or '1' if solution error is above
 * the acceptable tolerances.
 */
int main()
{
    // Solution time 0 <= t <= pi
    double t1 = 3.14;

    // Define the state vector
    state_type x(2);

    // Initial conditions
    x[0] = 0;
    x[1] = 1;

    // Set up the RHS of the problem.
    // Class MyRHS inherits abstract RHS class from dae-cpp library.
    MyRHS rhs;

    // Set up the Mass Matrix of the problem.
    // MyMassMatrix inherits abstract MassMatrix class from dae-cpp library.
    MyMassMatrix mass;

    // Create an instance of the solver options and update some of the solver
    // parameters defined in solver_options.h
    SolverOptions opt;

    opt.dt_init = 1.0e-2;   // Change the initial time step.
                            // It should be relatively small, because the first
                            // step in time is first order accuracy.
                            // Reducing dt_init decreases the error (2)
    opt.time_stepping = 1;  // Use simple stability-based adaptive time stepping
                            // algorithm.
    opt.bdf_order = 6;      // Use BDF-6
    opt.verbosity = 0;      // Suppress output to screen (we have our own output
                            // defined in Observer function above)

    // We can override Jacobian class from dae-cpp library and provide
    // analytical Jacobian.
    MyJacobian jac(rhs);

    // Or we can use numerically estimated Jacobian with the given tolerance.
    // Jacobian jac_est(rhs, 1e-6);

    // This will decrease the error (1) for x*x + y*y = 1
#ifdef DAE_SINGLE
    opt.atol = 1e-6;  // Redefine absolute tolerance for single precision
    opt.rtol = 1e-6;  // Redefine relative tolerance for single precision
#else
    opt.atol = 1e-8;  // Redefine absolute tolerance for double precision
    opt.rtol = 1e-8;  // Redefine relative tolerance for double precision
#endif

    // Create an instance of the solver with particular RHS, Mass matrix,
    // Jacobian and solver options
    MySolver solve(rhs, jac, mass, opt);

    // Now we are ready to solve the set of DAEs
    std::cout << "\nStarting DAE solver...\n";
    std::cout << "time\tx\ty\terror1\terror2\n";

    // Solve the system
    int status = solve(x, t1);

    // Check errors
    double max_err1 = *std::max_element(solve.err1.begin(), solve.err1.end());
    double max_err2 = *std::max_element(solve.err2.begin(), solve.err2.end());

    std::cout << "\nMaximum absoulte error (1) x*x + y*y = 1: " << max_err1
              << '\n';
    std::cout << "Maximum absolute error (2) x(t) - sin(t) = 0 for t <= pi/2 "
                 "or x(t) = 1 for t > pi/2: "
              << max_err2 << '\n';

    // Plot the solution
#ifdef PLOTTING
    plt::figure();
    plt::figure_size(640, 480);
    plt::named_semilogx("x", solve.x_axis, solve.x0);
    plt::named_semilogx("y", solve.x_axis, solve.x1);
    plt::xlabel("time");
    plt::title("Simple 2x2 DAE system");
    plt::grid(true);
    plt::legend();

    // Save figure
    const char *filename = "simple_dae.png";
    std::cout << "Saving result to " << filename << "...\n";
    plt::save(filename);
#endif

#ifdef DAE_SINGLE
    const bool check_result = (max_err1 > 1e-6 || max_err2 > 1e-6 || status);
#else
    const bool check_result = (max_err1 > 1e-15 || max_err2 > 1e-6 || status);
#endif

    if(check_result)
        std::cout << "...Test FAILED\n\n";
    else
        std::cout << "...done\n\n";

    return check_result;
}
