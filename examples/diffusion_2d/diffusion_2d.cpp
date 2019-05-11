/*
 * This example solves the system of ODEs that describe diffusion in 2D plane:
 *
 * dC/dt = D*(d/dx(dC/dx) + d/dy(dC/dy)),
 *
 * where C is the concentration (dimensionless) one the square unit domain,
 * 0 <= x <= 1 and 0 <= y <= 1. D is the diffusion coefficient.
 *
 * Initial condition is: C(x,y,t=0) = delta_function(x-1/2,y-1/2).
 * Boundary conditions: dC/dx = dC/dy = 0 on the boundaries.
 *
 * The system will be resolved using Finite Volume approach in the time
 * interval 0 <= t <= 10, and compared with the analytical solution.
 *
 * Keywords: diffusion equation, 2D, finite volume method.
 */

#include <iostream>
#include <chrono>
#include <cmath>

#include "../../src/solver.h"  // the main header of dae-cpp library solver

#include "diffusion_2d_RHS.h"  // RHS of the problem

namespace dae = daecpp;  // A shortcut to dae-cpp library namespace

// python3 + numpy + matplotlib should be installed in order to enable plotting
// #define PLOTTING

#ifdef PLOTTING
#include "../../src/external/matplotlib-cpp/matplotlibcpp.h"
namespace plt = matplotlibcpp;
#endif

// To compare dae-cpp solution with the analytical solution
int solution_check(dae::state_type &x);

/*
 * MAIN FUNCTION
 * =============================================================================
 * Returns '0' if solution comparison is OK or '1' if solution error is above
 * acceptable tolerance.
 */
int main()
{
    // These parameters can be obtained from a parameter file or as command line
    // options. Here for simplicity we define them as constants.
    const MKL_INT N  = 10;    // Number of cells along one axis
    const double  D  = 1.0;   // Diffusion coefficient
    const double  t1 = 10.0;  // Integration time (0 < t < t1)

    std::cout << "N = " << N << "; D = " << D << "; t = " << t1 << '\n';

    using clock     = std::chrono::high_resolution_clock;
    using time_unit = std::chrono::milliseconds;

    // Define state vectors. Here N*N is the total number of the equations.
    // We are going to carry out two independent simulations: with analytical
    // Jacobian and with numerically estimated one, hence two vectors.
    dae::state_type x1(N * N);
    dae::state_type x2(N * N);

    // Initial conditions
    for(MKL_INT i = 0; i < N * N; i++)
    {
        x1[i] = 0.0;
    }
    x1[N * N / 2] = N * N;  // 1/(h*h)

    // x1 and x2 will be overwritten by the solver
    x2 = x1;

    // Set up the RHS of the problem.
    // Class MyRHS inherits abstract RHS class from dae-cpp library.
    MyRHS rhs(N, D);

    // We can override Jacobian class from dae-cpp library
    // and provide analytical Jacobian
    //    MyJacobian jac(rhs, p);

    // Set up the Mass Matrix of the problem. In this case this matrix is
    // identity, so we can use a helper class provided by dae-cpp library.
    dae::MassMatrixIdentity mass(N * N);

    // Create an instance of the solver options and update some of the solver
    // parameters defined in solver_options.h
    dae::SolverOptions opt;

    opt.dt_init         = 0.1;    // Change initial time step;
    opt.fact_every_iter = false;  // Gain some speed. The matrices will be
                                  // factorized only once each time step.

    // Create an instance of the solver with particular RHS, Mass matrix,
    // Jacobian and solver options
    //    dae::Solver solve(rhs, jac, mass, opt, t1);

    // Now we are ready to solve the set of DAEs
    std::cout << "\nStarting DAE solver...\n";

    {
        auto tic0 = clock::now();
        //        solve(x1);
        auto tic1 = clock::now();

        std::cout
            << "Solver execution time: "
            << std::chrono::duration_cast<time_unit>(tic1 - tic0).count() /
                   1000.0
            << " sec." << '\n';
    }

    // Compare solution with an alternative solver (e.g. MATLAB)
    int check_result = solution_check(x1);

    // If we don't provide analytical Jacobian we need to estimate it
    // with a given tolerance:
    dae::Jacobian jac_est(rhs, opt.atol);

    // Create a new instance of the solver for estimated Jacobian
    dae::Solver solve_slow(rhs, jac_est, mass, opt, t1);

    // Solve the set of DAEs again
    std::cout << "\nStarting DAE solver with estimated Jacobian...\n";

    {
        auto tic0 = clock::now();
        solve_slow(x2);
        auto tic1 = clock::now();

        std::cout
            << "Solver execution time: "
            << std::chrono::duration_cast<time_unit>(tic1 - tic0).count() /
                   1000.0
            << " sec." << '\n';
    }

    // Compare solution
    check_result += solution_check(x2);

    // Plot the results
#ifdef PLOTTING
    dae::state_type x_axis(N), P(N), Phi(N);

    for(MKL_INT i = 0; i < N; i++)
    {
        x_axis[i] = (double)(i) / (N - 1);
        P[i]      = x1[i];
        Phi[i]    = x1[i + N];
    }

    plt::figure();
    plt::figure_size(800, 600);
    plt::named_plot("P(x)", x_axis, P, "b-");
    plt::named_plot("Phi(x)", x_axis, Phi, "r-");
    plt::xlabel("x");
    plt::ylabel("P and Phi");
    plt::xlim(0.0, 1.0);
    plt::grid(true);
    plt::legend();

    // Save figure
    const char *filename = "figure.png";
    std::cout << "Saving result to " << filename << "...\n";
    plt::save(filename);
#endif

    if(check_result)
        std::cout << "...Test FAILED\n\n";
    else
        std::cout << "...done\n\n";

    return check_result;
}

/*
 * Returns '0' if solution comparison is OK or '1' if the error is above
 * acceptable tolerance
 */
int solution_check(dae::state_type &x)
{
    std::cout << "Solution check:\n";

    const MKL_INT N = (MKL_INT)(x.size()) / 2;

    const int N_sol = 9;

    double sol[N_sol];

    // MATLAB ode15s solution at different x, Finite Elements, N = 4000 points
    const double ode15s_MATLAB[N_sol] = {19.9949, 2.72523,  0.382148,
                                         -10.0,   -6.04056, -2.08970,
                                         1.90021, 5.93011,  10.0};

    // dae-cpp solution at the same coordinates x:
    // clang-format off
    sol[0] = x[0];                                           // P(x = 0)
    sol[1] = x[(N-1)/10] * 0.1 + x[(N-1)/10+1] * 0.9;        // P(x = 0.1)
    sol[2] = x[(N-1)/5] * 0.2 + x[(N-1)/5+1] * 0.8;          // P(x = 0.2)
    sol[3] = x[N];                                           // Phi(x = 0)
    sol[4] = x[N+(N-1)/5*1] * 0.2 + x[N+(N-1)/5*1+1] * 0.8;  // Phi(x = 0.2)
    sol[5] = x[N+(N-1)/5*2] * 0.4 + x[N+(N-1)/5*2+1] * 0.6;  // Phi(x = 0.4)
    sol[6] = x[N+(N-1)/5*3] * 0.6 + x[N+(N-1)/5*3+1] * 0.4;  // Phi(x = 0.6)
    sol[7] = x[N+(N-1)/5*4] * 0.8 + x[N+(N-1)/5*4+1] * 0.2;  // Phi(x = 0.8)
    sol[8] = x[2*N-1];                                       // Phi(x = 1)
    // clang-format on

    std::cout << "  MATLAB ode15s\t<->  dae-cpp\t(rel. error)\n";

    double err_max = 0;

    for(int i = 0; i < N_sol; i++)
    {
        double error = (sol[i] - ode15s_MATLAB[i]) / ode15s_MATLAB[i] * 100.0;

        if(std::fabs(error) > err_max)
        {
            err_max = std::fabs(error);
        }

        std::cout << "      " << ode15s_MATLAB[i] << "\t<->  " << sol[i]
                  << " \t(" << error << "%)\n";
    }

    std::cout << "Maximum relative error: " << err_max << "%\n";

    if(err_max < 1.0)
        return 0;
    else
        return 1;
}