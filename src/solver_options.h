/*
 * Contains public solver parameters
 */

#pragma once

#include "typedefs.h"

namespace daecpp_namespace_name
{

#define BDF_MAX_ORDER 6  // Max BDF scheme order currently implemented

class SolverOptions
{
public:
    // You can control the parallel execution of the solver by explicitly
    // setting the MKL_NUM_THREADS environment variable. If fewer OpenMP threads
    // are available than specified, the execution may slow down instead of
    // speeding up. If MKL_NUM_THREADS is not defined, then the solver uses all
    // available processors.
    // In order to control the parallel execution of the numerical Jacobian, use
    // OMP_NUM_THREADS environment variable (uses all cores by default).

    // Perform Jacobian update, Reordering, Symbolic and Numerical Factorization
    // every Newton iteration. Changing to 'false' can increase speed but also
    // can lead to instability.
    bool fact_every_iter = true;

    // If fact_every_iter = false, update Jacobian every fact_iter Newton
    // iterations
    int fact_iter = 15;

    // Order of BDF implicit numerical integration method:
    // 1 - first order BDF, 2 - BDF-2, ..., 6 - BDF-6
    // Default is BDF-2 since it fully supports variable time stepping
    int bdf_order = 2;

    // Time stepping algorithm:
    // 1 - Stability-based Simple Adaptive Time Stepping (S-SATS),
    // 2 - Variability-based Simple Adaptive Time Stepping (V-SATS) from
    // https://www.sciencedirect.com/science/article/pii/S0377042705005534
    int time_stepping = 2;

    // Maximum number of Newton iterations. If the Newton method fails to
    // converge after max_Newton_iter iterations, the solver reduces time step
    // and tries to make the current step again.
    int max_Newton_iter = 15;

#ifdef DAE_SINGLE
    double atol      = 1.0e-3;  // Absolute tolerance for the Newton algorithm
    double rtol      = 1.0e-6;  // Relative tolerance for the Newton algorithm
    double dt_eps_m  = 1.0e-6;  // The order of the rounding unit
    double value_max = 1.0e20;  // Solution shouldn't be higher than this
#else
    double atol      = 1.0e-6;   // Absolute tolerance for the Newton algorithm
    double rtol      = 1.0e-6;   // Relative tolerance for the Newton algorithm
    double dt_eps_m  = 1.0e-14;  // The order of the rounding unit
    double value_max = 1.0e100;  // Solution shouldn't be higher than this
#endif

    // Initial time step
    double dt_init = 0.1;

    // Initial integration time t0
    double t0 = 0.0;

    // Minimum time step
    double dt_min = dt_eps_m;

    // Maximum time step
    double dt_max = 1.0 / dt_eps_m;

    // Verbosity level of the solver:
    // 0 - silent, 1 - basic information, 2 - time stepping info,
    // 3 - time integrator info, 4 - all info
    int verbosity = 2;

    // Simple Adaptive Time Stepping options
    int dt_increase_threshold = 2;    // Time step amplification threshold
                                      // (S-SATS only)
    int dt_decrease_threshold = 7;    // Time step reduction threshold
                                      // (S-SATS only)
    double dt_increase_factor = 2.0;  // Time step amplification factor
    double dt_decrease_factor = 2.0;  // Time step reduction factor
    double dt_eta_min = 0.05;  // Monitor function lower threshold (V-SATS only)
    double dt_eta_max = 0.5;  // Monitor function higher threshold (V-SATS only)

    // Try to roll back and reduce the time step if Newton iterations diverged.
    // Otherwise stop with error message.
    bool redo_newton = false;

    // If Newton method fails to converge within 'max_Newton_iter' iterations
    // 'newton_failed_attempts' times in a row, the solver will try to update
    // Jacobian every single iteration next time step.
    int newton_failed_attempts = 3;

    // 1 - V-SATS will use NORM_infinity to estimate solution variability,
    // 2 - V-SATS will use NORM_2 (default)
    int vsats_norm = 2;

    // Intel MKL PARDISO parameters (iparam). More about iparam:
    // https://software.intel.com/en-us/mkl-developer-reference-c-pardiso-iparm-parameter

    // iparm[3]: Controls preconditioned CGS.
    // 0  - The factorization is always computed as required.
    // 31 - LU-preconditioned CGS iteration with a stopping criterion of 1.0E-3.
    // 61 - LU-preconditioned CGS iteration with a stopping criterion of 1.0E-6.
    MKL_INT preconditioned_CGS = 0;

    // iparm[7]: Maximum number of iterative refinement steps.
    // 0  - The solver automatically performs two steps of iterative refinement.
    // >0 - Maximum number of iterative refinement steps that the solver
    //      performs.
    MKL_INT refinement_steps = 2;

    // iparm[12]: Maximum weighted matching algorithm.
    // 0 - OFF
    // 1 - ON (may conflict with parallel_fact_control)
    MKL_INT matching_alg = 1;

    // iparm[23]: Parallel factorization control.
    // 0  - Intel MKL PARDISO uses the classic algorithm for factorization.
    // 1  - Two-level factorization algorithm. This algorithm generally improves
    //      scalability in case of parallel factorization on many OpenMP threads
    //      (more than eight).
    // 10 - Improved two-level factorization algorithm.
    MKL_INT parallel_fact_control = 10;

    // iparm[24]: Parallel forward/backward solve control.
    // 0 - Intel MKL PARDISO uses the sequential forward and backward solve.
    // 1 - Parallel algorithm for the solve step (in-core mode only).
    // 2 - If there is only one right hand side, Intel MKL PARDISO performs
    //     parallelization. If there are multiple right hand sides, PARDISO
    //     performs parallel forward and backward substitution via matrix.
    MKL_INT parallel_solve_control = 0;

    SolverOptions() = default;

    virtual ~SolverOptions() = default;

    // Initialises Intel MKL PARDISO parameters (iparam) array
    void set_iparm_for_pardiso(MKL_INT *iparm);

    // Checks correctness of the solver parameters
    // TODO: should return error code
    void check_options();
};

}  // namespace daecpp_namespace_name
