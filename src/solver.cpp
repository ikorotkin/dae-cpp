/*
 * The main solver
 */

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <chrono>

#include <mkl_types.h>
#include <mkl_pardiso.h>

#include "solver.h"

namespace daecpp_namespace_name
{

/*
 * The solver constructor
 */
Solver::Solver(RHS &rhs, Jacobian &jac, MassMatrix &mass, SolverOptions &opt)
    : m_rhs(rhs), m_jac(jac), m_mass(mass), m_opt(opt)
{
    // Initial output
    if(m_opt.verbosity > 1)
    {
        std::cout << "Float precision:     " << 8 * sizeof(float_type)
                  << " bit\n";
        std::cout << "Integer precision:   " << 8 * sizeof(MKL_INT) << " bit\n";
        std::cout << "Numerical algorithm: BDF-" << m_opt.bdf_order
                  << std::endl;
    }

    // Initialises the internal solver memory pointer. This is only
    // necessary for the FIRST call of the PARDISO solver.
    for(MKL_INT i = 0; i < 64; i++)
    {
        m_pt[i] = 0;
    }

    // Initialises current integration time
    m_iterator_state.t = opt.t0;

    // Sets initial time steps
    m_iterator_state.dt[0] = opt.dt_init;
    m_iterator_state.dt[1] = 0.0;

    // Solver starts the first time step using BDF-1 method
    // (since it doesn't have enough history yet)
    m_iterator_state.current_scheme = 1;

    // Initialises the time integrator
    m_ti = new TimeIntegrator(rhs, jac, mass, opt);

    // Initialises solution history for the time integrator
    m_x_prev.resize(opt.bdf_order);

    // Loads Intel MKL PARDISO iparm parameter from solver_options class
    opt.set_iparm_for_pardiso(m_iparm);
}

/*
 * The main solver
 * =============================================================================
 */
int Solver::operator()(state_type &x, double &t1)
{
    // Set system size
    m_size = (MKL_INT)(x.size());

    // Check user-defined solver options
    m_opt.check_options();

    // We don't need to do anything if t1 == t0. Return initial conditions.
    if(t1 == m_iterator_state.t)
        return 0;

    // Assert t1 > t0
    if(t1 < m_iterator_state.t)
    {
        std::cout << "ERROR: Integration time t1 = " << t1
                  << " cannot be less than the initial time t0 = "
                  << m_iterator_state.t << std::endl;
        return 1;
    }

    // Assert dt > dt_min
    if(m_iterator_state.dt[0] < m_opt.dt_min)
        m_iterator_state.dt[0] = m_opt.dt_init;

    // Check initial time steps
    m_iterator_state.dt_eval =
        (m_iterator_state.dt[0] > (t1 - m_iterator_state.t))
            ? (t1 - m_iterator_state.t)
            : m_iterator_state.dt[0];
    m_iterator_state.dt[0] = m_iterator_state.dt_eval;

    // Roll back to BDF-2 if the time step has changed
    if((m_iterator_state.dt[0] != m_iterator_state.dt[1]) &&
       (m_iterator_state.current_scheme > 1))
    {
        m_iterator_state.current_scheme = 2;
    }

    // Initialize the time integrator state structure.
    m_iterator_state.step_counter_local = 0;
    m_iterator_state.final_time_step    = false;

    // Initial output
    if(m_opt.verbosity > 0)
    {
        std::cout << "Number of equations: " << m_size << std::endl;
    }

    // Reserve memory for the solution history. This will be done only once
    if(m_x_prev[0].size() == 0)
    {
        for(int i = 0; i < m_opt.bdf_order; i++)
            m_x_prev[i].resize(m_size);
    }

    // Copy current state vector into the history vector
    m_x_prev[0] = x;

    // Full Jacobian matrix holder
    sparse_matrix_holder J;

    // Full RHS vector
    state_type b(m_size);

    // Solution vector used for Newton iterations
    state_type xk(m_size);

    // Reset PARDISO pointers
    m_mkl_b = b.data();
    m_mkl_x = xk.data();

    // Memory control variables
    int    peak_mem1 = 0, peak_mem2 = 0, peak_mem3 = 0;
    double total_peak_mem = 0.0;

    // Reset time integrator timers
    m_ti->reset_jac_time();
    m_ti->reset_rhs_time();

    // Reset linear algebra solver timer
    double lin_alg_time = 0.0;

    // Initialise clock
    using clock     = std::chrono::high_resolution_clock;
    using time_unit = std::chrono::microseconds;

    // Counts linear solver calls
    std::size_t calls = 0;

    // Counts how many times the Newton iterator failed to converge within
    // max_Newton_iter iterations in a row.
    int n_iter_failed = 0;

    if(m_opt.verbosity == 1)
    {
        std::cout << "Calculating...";
        std::cout.flush();
    }

    /*
     * Start the solver
     * =========================================================================
     */

    // Timer starts here
    auto tic0 = clock::now();

    while(m_iterator_state.t < (t1 + m_iterator_state.dt[0] * 0.5))
    {
        m_iterator_state.t += m_iterator_state.dt[0];  // Time step lapse

        m_iterator_state.step_counter_local++;
        m_steps++;

        if(m_opt.verbosity > 1)
        {
            std::cout << std::left;
            std::cout << "\nStep " << std::setw(7) << m_steps
                      << " :: t = " << std::setw(12) << m_iterator_state.t
                      << " :: ";
            // std::cout.flush();  // This degrades performance in some cases
        }

        if(m_opt.verbosity > 2)
        {
            std::cout << "BDF-" << m_iterator_state.current_scheme << ": ";
            std::cout << "dt=" << m_iterator_state.dt[0]
                      << ", dt_prev=" << m_iterator_state.dt[1] << ": ";
        }

        m_ti->set_scheme(m_iterator_state.current_scheme);

        // Can be set to true by the solver if it fails to converge
        bool fact_every_iter = (n_iter_failed >= m_opt.newton_failed_attempts)
                                   ? true
                                   : m_opt.fact_every_iter;

        int iter;  // Loop index. We need this value later

        for(iter = 0; iter < m_opt.max_Newton_iter; iter++)
        {
            // Reordering, Symbolic and Numerical Factorization
            if(fact_every_iter || iter == 0 || !(iter % m_opt.fact_iter))
            {
                // Time Integrator with updated Jacobian
                m_ti->integrate(J, b, x, m_x_prev, m_iterator_state.t,
                                m_iterator_state.dt, true);

                // Jacobian can change its size and can be re-allocated.
                // Catch up new array addresses.
                m_mkl_a = J.A.data();
                m_ia    = J.ia.data();
                m_ja    = J.ja.data();

                auto tic_phase1 = clock::now();

                // PHASE 1.
                // Reordering and Symbolic Factorization. This step also
                // allocates all memory that is necessary for the factorization
                m_phase = 11;
                PARDISO(m_pt, &m_maxfct, &m_mnum, &m_mtype, &m_phase, &m_size,
                        m_mkl_a, m_ia, m_ja, &m_idum, &m_nrhs, m_iparm,
                        &m_msglvl, &m_ddum, &m_ddum, &m_error);

                if(m_opt.verbosity > 3)
                {
                    if(m_iparm[14] > peak_mem1 || m_iparm[15] > peak_mem2 ||
                       m_iparm[16] > peak_mem3)
                    {
                        peak_mem1 = m_iparm[14];
                        peak_mem2 = m_iparm[15];
                        peak_mem3 = m_iparm[16];

                        std::cout << "\nPeak memory on symbolic factorization: "
                                  << "              "
                                  << (double)peak_mem1 / 1024.0 << " Mb";
                        std::cout
                            << "\nPermanent memory on symbolic factorization: "
                            << "         " << (double)peak_mem2 / 1024.0
                            << " Mb";
                        std::cout << "\nPeak memory on numerical factorization "
                                     "and solution: "
                                  << (double)peak_mem3 / 1024.0 << " Mb"
                                  << std::endl;
                    }
                }
                if(m_opt.verbosity > 0)
                {
                    total_peak_mem =
                        (double)(m_iparm[14] + m_iparm[16]) / 1024.0;
                }

                if(m_error != 0)
                {
                    std::cout << "\nERROR during symbolic factorization...\n";
                    m_check_pardiso_error(m_error);
                    return 11;
                }

                // PHASE 2.
                // Numerical factorization
                m_phase = 22;
                PARDISO(m_pt, &m_maxfct, &m_mnum, &m_mtype, &m_phase, &m_size,
                        m_mkl_a, m_ia, m_ja, &m_idum, &m_nrhs, m_iparm,
                        &m_msglvl, &m_ddum, &m_ddum, &m_error);

                if(m_error != 0)
                {
                    std::cout << "\nERROR during numerical factorization...\n";
                    m_check_pardiso_error(m_error);
                    return 22;
                }

                lin_alg_time += std::chrono::duration_cast<time_unit>(
                                    clock::now() - tic_phase1)
                                    .count() *
                                1e-6;
            }
            else
            {
                // Time Integrator with the previous Jacobian
                m_ti->integrate(J, b, x, m_x_prev, m_iterator_state.t,
                                m_iterator_state.dt, false);
            }

            auto tic_phase3 = clock::now();

            // PHASE 3.
            // Back substitution and iterative refinement
            m_phase = 33;
            PARDISO(m_pt, &m_maxfct, &m_mnum, &m_mtype, &m_phase, &m_size,
                    m_mkl_a, m_ia, m_ja, &m_idum, &m_nrhs, m_iparm, &m_msglvl,
                    m_mkl_b, m_mkl_x, &m_error);

            if(m_error != 0)
            {
                std::cout << "\nERROR during solution...\n";
                m_check_pardiso_error(m_error);
                return 33;
            }

            lin_alg_time +=
                std::chrono::duration_cast<time_unit>(clock::now() - tic_phase3)
                    .count() *
                1e-6;

            calls++;
            m_calls++;

            bool is_converged = true;

            for(MKL_INT i = 0; i < m_size; i++)
            {
                double adiff = std::abs(m_mkl_x[i]);

                if(adiff > m_opt.value_max || std::isnan(m_mkl_x[i]))
                {
                    if(!m_opt.redo_newton || m_opt.verbosity > 1)
                    {
                        std::cout << "\nNewton iterations diverged. "
                                  << "Review the solver options.\n";
                    }

                    if(m_opt.redo_newton)
                    {
                        if(m_opt.verbosity > 1)
                            std::cout << "Trying to recover...\n";
                        iter = m_opt.max_Newton_iter;
                        break;
                    }
                    else
                    {
                        return 2;  // Newton iterations diverged
                    }
                }

                if(is_converged)
                {
                    if(x[i] != 0.0)
                    {
                        double rdiff = adiff / std::abs(x[i]);
                        if(adiff > m_opt.atol && rdiff > m_opt.rtol)
                        {
                            is_converged = false;
                        }
                    }
                    else
                    {
                        if(adiff > m_opt.atol)
                        {
                            is_converged = false;
                        }
                    }
                }

                x[i] -= m_mkl_x[i];
            }

            if(m_opt.verbosity > 1)
            {
                std::cout << "#";
                // std::cout.flush();  // This degrades performance
            }

            if(is_converged)
            {
                break;
            }

        }  // for iter

        // Newton iterator failed to converge within max_Newton_iter iterations.
        // Trying to reduce the time step.
        if(iter >= m_opt.max_Newton_iter)
        {
            n_iter_failed++;
            if(m_opt.verbosity > 1)
                std::cout << " <- redo";
            if(m_reset_ti_state(x, m_x_prev) == -1)
                return 3;  // The time step was reduced to dt_min,
                           // but the solver still cannot converge
            continue;
        }
        else
        {
            n_iter_failed = 0;
        }

        // The solver has reached the target time t1 or the stop condition
        // triggered.
        if(m_iterator_state.final_time_step ||
           m_rhs.stop_condition(x, m_iterator_state.t))
        {
            break;
        }

        // Adaptive time stepping algorithm
        int status = m_adaptive_time_stepping(x, m_x_prev, iter);
        if(status < 0)
            return 4;  // The algorithm failed to converge
        else if(status > 0)
            continue;  // Re-run the current time step

        // Call Observer to provide a user with intermediate results
        observer(x, m_iterator_state.t);

        // Check the time step does not exceed dt_max
        if(m_iterator_state.dt_eval > m_opt.dt_max)
        {
            m_iterator_state.dt_eval = m_opt.dt_max;
        }

        // Looks like the solver has reached the target time t1
        if((m_iterator_state.t + m_iterator_state.dt_eval) >=
           (t1 - m_opt.dt_min))
        {
            // Adjust the last time step size
            double dt_max = t1 - m_iterator_state.t;

            if(std::abs(dt_max) < m_opt.dt_eps_m)  // Should never happen
            {
                break;  // The solver has reached t1
            }
            else
            {
                m_iterator_state.final_time_step = true;
                m_iterator_state.dt_eval         = dt_max;
            }
        }

        // Rewrite solution history
        for(int d = m_opt.bdf_order - 1; d > 0; d--)
        {
            m_x_prev[d] = m_x_prev[d - 1];
        }
        m_x_prev[0] = x;

        // Update the BDF order
        if(m_iterator_state.current_scheme < m_opt.bdf_order)
        {
            if((m_iterator_state.current_scheme == 1) ||
               (m_iterator_state.dt[0] == m_iterator_state.dt[1]))
            {
                m_iterator_state.current_scheme++;
            }
        }

        // Update time step history
        m_iterator_state.dt[1] = m_iterator_state.dt[0];
        m_iterator_state.dt[0] = m_iterator_state.dt_eval;

        // Roll back to BDF-2 if the time step has changed
        if((m_iterator_state.dt[0] != m_iterator_state.dt[1]) &&
           (m_opt.bdf_order > 1))
        {
            m_iterator_state.current_scheme = 2;
        }

    }  // while t

    // Stop timer
    auto tic1 = clock::now();

    // Update solution time
    t1 = m_iterator_state.t;

    // Catch up the last time step
    observer(x, m_iterator_state.t);

    // Update solution history
    for(int d = m_opt.bdf_order - 1; d > 0; d--)
    {
        m_x_prev[d] = m_x_prev[d - 1];
    }
    m_x_prev[0] = x;

    // Restore the previous time step size
    m_iterator_state.dt_eval = m_iterator_state.dt[1];
    m_iterator_state.dt[1]   = m_iterator_state.dt[0];
    m_iterator_state.dt[0]   = m_iterator_state.dt_eval;

    // Final output
    if(m_opt.verbosity > 0)
    {
        double solver_time =
            std::chrono::duration_cast<time_unit>(tic1 - tic0).count() * 1e-6;
        double jac_time   = m_ti->get_jac_time();
        double rhs_time   = m_ti->get_rhs_time();
        double other_time = solver_time - (lin_alg_time + rhs_time + jac_time);
        double jac_time_rel     = jac_time / solver_time * 100.0;
        double rhs_time_rel     = rhs_time / solver_time * 100.0;
        double lin_alg_time_rel = lin_alg_time / solver_time * 100.0;
        double other_time_rel =
            100.0 - (lin_alg_time_rel + rhs_time_rel + jac_time_rel);

        m_timer_lin += lin_alg_time;
        m_timer_rhs += rhs_time;
        m_timer_jac += jac_time;
        m_timer_tot += solver_time;

        double timer_other =
            m_timer_tot - (m_timer_lin + m_timer_rhs + m_timer_jac);

        std::stringstream ss;

        ss << std::fixed << std::setprecision(3);
        ss << "\nLinear algebra solver calls:       " << calls;
        if(m_dae_solver_calls)
            ss << " (" << m_calls << " total)";
        ss << "\nPeak memory for the linear solver: " << total_peak_mem
           << " Mb\n";
        ss << "Time spent:\n  by linear algebra solver:     " << lin_alg_time
           << " sec. (" << lin_alg_time_rel << "%)";
        if(m_dae_solver_calls)
            ss << " --> " << m_timer_lin << " sec. ("
               << m_timer_lin / m_timer_tot * 100.0 << "%)";
        ss << "\n  to calculate the RHS:         " << rhs_time << " sec. ("
           << rhs_time_rel << "%)";
        if(m_dae_solver_calls)
            ss << " --> " << m_timer_rhs << " sec. ("
               << m_timer_rhs / m_timer_tot * 100.0 << "%)";
        ss << "\n  to calculate Jacobian:        " << jac_time << " sec. ("
           << jac_time_rel << "%)";
        if(m_dae_solver_calls)
            ss << " --> " << m_timer_jac << " sec. ("
               << m_timer_jac / m_timer_tot * 100.0 << "%)";
        ss << "\n  other calculations:           " << other_time << " sec. ("
           << other_time_rel << "%)";
        if(m_dae_solver_calls)
            ss << " --> " << timer_other << " sec. ("
               << timer_other / m_timer_tot * 100.0 << "%)";
        ss << "\nTotal time spent by the solver: " << solver_time
           << " sec. (100.0%)";
        if(m_dae_solver_calls)
            ss << " --> " << m_timer_tot << " sec. (100.0%)";
        ss << "\n\n";

        std::cout << ss.str();
        std::cout.flush();
    }

    // Success
    m_dae_solver_calls++;
    return 0;
}

/*
 * Releases memory
 */
Solver::~Solver()
{
    m_phase = -1;  // Termination and release of memory
    PARDISO(m_pt, &m_maxfct, &m_mnum, &m_mtype, &m_phase, &m_size, &m_ddum,
            m_ia, m_ja, &m_idum, &m_nrhs, m_iparm, &m_msglvl, &m_ddum, &m_ddum,
            &m_error);

    delete m_ti;
}

/*
 * Checks PARDISO solver error messages
 */
void Solver::m_check_pardiso_error(MKL_INT err)
{
    if(!err)
    {
        return;  // All good
    }
    else
    {
        std::cout << "Linear algebra solver ERROR: ";

        switch(err)
        {
        case -1:
            std::cout << "input inconsistent.\n";
            break;
        case -2:
            std::cout << "not enough memory.\n";
            break;
        case -3:
            std::cout << "reordering problem.\n";
            break;
        case -4:
            std::cout << "zero pivot, numerical factorization or iterative "
                         "refinement problem.\n";
            break;
        case -5:
            std::cout << "unclassified (internal) error.\n";
            break;
        case -6:
            std::cout << "reordering failed(matrix types 11 and 13 only).\n";
            break;
        case -7:
            std::cout << "diagonal matrix is singular.\n";
            break;
        case -8:
            std::cout << "32-bit integer overflow problem.\n";
            break;
        case -9:
            std::cout << "not enough memory for OOC.\n";
            break;
        case -10:
            std::cout << "error opening OOC files.\n";
            break;
        case -11:
            std::cout << "read/write error with OOC files.\n";
            break;
        case -12:
            std::cout << "pardiso_64 called from 32-bit library.\n";
            break;
        case -13:
            std::cout
                << "interrupted by the (user-defined) mkl_progress function.\n";
            break;
        case -15:
            std::cout << "internal error which can appear for iparm[23] "
                         "(parallel_fact_control) = 10 and iparm[12] = 1. Try "
                         "switch matching off (set iparm[12] = 0 and rerun).\n";
            break;
        default:
            std::cout << "Unknown error.\n";
            break;
        }

        std::cout << std::endl;
    }
}

}  // namespace daecpp_namespace_name
