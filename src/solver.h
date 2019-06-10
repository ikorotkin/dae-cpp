/*
 * The main solver class definition
 */

#pragma once

#include "typedefs.h"
#include "RHS.h"
#include "jacobian.h"
#include "mass_matrix.h"
#include "solver_options.h"

namespace daecpp_namespace_name
{

class Solver
{
    RHS &m_rhs;  // RHS

    Jacobian &m_jac;  // Jacobian matrix

    MassMatrix &m_mass;  // Mass matrix

    SolverOptions &m_opt;  // Solver options

    MKL_INT m_size;  // System size

    size_t m_steps = 0;  // Internal time iteration counter
    size_t m_calls = 0;  // Internal solver calls counter

    // Intel MKL PARDISO control parameters
    MKL_INT m_phase;        // Current phase of the solver
    MKL_INT m_maxfct = 1;   // Maximum number of numerical factorizations
    MKL_INT m_mnum   = 1;   // Which factorization to use
    MKL_INT m_mtype  = 11;  // Real unsymmetric matrix
    MKL_INT m_nrhs   = 1;   // Number of right hand sides
    MKL_INT m_msglvl = 0;   // Print statistical information
    MKL_INT m_error  = 0;   // Error flag

    // Intel MKL PARDISO sparse matrix indeces
    MKL_INT *m_ia = nullptr;
    MKL_INT *m_ja = nullptr;

    // Intel MKL PARDISO vectors and sparse matrix non-zero elements
    float_type *m_mkl_a = nullptr;
    float_type *m_mkl_b = nullptr;
    float_type *m_mkl_x = nullptr;

    // Intel MKL PARDISO internal solver memory pointer pt,
    // 32-bit: int pt[64]; 64-bit: long int pt[64]
    // or void *pt[64] should be OK on both architectures
    void *m_pt[64];

    // Intel MKL PARDISO auxiliary variables
    double  m_ddum;  // Double dummy
    MKL_INT m_idum;  // Integer dummy

    // Intel MKL PARDISO iparm parameter
    MKL_INT m_iparm[64];

    // Updates time integrator scheme when the time step changes
    int m_reset_ti_scheme(SolverOptions &m_opt, const int step_counter);

    // Checks PARDISO solver error messages
    void m_check_pardiso_error(MKL_INT err);

public:
    /*
     * Receives user-defined RHS, Jacobian, Mass matrix and solver options
     */
    Solver(RHS &rhs, Jacobian &jac, MassMatrix &mass, SolverOptions &opt)
        : m_rhs(rhs), m_jac(jac), m_mass(mass), m_opt(opt)
    {
        // Initialise the internal solver memory pointer. This is only
        // necessary for the FIRST call of the PARDISO solver.
        for(MKL_INT i = 0; i < 64; i++)
        {
            m_pt[i] = 0;
        }
    }

    /*
     * Releases memory. Defined in solver.cpp.
     */
    ~Solver();

    /*
     * Integrates the system of DAEs on the interval t = [t0; t1] and returns
     * result in the array x. Parameter t0 can be overriden in the solver
     * options (t0 = 0 by default).
     * The data stored in x (initial conditions) will be overwritten.
     * Returns 0 in case of success or error code if integration failed.
     */
    int operator()(state_type &x, const double t1);

    /*
     * Virtual Observer. Called by the solver every time step.
     * Receives current solution vector and the current time.
     * Can be overriden by a user to get intermediate results (for example,
     * for plotting).
     */
    virtual void observer(state_type &x, const double t)
    {
        return;  // It does nothing by deafult
    }
};

}  // namespace daecpp_namespace_name
