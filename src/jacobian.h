/*
 * Calculates Jacobian matrix.
 * If not overridden by a user, performs numerical differentiation of the RHS
 * with the given tolerance to estimate numerical Jacobian matrix.
 */

#pragma once

#include "typedefs.h"
#include "RHS.h"

namespace daecpp_namespace_name
{

class Jacobian
{
    RHS &m_rhs;

#ifdef DAE_SINGLE
    const double m_tol = 1.0e-3;  // Default tolerance
    const double m_eps = 1.0e-6;  // The order of the rounding unit
#else
    const double m_tol = 1.0e-6;   // Default tolerance
    const double m_eps = 1.0e-13;  // The order of the rounding unit
#endif

    std::size_t m_dump_file_counter    = 0;
    std::size_t m_compare_file_counter = 0;

    int m_jac_type = 0;  // This will be changed to 1
                         // if numerical Jacobian is used

    /*
     * Sparse matrix converter from simple three-array format to Intel MKL
     * three array format.
     * Input: matrix holder M with simple three-array format
     * Output: matrix holder M with Intel MKL three-array format
     */
    void m_matrix_converter(daecpp::sparse_matrix_holder &M);

public:
    explicit Jacobian(RHS &rhs) : m_rhs(rhs) {}

    Jacobian(RHS &rhs, const double tol) : m_rhs(rhs), m_tol(tol)
    {
        // TODO: Check user's tol parameter. Too small tol may lead to crash.
    }

    virtual ~Jacobian() {}

    /*
     * Can be overriden to provide analytical Jacobian
     */
    virtual void operator()(sparse_matrix_holder &J, const state_type &x,
                            const double t);

    /*
     * Helper function to show Jacobian structure on screen (in sparse format)
     */
    void print(const state_type &x, const double t);

    /*
     * Helper function to write Jacobian matrix to a file (in dense format)
     */
    void dump(const state_type &x, const double t);

    /*
     * Helper function to compare two Jacobians and write the differences.
     * Comparison will be made with the external Jacobian jac (usually,
     * numerical Jacobian) using vector x at time t with the given tolerance.
     * Returns the number of differences found.
     */
    int compare(Jacobian jac, const state_type &x, const double t,
                const double tol);
};

}  // namespace daecpp_namespace_name
