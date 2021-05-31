/*
 * Mass matrix classes
 */

#pragma once

#include "typedefs.h"

namespace daecpp_namespace_name
{

/*
 * Parent class. This class is abstract and must be inherited.
 */
class MassMatrix
{
    /*
     * Sparse matrix converter from simple three-array format to Intel MKL
     * three array format.
     * Input: matrix holder M with simple three-array format
     * Output: matrix holder M with Intel MKL three-array format
     */
    void m_matrix_converter(daecpp::sparse_matrix_holder &M);

public:
    virtual ~MassMatrix() {}

    /*
     * The matrix should be defined in sparse format,
     * see three array sparse format decription on
     * https://software.intel.com/en-us/mkl-developer-reference-c-sparse-blas-csr-matrix-storage-format
     *
     * The Mass matrix is static, i.e. it does not depend on time t and
     * the vector x
     *
     * This function is pure virtual and must be overriden.
     */
    virtual void operator()(sparse_matrix_holder &M) = 0;

    /*
     * Helper function to write the Mass matrix to a file
     */
    void dump();
};

/*
 * Helper child class to create identity Mass matrix of size N
 */
class MassMatrixIdentity : public MassMatrix
{
    const MKL_INT m_N;

public:
    MassMatrixIdentity(const MKL_INT N) : daecpp::MassMatrix(), m_N(N) {}

    void operator()(daecpp::sparse_matrix_holder &M);
};

}  // namespace daecpp_namespace_name
