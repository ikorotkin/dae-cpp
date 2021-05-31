/*
 * A set of helper functions to print on screen or write to files
 * Jacobian matrix, Mass matrix, the RHS for debugging purposes.
 */

#include <iostream>  // std::cout
#include <iomanip>   // std::setw etc.
#include <fstream>   // File output
#include <string>    // std::string, std::to_string
#include <cmath>     // std::abs

#include "RHS.h"
#include "mass_matrix.h"
#include "jacobian.h"

namespace daecpp_namespace_name
{

const char delimiter = '\t';  // Delimiter of columns in output text files

/*
 * Helper function to write the RHS vector to a file
 */
void RHS::dump(const state_type &x, const double t)
{
    std::cout << "RHS::dump()         -- INFO: Writing the RHS at time t = "
              << t << "...\n";

    const MKL_INT size = x.size();

    state_type f(size);  // the vector to be saved

    this->operator()(x, f, t);  // calls the RHS

    std::ofstream outFile;

    outFile.open("dump_RHS_" + std::to_string(m_dump_file_counter++) + ".txt");
    outFile << "t = " << t << ":\n";
    outFile << "i" << delimiter << "x[i]" << delimiter << "RHS[i]" << '\n';
    for(MKL_INT i = 0; i < size; i++)
        outFile << i << delimiter << x[i] << delimiter << f[i] << '\n';
    outFile.close();
}

/*
 * Helper function to write the Mass matrix to a file
 */
void MassMatrix::dump()
{
    std::cout << "MassMatrix::dump()  -- INFO: Writing the Mass matrix...\n";

    sparse_matrix_holder M;

    this->operator()(M);  // calls the Mass matrix operator

    m_matrix_converter(M);  // converts the matrix if it is in simple form

    const MKL_INT size =
        M.ia.size() - 1;  // derive the matrix size from ia index

    if(size > 10000)
    {
        std::cout << "MassMatrix::dump()  -- WARNING: the size of the Mass "
                     "matrix for writting is bigger than 10000x10000.\n";
    }

    std::ofstream outFile;

    MKL_INT ja = 0;

    outFile.open("dump_Mass_matrix.txt");  // Mass matrix is static - one file
    outFile << "i,j:";
    for(MKL_INT i = 0; i < size; i++)
    {
        outFile << delimiter << "i=" << i;
    }
    outFile << '\n';
    for(MKL_INT j = 0; j < size; j++)
    {
        MKL_INT ent = M.ia[j + 1] - M.ia[j];  // Number of entries in line j

        outFile << "j=" << j << delimiter;

        for(MKL_INT i = 0; i < size; i++)
        {
            if(M.ja[ja] == i)
            {
                outFile << M.A[ja++];
                if(!(--ent))
                    break;
            }
            outFile << delimiter;
        }

        outFile << '\n';
    }
    outFile.close();
}

/*
 * Helper function to write the Jacbian matrix to a file (in dense format)
 */
void Jacobian::dump(const state_type &x, const double t)
{
    std::cout << "Jacobian::dump()    -- INFO: ";

    sparse_matrix_holder M;

    this->operator()(M, x, t);  // calls the Jacobian matrix operator

    m_matrix_converter(M);  // converts the matrix if it is in simple form

    if(m_jac_type)
        std::cout << "Writing numerically estimated ";
    else
        std::cout << "Writing user-defined ";
    std::cout << "Jacobian matrix at time t = " << t << "...\n";

    const MKL_INT size =
        M.ia.size() - 1;  // derive the matrix size from ia index

    if(size > 10000)
    {
        std::cout << "Jacobian::dump()    -- WARNING: the size of the Jacobian "
                     "matrix for writting is bigger than 10000x10000.\n";
    }

    std::ofstream outFile;

    MKL_INT ja = 0;

    if(m_jac_type)
        outFile.open("dump_Jacobian_" + std::to_string(m_dump_file_counter++) +
                     "_numerical.txt");
    else
        outFile.open("dump_Jacobian_" + std::to_string(m_dump_file_counter++) +
                     ".txt");

    outFile << "t=" << t;
    for(MKL_INT i = 0; i < size; i++)
    {
        outFile << delimiter << "i=" << i;
    }
    outFile << '\n';
    for(MKL_INT j = 0; j < size; j++)
    {
        MKL_INT ent = M.ia[j + 1] - M.ia[j];  // Number of entries in line j

        outFile << "j=" << j << delimiter;

        for(MKL_INT i = 0; i < size; i++)
        {
            if(M.ja[ja] == i)
            {
                outFile << M.A[ja++];
                if(!(--ent))
                    break;
            }
            outFile << delimiter;
        }

        outFile << '\n';
    }
    outFile.close();
}

/*
 * Helper function to show Jacobian structure (in sparse format)
 */
void Jacobian::print(const state_type &x, const double t)
{
    if(x.size() > 1000)
    {
        std::cout << "\nJacobian::print() -- too much output. Skipped.\n";
        return;
    }

    sparse_matrix_holder J;

    this->operator()(J, x, t);

    m_matrix_converter(J);  // converts the matrix if it is in simple form

    std::cout << std::right;
    std::cout << "\nJacobian matrix at time t = " << t << ':';
    std::cout << "\n-----------------------------------------\n";
    std::cout << std::setw(7) << "i" << std::setw(16) << "J.A |"
              << std::setw(10) << "J.ja |" << std::setw(8) << "J.ia";
    std::cout << "\n-----------------------------------------\n";

    std::size_t size = (J.A.size() > J.ia.size()) ? J.A.size() : J.ia.size();

    for(std::size_t i = 0; i < size; i++)
    {
        std::cout << std::setw(7) << i << ": ";
        std::cout << std::setw(12);

        if(i < J.A.size())
            std::cout << J.A[i];
        else
            std::cout << ' ';

        std::cout << " | " << std::setw(7);

        if(i < J.ja.size())
            std::cout << J.ja[i];
        else
            std::cout << ' ';

        std::cout << " | ";

        if(i < J.ia.size())
            std::cout << std::setw(7) << J.ia[i];

        std::cout << std::endl;
    }
}

/*
 * Helper function to compare two Jacobians and write the difference
 */
int Jacobian::compare(Jacobian jac, const state_type &x, const double t,
                      const double tol)
{
    // std::cout << "Jacobian::compare() -- INFO: Trying to compare two "
    //              "Jacobians at time t = "
    //           << t << " and the tolerance tol = " << tol << "...\n";

    sparse_matrix_holder M, J;

    this->operator()(M, x, t);  // calls the Jacobian matrix operator
    jac(J, x, t);               // external Jacobian to compare with

    m_matrix_converter(M);  // converts the matrix M if it is in simple form
    m_matrix_converter(J);  // converts the matrix J if it is in simple form

    const MKL_INT size =
        M.ia.size() - 1;  // derive the matrix size from ia index

    if((std::size_t)(size) != (J.ia.size() - 1))
    {
        std::cout << "Jacobian::compare() -- ERROR: the sizes of the "
                     "matrices do not match ('ia' indexes are different).\n";
        return -1;
    }

    std::ofstream outFile;

    MKL_INT ja_M = 0;
    MKL_INT ja_J = 0;

    outFile.open("dump_Jacobians_compare_" +
                 std::to_string(m_compare_file_counter++) + ".txt");

    outFile << "List of differences in Jacobians for t = " << t
            << " and the tolerance tol = " << tol << ":\n";
    outFile << "i" << delimiter << "j" << delimiter << "Jac_original"
            << delimiter << "Jac_reference" << delimiter << "Rel_difference"
            << '\n';

    int ndiff = 0;  // counts differences

    for(MKL_INT j = 0; j < size; j++)
    {
        MKL_INT ent_M = M.ia[j + 1] - M.ia[j];
        MKL_INT ent_J = J.ia[j + 1] - J.ia[j];

        for(MKL_INT i = 0; i < size; i++)
        {
            double MA = 0.0;
            double JA = 0.0;
            double diff;

            if((!ent_M) && (!ent_J))
                break;

            if((M.ja[ja_M] == i) && ent_M)
            {
                MA = M.A[ja_M++];
                ent_M--;
            }
            if((J.ja[ja_J] == i) && ent_J)
            {
                JA = J.A[ja_J++];
                ent_J--;
            }

            if(JA != 0.0)
            {
                diff = (MA - JA) / std::abs(JA);
            }
            else
            {
                diff = (MA - JA);
            }

            if(std::abs(diff) > tol)
            {
                outFile << i << delimiter << j << delimiter << MA << delimiter
                        << JA << delimiter << diff << '\n';
                ndiff++;
            }
        }
    }

    outFile << "Total number of differences found: " << ndiff << '\n';
    // std::cout << "Jacobian::compare() -- INFO: Found " << ndiff
    //           << " difference(s).\n";

    outFile.close();

    return ndiff;
}

}  // namespace daecpp_namespace_name
