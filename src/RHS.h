/*
 * The RHS class.
 * This class is abstract and must be inherited.
 */

#pragma once

#include "typedefs.h"

namespace daecpp_namespace_name
{

class RHS
{
    std::size_t m_dump_file_counter = 0;

public:
    virtual ~RHS() {}

    /*
     * Takes vector x and time t and returns vector f.
     * This function is pure virtual and must be overriden.
     */
    virtual void operator()(const state_type &x, state_type &f,
                            const double t) = 0;

    /*
     * User-defined condition, when the solver should stop and return the
     * solution at the current time step
     */
    virtual bool stop_condition(const state_type &x, const double t)
    {
        return false;
    }

    /*
     * Helper function to write the RHS vector to a file
     */
    void dump(const state_type &x, const double t);
};

}  // namespace daecpp_namespace_name
