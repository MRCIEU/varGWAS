
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef JLST_CPP_PROBSETTER_H
#define JLST_CPP_PROBSETTER_H

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <memory>
#include "genfile/bgen/bgen.hpp"

// ProbSetter is a callback object appropriate for passing to bgen::read_genotype_data_block() or
// the synonymous method of genfile::bgen::View. See the comment in bgen.hpp above
// bgen::read_genotype_data_block(), or the bgen wiki for a description of the API.
// The purpose of this object is to store genotype probability values in the desired
// data structure (which here is a vector of vectors of doubles).
namespace genfile {
    namespace bgen {
        struct ProbSetter {
            typedef std::vector<std::vector<double> > Data;

            ProbSetter(Data *result) :
                    m_result(result),
                    m_sample_i(0) {}

            void initialise(std::size_t number_of_samples, std::size_t number_of_alleles);

            void
            set_min_max_ploidy(uint32_t min_ploidy, uint32_t max_ploidy, uint32_t min_entries, uint32_t max_entries);

            bool set_sample(std::size_t i);

            void set_number_of_entries(std::size_t ploidy, std::size_t number_of_entries, genfile::OrderType order_type,
                                       genfile::ValueType value_type);

            void set_value(uint32_t, double value);

            void set_value(uint32_t, genfile::MissingValue value);

            void finalise();

        private:
            Data *m_result;
            std::size_t m_sample_i;
            std::size_t m_entry_i;
        };
    };
}

#endif