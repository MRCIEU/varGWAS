
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


namespace genfile {
    namespace bgen {

        struct ProbSetter {
            typedef std::vector<std::vector<double> > Data;

            explicit ProbSetter(Data *result);

            void initialise(std::size_t number_of_samples, std::size_t number_of_alleles);

            void
            set_min_max_ploidy(uint32_t min_ploidy, uint32_t max_ploidy, uint32_t min_entries, uint32_t max_entries);

            bool set_sample(std::size_t i);

            void set_number_of_entries(std::size_t ploidy, std::size_t number_of_entries, genfile::OrderType order_type,
                                       genfile::ValueType value_type);

            void set_value(uint32_t, double value);

            void set_value(uint32_t, genfile::MissingValue value);

            static void finalise();

        private:
            Data *m_result;
            std::size_t m_sample_i;
            std::size_t m_entry_i{};
        };

    };
}

#endif //JLST_CPP_PROBSETTER_H