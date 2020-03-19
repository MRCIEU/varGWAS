
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <memory>
#include "genfile/bgen/bgen.hpp"
#include "ProbSetter.h"

// Called once allowing us to set storage.
void genfile::bgen::ProbSetter::initialise(std::size_t number_of_samples, std::size_t number_of_alleles) {
    m_result->clear();
    m_result->resize(number_of_samples);
}

// If present with this signature, called once after initialise()
// to set the minimum and maximum ploidy and numbers of probabilities among samples in the data.
// This enables us to set up storage for the data ahead of time.
void genfile::bgen::ProbSetter::set_min_max_ploidy(uint32_t min_ploidy, uint32_t max_ploidy, uint32_t min_entries,
                                                   uint32_t max_entries) {
    for (std::size_t i = 0; i < m_result->size(); ++i) {
        m_result->at(i).reserve(max_entries);
    }
}

// Called once per sample to determine whether we want data for this sample
bool genfile::bgen::ProbSetter::set_sample(std::size_t i) {
    m_sample_i = i;
    // Yes, here we want info for all samples.
    return true;
}

// Called once per sample to set the number of probabilities that are present.
void genfile::bgen::ProbSetter::set_number_of_entries(std::size_t ploidy, std::size_t number_of_entries,
                                                      genfile::OrderType order_type, genfile::ValueType value_type) {
    assert(value_type == genfile::eProbability);
    m_result->at(m_sample_i).resize(number_of_entries);
    m_entry_i = 0;
}

// Called once for each genotype (or haplotype) probability per sample.
void genfile::bgen::ProbSetter::set_value(uint32_t, double value) {
    m_result->at(m_sample_i).at(m_entry_i++) = value;
}

// Ditto, but called if data is missing for this sample.
void genfile::bgen::ProbSetter::set_value(uint32_t, genfile::MissingValue value) {
    // Here we encode missing probabilities with -1
    m_result->at(m_sample_i).at(m_entry_i++) = -1;
}

// If present with this signature, called once after all data has been set.
void genfile::bgen::ProbSetter::finalise() {
    // nothing to do in this implementation.
}
