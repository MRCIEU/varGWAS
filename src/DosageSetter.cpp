
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <cassert>
#include "genfile/bgen/bgen.hpp"

#include "DosageSetter.h"

// DosageSetter is a callback object appropriate for passing to bgen::read_genotype_data_block() or
// the synonymous method of genfile::bgen::View.
// Its job is to compute genotype dosage (dosage of the 2nd allele).
// This implementation assumes samples are diploid, data unphased, and variants are biallelic
// and will assert() if not.

// Called once allowing us to set storage.
void genfile::bgen::DosageSetter::initialise(std::size_t number_of_samples, std::size_t number_of_alleles) {
    m_result = 0;
    m_n = 0;
}

// If present with this signature, called once after initialise()
// to set the minimum and maximum ploidy and numbers of probabilities among samples in the data.
// This enables us to set up storage for the data ahead of time.
void genfile::bgen::DosageSetter::set_min_max_ploidy(uint32_t min_ploidy, uint32_t max_ploidy, uint32_t min_entries,
                                                     uint32_t max_entries) {
    assert(min_ploidy == 2);
    assert(max_ploidy == 2);
}

// Called once per sample to determine whether we want data for this sample
bool genfile::bgen::DosageSetter::set_sample(std::size_t i) {
    m_sample_i = i;
    // Yes, here we want info for all samples.
    return true;
}

// Called once per sample to set the number of probabilities that are present.
void genfile::bgen::DosageSetter::set_number_of_entries(
        std::size_t ploidy,
        std::size_t number_of_entries,
        genfile::OrderType order_type,
        genfile::ValueType value_type
) {
    assert(value_type == genfile::eProbability);
    assert(order_type == genfile::ePerUnorderedGenotype);
    assert(ploidy == 2);
    assert(number_of_entries == 3);

    ++m_n;
}

// Called once for each genotype (or haplotype) probability per sample.
void genfile::bgen::DosageSetter::set_value(uint32_t entry_i, double value) {
    m_result += entry_i * value;
}

// Ditto, but called if data is missing for this sample.
void genfile::bgen::DosageSetter::set_value(uint32_t, genfile::MissingValue value) {
}

// If present with this signature, called once after all data has been set.
void genfile::bgen::DosageSetter::finalise() {
    // nothing to do in this implementation.
    m_result /= m_n;
}

double const &genfile::bgen::DosageSetter::result() const { return m_result; }