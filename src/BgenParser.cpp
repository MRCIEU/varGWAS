
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <cassert>
#include <stdexcept>
#include "genfile/bgen/bgen.hpp"
#include "BgenParser.h"
#include "ProbSetter.h"

// BgenParser is a thin wrapper around the core functions in genfile/bgen/bgen.hpp.
// This class tracks file state and handles passing the right callbacks.

std::ostream &genfile::bgen::BgenParser::summarise(std::ostream &o) const {
    o << "BgenParser: bgen file ("
      << (m_context.flags & genfile::bgen::e_Layout2 ? "v1.2 layout" : "v1.1 layout")
      << ", "
      << (m_context.flags & genfile::bgen::e_CompressedSNPBlocks ? "compressed" : "uncompressed") << ")"
      << " with "
      << m_context.number_of_samples << " " << (m_have_sample_ids ? "named" : "anonymous") << " samples and "
      << m_context.number_of_variants << " variants.\n";
    return o;
};

// Attempt to read identifying information about a variant from the bgen file, returning
// it in the given fields.
// If this method returns true, data was successfully read, and it should be safe to call read_probs()
// or ignore_probs().
// If this method returns false, data was not successfully read indicating the end of the file.
bool genfile::bgen::BgenParser::read_variant(
        std::string *chromosome,
        uint32_t *position,
        std::string *rsid,
        std::vector<std::string> *alleles
) {
    assert(m_state == e_ReadyForVariant);
    std::string SNPID; // read but ignored in this toy implementation

    if (
            genfile::bgen::read_snp_identifying_data(
                    *m_stream, m_context,
                    &SNPID, rsid, chromosome, position,
                    [&alleles](std::size_t n) { alleles->resize(n); },
                    [&alleles](std::size_t i, std::string const &allele) { alleles->at(i) = allele; }
            )
            ) {
        m_state = e_ReadyForProbs;
        return true;
    } else {
        return false;
    }
};

std::size_t genfile::bgen::BgenParser::number_of_samples() const {
    return m_context.number_of_samples;
};

// Read genotype probability data for the SNP just read using read_variant()
// After calling this method it should be safe to call read_variant() to fetch
// the next variant from the file.
void genfile::bgen::BgenParser::read_probs(std::vector<std::vector<double> > *probs) {
    assert(m_state == e_ReadyForProbs);
    ProbSetter setter(probs);
    genfile::bgen::read_and_parse_genotype_data_block<ProbSetter>(
            *m_stream,
            m_context,
            setter,
            &m_buffer1,
            &m_buffer2
    );
    m_state = e_ReadyForVariant;
};

// Ignore genotype probability data for the SNP just read using read_variant()
// After calling this method it should be safe to call read_variant()
// to fetch the next variant from the file.
void genfile::bgen::BgenParser::ignore_probs() {
    genfile::bgen::ignore_genotype_data_block(*m_stream, m_context);
    m_state = e_ReadyForVariant;
};