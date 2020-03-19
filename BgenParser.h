
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef JLST_CPP_BGENPARSER_H
#define JLST_CPP_BGENPARSER_H

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <memory>
#include "genfile/bgen/bgen.hpp"

// BgenParser is a thin wrapper around the core functions in genfile/bgen/bgen.hpp.
// This class tracks file state and handles passing the right callbacks.
namespace genfile {
    namespace bgen {

        struct BgenParser {
            BgenParser(std::string const &filename) :
                    m_filename(filename),
                    m_state(e_NotOpen),
                    m_have_sample_ids(false) {
                // Open the stream
                m_stream.reset(
                        new std::ifstream(filename, std::ifstream::binary)
                );
                if (!*m_stream) {
                    throw std::invalid_argument(filename);
                }
                m_state = e_Open;

                // Read the offset, header, and sample IDs if present.
                genfile::bgen::read_offset(*m_stream, &m_offset);
                genfile::bgen::read_header_block(*m_stream, &m_context);
                if (m_context.flags & genfile::bgen::e_SampleIdentifiers) {
                    genfile::bgen::read_sample_identifier_block(
                            *m_stream, m_context,
                            [this](std::string id) { m_sample_ids.push_back(id); }
                    );
                    m_have_sample_ids = true;
                }

                // Jump to the first variant data block.
                m_stream->seekg(m_offset + 4);

                // We keep track of state (though it's not really needed for this implementation.)
                m_state = e_ReadyForVariant;
            }


            std::ostream &summarise(std::ostream &o) const;

            std::size_t number_of_samples() const;

            // Report the sample IDs in the file using the given setter object
            // (If there are no sample IDs in the file, we report a dummy identifier).
            template<typename Setter>
            void get_sample_ids(Setter setter) {
                if (m_have_sample_ids) {
                    for (std::size_t i = 0; i < m_context.number_of_samples; ++i) {
                        setter(m_sample_ids[i]);
                    }
                } else {
                    for (std::size_t i = 0; i < m_context.number_of_samples; ++i) {
                        setter("S" + std::to_string(i + 1) + "");
                    }
                }
            };

            // Attempt to read identifying information about a variant from the bgen file, returning
            // it in the given fields.
            // If this method returns true, data was successfully read, and it should be safe to call read_probs()
            // or ignore_probs().
            // If this method returns false, data was not successfully read indicating the end of the file.
            bool read_variant(
                    std::string *chromosome,
                    uint32_t *position,
                    std::string *rsid,
                    std::vector<std::string> *alleles
            );

            // Read genotype probability data for the SNP just read using read_variant()
            // After calling this method it should be safe to call read_variant() to fetch
            // the next variant from the file.
            void read_probs(std::vector<std::vector<double> > *probs);

            // Ignore genotype probability data for the SNP just read using read_variant()
            // After calling this method it should be safe to call read_variant()
            // to fetch the next variant from the file.
            void ignore_probs();

        private:
            std::string const m_filename;
            std::unique_ptr<std::istream> m_stream;

            // bgen::Context object holds information from the header block,
            // including bgen flags
            genfile::bgen::Context m_context;

            // offset byte from top of bgen file.
            uint32_t m_offset;

            // We keep track of our state in the file.
            // Not strictly necessary for this implentation but makes it clear that
            // calls must be read_variant() followed by read_probs() (or ignore_probs())
            // repeatedly.
            enum State {
                e_NotOpen = 0, e_Open = 1, e_ReadyForVariant = 2, e_ReadyForProbs = 3, eComplete = 4
            };
            State m_state;

            // If the BGEN file contains samples ids, they will be read here.
            bool m_have_sample_ids;
            std::vector<std::string> m_sample_ids;

            // Buffers, these are used as working space by bgen implementation.
            std::vector<genfile::byte_t> m_buffer1, m_buffer2;
        };

    }
}

#endif