//
// Created by Matthew Lyon on 19/03/2020.
//

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <memory>
#include "genfile/bgen/bgen.hpp"
#include "genfile/bgen/View.hpp"
#include "DosageSetter.h"
#include "JLSP.h"
#include "BgenParser.h"

int main(int argc, char **argv) {
    std::ios_base::sync_with_stdio(false);
    std::string const filename = "../includes/bgen/example/example.v11.bgen";

    // vcf example
    try {
        genfile::bgen::BgenParser bgenParser(filename);

        bgenParser.summarise(std::cerr);

        // Output header
        std::cout << "CHROM\tPOS\tID\tREF\tALT";
        bgenParser.get_sample_ids(
                [](std::string const &id) { std::cout << "\t" << id; }
        );
        std::cout << "\n";


        // Output variants
        std::string chromosome;
        uint32_t position;
        std::string rsid;
        std::vector<std::string> alleles;
        std::vector<std::vector<double> > probs;

        while (bgenParser.read_variant(&chromosome, &position, &rsid, &alleles)) {
            // only support bi-allelic variants
            assert(alleles.size() == 2);

            // print variant
            std::cout << chromosome << '\t'
                      << position << '\t'
                      << rsid << '\t'
                      << alleles[0] << '\t'
                      << alleles[1] << '\t';

            bgenParser.read_probs(&probs);
            for (std::size_t i = 0; i < probs.size(); ++i) {
                std::cout << '\t';

                // only support bi-allelic variants [0, 1, 2 copies of alt]
                assert(probs[i].size() == 3);

                // convert genotype probabilities to copies of alt
                std::cout << probs[i][1] + (2 * probs[i][2]);

            }

            std::cout << "\n";
        }
        return 0;
    }
    catch (genfile::bgen::BGenError const &e) {
        std::cerr << "!! Error parsing bgen file.\n";
        return -1;
    }


    return 0;

    // dosage example
    try {
        std::string SNPID, rsid, chromosome;
        uint32_t position;
        std::vector<std::string> alleles;
        genfile::bgen::View view(filename);

        for (std::size_t i = 0; i < view.number_of_variants(); ++i) {
            bool success = view.read_variant(
                    &SNPID, &rsid, &chromosome, &position, &alleles
            );
            assert(success);

            genfile::bgen::DosageSetter dc;
            view.read_genotype_data_block(dc);
            std::cout << chromosome + "-" + std::to_string(position) + "-" + alleles[0] + "-" + alleles[1] << ":"
                      << dc.result() << "\n";

        }

    } catch (std::invalid_argument const &e) {
        std::cerr << "!! Error: " << e.what() << ".\n";
        return -1;
    } catch (genfile::bgen::BGenError const &e) {
        std::cerr << "!! Error parsing bgen file.\n";
        return -1;
    }
}