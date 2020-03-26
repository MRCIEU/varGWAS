//
// Created by Matthew Lyon on 19/03/2020.
//

#include <iostream>
#include <cassert>
#include <stdexcept>
#include "genfile/bgen/bgen.hpp"
#include "BgenParser.h"

int main(int argc, char **argv) {
    // TODO check file exists before passing to BgenParser
    std::string const filename = "/Users/matthewlyon/projects/jlst_cpp/lib/gavinband-bgen-44fcabbc5c38/example/example.v11.bgen";

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
            for (auto &prob : probs) {
                std::cout << '\t';

                // only support bi-allelic variants [0, 1, 2 copies of alt]
                assert(prob.size() == 3);

                // convert genotype probabilities to copies of alt
                std::cout << prob[1] + (2 * prob[2]);
            }

            std::cout << "\n";
        }

        return 0;
    } catch (genfile::bgen::BGenError const &e) {
        std::cerr << "!! Error parsing bgen file.\n";
        return -1;
    }

}