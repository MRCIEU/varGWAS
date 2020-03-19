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

int main(int argc, char **argv) {
    std::ios_base::sync_with_stdio(false);
    std::string const filename = "../includes/bgen/example/example.v11.bgen";
    std::string mode = "default";

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

            std::cout << rsid << ": " << dc.result() << "\n";
        }

    } catch (std::invalid_argument const &e) {
        std::cerr << "!! Error: " << e.what() << ".\n";
        return -1;
    } catch (genfile::bgen::BGenError const &e) {
        std::cerr << "!! Error parsing bgen file.\n";
        return -1;
    }
}