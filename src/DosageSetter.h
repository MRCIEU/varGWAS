//
// Created by Matthew Lyon on 19/03/2020.
//

#ifndef JLST_CPP_DOSAGESETTER_H
#define JLST_CPP_DOSAGESETTER_H

namespace genfile {
    namespace bgen {
        struct DosageSetter {
            DosageSetter() :
                    m_sample_i(0) {}

            void initialise(std::size_t number_of_samples, std::size_t number_of_alleles);

            void
            set_min_max_ploidy(uint32_t min_ploidy, uint32_t max_ploidy, uint32_t min_entries, uint32_t max_entries);

            bool set_sample(std::size_t i);

            void set_number_of_entries(std::size_t ploidy, std::size_t number_of_entries, genfile::OrderType order_type,
                                       genfile::ValueType value_type);

            void set_value(uint32_t entry_i, double value);

            void set_value(uint32_t, genfile::MissingValue value);

            void finalise();

            double const &result() const;

        private:
            double m_result;
            std::size_t m_n;
            std::size_t m_sample_i;
            std::size_t m_entry_i;
            typedef std::vector<double> Data;
        };
    }
}

#endif //JLST_CPP_DOSAGESETTER_H
