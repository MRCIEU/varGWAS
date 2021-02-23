//
// Created by Matt Lyon on 18/02/2021.
//

#include "SynchronizedFile.h"
#include "Result.h"
#include <iostream>
#include <mutex>

namespace jlst {
void SynchronizedFile::write(const jlst::Result &res) {
  std::lock_guard<std::mutex> lock(writer_mutex);
  file << res.chromosome << "\t" << res.position << "\t" << res.rsid << "\t" << res.other_allele << "\t"
       << res.effect_allele << "\t" << res.beta << "\t" << res.se << "\t" << res.pval << "\t" << res.n << "\t"
       << res.eaf << "\n";
}
void SynchronizedFile::close() {
  file.close();
}
}