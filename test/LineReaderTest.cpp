#include "csv.h"
#include "gtest/gtest.h"

TEST(LinReaderTest, test) {
    io::LineReader in("data/regression.csv");
    char *pch;
    while (char *line = in.next_line()) {
        pch = strtok(line, ",");
        while (pch != nullptr) {
            printf("%s\n", pch);
            pch = strtok(nullptr, ",");
        }
    }
}