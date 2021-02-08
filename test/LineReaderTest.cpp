#include "csv.h"
#include "gtest/gtest.h"

TEST(LinReaderTest, test) {
    io::LineReader in("data/regression.csv");
    char *pch;
    while (char *line = in.next_line()) {
        pch = strtok(line, ",");
        printf("%s\n", pch);
        break;
        while (pch != nullptr) {
            pch = strtok(nullptr, ",");
            //printf("%s\n", pch);
            //break;
        }
    }
}