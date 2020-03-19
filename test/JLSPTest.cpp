#include "gtest/gtest.h"
#include "JLSP.h"

TEST(blaTest, test1) {
    jlst::JLSP t;
    t.get_linear_estimate();

    EXPECT_EQ (true, true);
}