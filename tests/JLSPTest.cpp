//
// Created by Matthew Lyon on 19/03/2020.
//

#include "gtest/gtest.h"
#include "JLSP.h"


TEST(JLSPTestSuite, TestName) {
    jlst::JLSP t;
    t.get_abs_residuals_from_median();

    EXPECT_EQ(true, true);
}