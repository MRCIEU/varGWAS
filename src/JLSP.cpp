//
// Created by Matthew Lyon on 19/03/2020.
//

#include "JLSP.h"
#include "libscl.h"
#include <iostream>

void jlst::JLSP::get_linear_estimate() {
    INTEGER n = 100;
    INTEGER p = 4;
    INT_32BIT seed = 100542;

    scl::realmat data(n, p);

    for (INTEGER i = 1; i <= n; i++) {
        data(i, 2) = 1.0;
        data(i, 3) = scl::ran(seed);
        data(i, 4) = scl::unsk(seed);
        data(i, 1) = data(i, 2) + data(i, 3) + data(i, 4) + 0.5 * scl::unsk(seed);
    }

    scl::realmat y = data("", 1);
    scl::realmat X = data("", scl::seq(2, data.ncol()));

    scl::realmat b = inv(T(X) * X) * T(X) * y;

    scl::realmat sse = T(y - X * b) * (y - X * b);

    scl::realmat V = sse[1] * inv(T(X) * X) / (y.nrow() - X.ncol());

    std::cout << scl::starbox("/Estimate of b and its variance//");

    std::cout << b << V;
}