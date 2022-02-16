#ifndef POLY_TESTS_H
#define POLY_TESTS_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <cmath>
#include <gtest/gtest.h>

#include "inc/polynomial.hpp"

TEST(PolySuite, AddOperator)
{
    adf::polynomial<int> poly1{1, 2, 3, 4, 5, 7, 9, 3, 5, 8};
    adf::polynomial<int> poly2{1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    adf::polynomial<int> result = poly1 + poly2;

    for(std::size_t i = 0; i < poly1.size(); ++i)
    {
        ASSERT_EQ(result[i], (poly1[i] + poly2[i]));
    }

    adf::polynomial<int> poly3{5, 7, 3, 5, 2, 1, 1, 7, 3, 4};
    adf::polynomial<int> poly4{1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    poly4 += poly3;

    for(std::size_t i = 0; i < poly1.size(); ++i)
    {
        ASSERT_EQ(poly4[i], (poly3[i] + 1));
    }

}

TEST(PolySuite, SubOperator)
{

}

TEST(PolySuite, Multoperator)
{

}

TEST(PolySuite, DivOperator)
{

}

TEST(PolySuite, RemOperator)
{

}

TEST(PolySuite, Equal)
{

}

TEST(PolySuite, NotEqual)
{

}

TEST(PolySuite, Index)
{

}

#endif //POLY_TESTS_H
