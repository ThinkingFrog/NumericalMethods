cmake_minimum_required(VERSION 3.0)
project(LUFactorization)

add_executable(LUFact src/main.c)
target_link_libraries(LUFact m)