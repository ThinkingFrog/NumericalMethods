cmake_minimum_required(VERSION 3.0)
project(SimpleIterations)

add_executable(iters src/main.c)
target_link_libraries(iters m)