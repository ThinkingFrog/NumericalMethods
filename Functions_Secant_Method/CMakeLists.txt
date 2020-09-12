cmake_minimum_required(VERSION 3.0)
project(SecantMethod)
set(LANGUAGE C)

add_executable(Secant src/main.c)
target_link_libraries(Secant m)