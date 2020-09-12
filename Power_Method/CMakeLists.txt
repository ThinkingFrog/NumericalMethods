cmake_minimum_required(VERSION 3.0)
project(PowerMethod)

add_executable(Power src/main.c)
target_link_libraries(Power m)