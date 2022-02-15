cmake_minimum_required(VERSION 3.5)
# Project name, version & main languages
set(BIN adf)
project(${BIN} VERSION 0.0.1 LANGUAGES CXX)

include(GNUInstallDirs)

add_subdirectory(adf)
