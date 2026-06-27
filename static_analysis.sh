#!/bin/bash
set -e

echo "=== Python: Ruff lint ==="
ruff check skfmm/

echo "=== Python: Ruff format check ==="
ruff format --check skfmm/

echo "=== Python: mypy ==="
mypy skfmm/

echo "=== C++: cppcheck ==="
cppcheck --enable=all --std=c++17 --suppress=missingIncludeSystem skfmm/

echo "=== C++: clang-tidy ==="
clang-tidy skfmm/*.cpp -- -std=c++17 -I skfmm/

echo "=== C++: clang-format check ==="
clang-format --dry-run --Werror skfmm/*.cpp skfmm/*.h
