#! /bin/bash

command -v clang-format-6.0 >/dev/null 2>&1 || \
  { echo >&2 "This script requires clang-format-6.0, but it is not installed!" \
             "Aborting."; exit 1; }

files=( src/*.cpp src/*.hpp src/*.cpp.in src/*.hpp.in test/*.cpp test/*.hpp )

echo "Formatting C++ files using clang-format-6.0..."
for f in "${files[@]}"
do clang-format-6.0 -style=file -i $f
done
echo "Done."

echo "Formatting Python files using black..."
python3 -m black -l 80 test/*.py
echo "Done."
