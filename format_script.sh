#! /bin/bash

################################################################################
 # This file is part of CosTuuM
 # Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 #
 # CosTuuM is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Affero General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option) any
 # later version.
 #
 # CosTuuM is distributed in the hope that it will be useful, but WITOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 # A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 # details.
 #
 # You should have received a copy of the GNU Affero General Public License
 # along with CosTuuM. If not, see <http://www.gnu.org/licenses/>.
 ###############################################################################

command -v clang-format-6.0 >/dev/null 2>&1 || \
  { echo >&2 "This script requires clang-format-6.0, but it is not installed!" \
             "Aborting."; exit 1; }

files=( src/*.cpp src/*.hpp src/*.cpp.in src/*.hpp.in test/*.cpp test/*.hpp \
        quicksched/*.[ch] )

echo "Formatting C++ files using clang-format-6.0..."
for f in "${files[@]}"
do clang-format-6.0 -style=file -i $f
done
echo "Done."

echo "Formatting Python files using black..."
python3 -m black -l 80 {test,tools,benchmark,scripts}/*.py
python3 -m black -l 80 setup.py.in
echo "Done."
