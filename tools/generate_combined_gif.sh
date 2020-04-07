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


# parse command line options
while getopts ":f:o:" o; do
  case "${o}" in
    f)
      input_file=${OPTARG}
      ;;
    o)
      output_file=${OPTARG}
      ;;
    *)
      echo "Usage: $0 -f INPUT_FILE -o OUTPUT_FILE" 1>&2; exit 1
      ;;
  esac
done
shift $((OPTIND-1))

# check all required options were provided
if [ -z "${input_file}" ] || [ -z "${output_file}" ]; then
  echo "Usage: $0 -f INPUT_FILE -o OUTPUT_FILE" 1>&2; exit 1
fi

# generate a random 10 character string for the temporary folder name
rstr=$( head /dev/urandom | tr -dc A-Za-z | head -c 10 )
# create the temporary folder name
tfolder="temp${rstr}"

echo "Creating ${output_file} from ${input_file}"

echo "Generating 100 sphere frames..."
# call blender with the script and reroute its output to /dev/null
blender -b -P rotating_sphere.py -- \
  --file $input_file --output-folder $tfolder --output-prefix "frame" \
  > /dev/null 2> /dev/null
echo "Done."
echo "Generating 100 dust particle frames..."
# call blender with the script and reroute its output to /dev/null
blender -b -P rotating_dust_grain.py -- \
  --output-folder $tfolder --output-prefix "dust" \
  > /dev/null 2> /dev/null
echo "Done."
echo "Combining frames..."
for i in {0001..0100}
do
  montage $tfolder/frame$i.png $tfolder/dust$i.png -tile 1x2 -geometry +0+0 \
    -background none $tfolder/combined$i.png
done
echo "Done."
echo "Producing animated GIF..."
convert -delay 10 -dispose background $tfolder/combined*.png -loop 0 \
  $output_file
echo "Done."
echo "Cleaning up."
rm -rf $tfolder
