#! /bin/bash

# parse command line options
while getopts ":c:f:o:" o; do
  case "${o}" in
    c)
      column=${OPTARG}
      ;;
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

# set default values for optional arguments that were not provided
if [ -z "${column}" ]; then
  column=2
  echo "Using default column ${column}"
fi

# generate a random 10 character string for the temporary folder name
rstr=$( head /dev/urandom | tr -dc A-Za-z | head -c 10 )
# create the temporary folder name
tfolder="temp${rstr}"

echo "Creating ${output_file} from ${input_file}"

echo "Generating 100 frames..."
# call blender with the script and reroute its output to /dev/null
blender -b -P rotating_sphere.py -- \
  --file $input_file --output-folder $tfolder --output-prefix "frame" \
  --column $column > /dev/null 2> /dev/null
echo "Done."
echo "Cropping frames..."
for f in $tfolder/frame*.png
do
  mogrify -trim +repage $f
done
echo "Done."
echo "Producing animated GIF..."
convert -delay 10 $tfolder/frame*.png -loop 0 $output_file
echo "Done."
echo "Cleaning up."
rm -rf $tfolder
