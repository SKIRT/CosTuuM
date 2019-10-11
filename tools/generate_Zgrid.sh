#! /bin/bash

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

echo "Generating 1 frame..."
for column in {2..17}
do
# call blender with the script and reroute its output to /dev/null
blender -b -P rotating_sphere.py -- \
  --file $input_file --output-folder $tfolder --output-prefix "frame${column}" \
  --column $column -v 10 -s 1 > /dev/null 2> /dev/null
done
echo "Done."
echo "Cropping frames..."
for f in $tfolder/frame*.png
do
  mogrify -trim +repage $f
done
echo "Done."
echo "Collating frames..."
for i in {0001..0100}
do
  montage $tfolder/frame2$i.png $tfolder/frame3$i.png $tfolder/frame4$i.png \
    $tfolder/frame5$i.png $tfolder/frame6$i.png $tfolder/frame7$i.png \
    $tfolder/frame8$i.png $tfolder/frame9$i.png $tfolder/frame10$i.png \
    $tfolder/frame11$i.png $tfolder/frame12$i.png $tfolder/frame13$i.png \
    $tfolder/frame14$i.png $tfolder/frame15$i.png $tfolder/frame16$i.png \
    $tfolder/frame17$i.png -tile 4x4 -geometry +0+0 \
    -background none $tfolder/combined$i.png
done
echo "Done."
echo "Producing animated GIF..."
convert -delay 10 $tfolder/combined*.png -loop 0 $output_file
echo "Done."
echo "Frames are in ${tfolder}"
#echo "Cleaning up."
#rm -rf $tfolder
