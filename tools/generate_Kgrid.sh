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

echo "Generating 100 frames..."
for column in {18..33}
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
  montage $tfolder/frame18$i.png $tfolder/frame19$i.png $tfolder/frame20$i.png \
    $tfolder/frame21$i.png $tfolder/frame22$i.png $tfolder/frame23$i.png \
    $tfolder/frame24$i.png $tfolder/frame25$i.png $tfolder/frame26$i.png \
    $tfolder/frame27$i.png $tfolder/frame28$i.png $tfolder/frame29$i.png \
    $tfolder/frame30$i.png $tfolder/frame31$i.png $tfolder/frame32$i.png \
    $tfolder/frame33$i.png -tile 4x4 -geometry +0+0 \
    -background none $tfolder/combined$i.png
done
echo "Done."
echo "Producing animated GIF..."
convert -delay 10 $tfolder/combined*.png -loop 0 $output_file
echo "Done."
echo "Frames are in ${tfolder}"
#echo "Cleaning up."
#rm -rf $tfolder
