#!/bin/bash

rootFolder=$1
outputFolder="convert"

if [ $(expr length "$rootFolder") -gt 0 ]; then
  outputFolder="$rootFolder/$outputFolder"
else
  rootFolder="./"
fi

echo "root folder: $rootFolder"
echo "output folder: $outputFolder"

# Create the output folder if it does not already exist
if [ ! -d "$outputFolder" ]; then
  mkdir "$outputFolder"
fi

# convert ppm files to jpg files
for i in "$rootFolder"/*.ppm; do
  name=$(basename "$i")
  name=${name::-4}
  echo "$name"
  convert "$i" "$outputFolder/$name.jpg"
done
