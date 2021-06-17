#!/bin/bash

outputLog="callgrind.log"
outputCallgrind="callgrind.profile"
outputStd="std.out"

SECONDS=0

#valgrind --tool=callgrind --log-file="${outputLog}" --callgrind-out-file="${outputCallgrind}" ./Debug/FoggySim
valgrind --tool=callgrind --collect-systime=yes --log-file="${outputLog}" --callgrind-out-file="${outputCallgrind}" ../Debug/FoggySim &>"${outputStd}"

echo $SECONDS

kcachegrind "${outputCallgrind}"
