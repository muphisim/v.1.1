#!/bin/sh
set -e -u

set -e -u
echo "--- Cleaning up OpenFOAM case in $(pwd)"
if [ -n "${WM_PROJECT:-}" ] || error "No OpenFOAM environment is active."; then
    # shellcheck disable=SC1090 # This is an OpenFOAM file which we don't need to check
    . "${WM_PROJECT_DIR}/bin/tools/CleanFunctions"
    cleanCase
    rm -rfv 0/uniform/functionObjects/functionObjectProperties
fi
rm -rfv ./preCICE-output/

echo "---- Cleaning up preCICE logs in $(pwd)"
rm -fv ./precice-*-iterations.log \
    ./precice-*-convergence.log \
    ./precice-*-events.json \
    ./precice-*-events-summary.log \
    ./precice-postProcessingInfo.log \
    ./precice-*-watchpoint-*.log \
    ./precice-*-watchintegral-*.log \
    ./core
