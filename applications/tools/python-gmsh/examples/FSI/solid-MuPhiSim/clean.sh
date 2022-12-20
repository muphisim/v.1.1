#!/bin/sh

echo 'remove  output precice-Solid*'
rm -rf output input perFlap.msh perFlap.inp

echo "---- Cleaning up preCICE logs in $(pwd)"
rm -fv ./precice-*-iterations.log \
    ./precice-*-convergence.log \
    ./precice-*-events.json \
    ./precice-*-events-summary.log \
    ./precice-postProcessingInfo.log \
    ./precice-*-watchpoint-*.log \
    ./precice-*-watchintegral-*.log \
    ./core
