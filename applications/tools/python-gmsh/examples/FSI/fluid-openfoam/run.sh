#!/bin/sh
set -e

python createMesh.py
touch fluid-openfoam.foam

# OpenFOAM run functions: getApplication, getNumberOfProcessors
# shellcheck disable=SC1090 # This is an OpenFOAM file which we don't need to check
. "${WM_PROJECT_DIR}/bin/tools/RunFunctions"

checkMesh
renumberMesh -overwrite

solver=$(getApplication)
if [ "${1:-}" = "-parallel" ]; then
    procs=$(getNumberOfProcessors)
    decomposePar -force
    mpirun -np "${procs}" "${solver}" -parallel
    reconstructPar
else
    ${solver}
fi

