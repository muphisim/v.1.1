#!/bin/sh

rm -rf input output && mkdir input output
python perFlap.py
cp perFlap.inp $PWD/input
MuPhiSim -input perFlap.inp
