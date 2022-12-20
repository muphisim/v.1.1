rm -rf input output
mkdir input output
python3 model.py
cp input.inp input
MuPhiSim -input input.inp
