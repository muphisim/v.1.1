reset

for v in {0..3}
do 
    outDir=FEM/plate-$v 
    rm -rf $outDir
    mkdir -p $outDir
    python3 plate.py -m plate-$v.msh -outFolder $outDir -FEM 10 11
done


for v in {0..3}
do 
outDir=MM/plate-$v 
rm -rf $outDir
mkdir -p $outDir
python3 plate.py -m plate-$v.msh \
    -MMScalingFactor 2 \
    -outFolder $outDir \
    -MM 10 11 \
    -MMIntegOrder 4 \
    -MMAdaptiveSupportRadius 0
done

