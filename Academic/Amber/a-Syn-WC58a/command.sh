for i in {2..15}; do
    cd State$i
    cp ../common/* .
    python process_complex_for_amber.py
    tleap -s -f tleap_complex.in > tleap_complex.log
    python write_decompose_residues.py
    cd ..
done
 
