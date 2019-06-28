#!/bin/bash
agp="$1"
scaf="$2"
nam=$(basename $scaf |cut -f 1 -d ".")
python -m jcvi.formats.agp build $agp $scaf ${nam}.pseudomolecules-v2.fasta
