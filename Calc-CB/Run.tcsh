#!/bin/tcsh

set echo

python ../ErrMetrics-v2.py ka ../Exp/CBClip.txt Null1.txt Null2.txt BAR-ab-initio.txt TI-ab-initio.txt BAR-dock.txt TI-dock.txt TIxBAR.txt BEDAM.txt MovTyp-1.txt MovTyp-2.txt SOMD-1.txt SOMD-2.txt SOMD-3.txt SOMD-4.txt >& CBClip-Abs.dat &

python ../ErrMetrics-v2.py ka CorrectCB ../Exp/CBClip.txt Null1.txt Null2.txt BAR-ab-initio.txt TI-ab-initio.txt BAR-dock.txt TI-dock.txt TIxBAR.txt BEDAM.txt MovTyp-1.txt MovTyp-2.txt SOMD-1.txt SOMD-2.txt SOMD-3.txt SOMD-4.txt > & CBClip-RelCor.dat &

wait

