#!/bin/tcsh

set echo

python ../ErrMetrics-v2.py ka ../Exp/OAAllAvg.txt Null1.txt Null2.txt APR-OPC.txt APR-TIP3P.txt BEDAM.txt Metadynamics.txt MMPBSA-GAFF.txt MovTyp-1.txt MovTyp-2.txt PERT-bound-c.txt PERT-bound.txt PERT-hrex-c1.txt PERT-hrex-c2.txt PERT-hrex-c.txt PERT-hrex.txt SOMD-1.txt SOMD-2.txt SOMD-3.txt SOMD-4.txt > & AllAvg-Abs.dat &

python ../ErrMetrics-v2.py ka ../Exp/OASaltDep.txt Null1.txt Null2.txt APR-OPC.txt APR-TIP3P.txt BEDAM.txt Metadynamics.txt MMPBSA-GAFF.txt MovTyp-1.txt MovTyp-2.txt PERT-bound-c.txt PERT-bound.txt PERT-hrex-c1.txt PERT-hrex-c2.txt PERT-hrex-c.txt PERT-hrex.txt SOMD-1.txt SOMD-2.txt SOMD-3.txt SOMD-4.txt > & SaltDep-Abs.dat &

python ../ErrMetrics-v2.py ka CorrectOA ../Exp/OAAllAvg.txt Null1.txt Null2.txt APR-OPC.txt APR-TIP3P.txt BEDAM.txt CCSD\(T\)-neutral.txt DFT-charged.txt DFT-neutral.txt Metadynamics.txt MMPBSA-GAFF.txt MovTyp-1.txt MovTyp-2.txt PERT-bound-c.txt PERT-bound.txt PERT-hrex-c1.txt PERT-hrex-c2.txt PERT-hrex-c.txt PERT-hrex.txt SOMD-1.txt SOMD-2.txt SOMD-3.txt SOMD-4.txt > & AllAvg-RelCor.dat &

python ../ErrMetrics-v2.py ka CorrectOA ../Exp/OASaltDep.txt Null1.txt Null2.txt APR-OPC.txt APR-TIP3P.txt BEDAM.txt CCSD\(T\)-neutral.txt DFT-charged.txt DFT-neutral.txt Metadynamics.txt MMPBSA-GAFF.txt MovTyp-1.txt MovTyp-2.txt PERT-bound-c.txt PERT-bound.txt PERT-hrex-c1.txt PERT-hrex-c2.txt PERT-hrex-c.txt PERT-hrex.txt SOMD-1.txt SOMD-2.txt SOMD-3.txt SOMD-4.txt > & SaltDep-RelCor.dat &

python ../ErrMetrics-v2.py ka CorrectOA OAHOnly ../Exp/OASaltDep.txt Null1.txt Null2.txt APR-OPC.txt APR-TIP3P.txt BEDAM.txt CCSD\(T\)-neutral.txt DFT-charged.txt DFT-neutral.txt Metadynamics.txt MMPBSA-GAFF.txt MovTyp-1.txt MovTyp-2.txt PERT-bound-c.txt PERT-bound.txt PERT-hrex-c1.txt PERT-hrex-c2.txt PERT-hrex-c.txt PERT-hrex.txt SOMD-1.txt SOMD-2.txt SOMD-3.txt SOMD-4.txt PERT-combo.txt > & SaltDep-RelCor-OAH.dat &

python ../ErrMetrics-v2.py ka CorrectOA OAMeOnly ../Exp/OASaltDep.txt Null1.txt Null2.txt APR-OPC.txt APR-TIP3P.txt BEDAM.txt CCSD\(T\)-neutral.txt DFT-charged.txt DFT-neutral.txt Metadynamics.txt MMPBSA-GAFF.txt MovTyp-1.txt MovTyp-2.txt PERT-bound-c.txt PERT-bound.txt PERT-hrex-c1.txt PERT-hrex-c2.txt PERT-hrex-c.txt PERT-hrex.txt SOMD-1.txt SOMD-2.txt SOMD-3.txt SOMD-4.txt MMPBSA-OPLS.txt > & SaltDep-RelCor-OAMe.dat &

python ../ErrMetrics-v2.py dh ../Exp/OAEnth10.txt Enth-Null1.txt Enth-Null2.txt Enth-APR_OPC-10.txt Enth-APR_TIP3P-10.txt > & Enthalpy.dat &

wait

