%mem = 8GB
%nproc = 8
#b3lyp aug-cc-pvtz density=current nosym output=wfn 

21

0 1
O    0.0	-0.05	0.054316
H    0.0	-0.763369	-0.580667
H    0.0	0.763369	-0.580667

21_B.wfn

