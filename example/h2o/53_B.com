%mem = 8GB
%nproc = 8
#b3lyp aug-cc-pvtz density=current nosym output=wfn 

53

0 1
O    0.0	0.0	0.004316
H    -0.05	-0.763369	-0.530667
H    0.0	0.763369	-0.580667

53_B.wfn

