%mem = 8GB
%nproc = 8
#b3lyp aug-cc-pvtz density=current nosym output=wfn 

75

0 1
O    0.0	0.0	0.004316
H    0.0	-0.763369	-0.6306670000000001
H    0.0	0.7133689999999999	-0.580667

75_D.wfn

