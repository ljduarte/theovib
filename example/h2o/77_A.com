%mem = 8GB
%nproc = 8
#b3lyp aug-cc-pvtz density=current nosym output=wfn 

77

0 1
O    0.0	0.0	0.004316
H    0.0	-0.763369	-0.580667
H    0.0	0.813369	-0.580667

77_A.wfn

