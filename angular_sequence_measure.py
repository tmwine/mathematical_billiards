'''

given a file of symbolic sequences, assigns a value to each according to
their symbol proportions (order does not matter); the more disproportionate
a sequence's symbols are, the lower the value, in range (0,pi/2]

this requires specification of 2 files in header area:
	in_file--a file of symbolic sequences, one on each line
	out_file--a file to save the angular proportion measure of each
bounce sequence in, one entry per line

'''

import numpy as np

#@@@ set I/O filenames:
in_file = ""
out_file = ""
#@@@

if in_file=="" or out_file=="":
    print("enter names of both input and output files in script header")
    raise ValueError


sym_vec = ['1','2','3']
ct_vc = np.array([0.0,0.0,0.0])
ref_vec = np.array([1/np.sqrt(3) for x in range(3)])
val_vec = []
with open(in_file) as fp:
    for line in fp:
        for ii,sym in enumerate(sym_vec):
            ct_vc[ii] = line.count(sym)
        val = np.arccos(np.dot(ct_vc/np.linalg.norm(ct_vc),ref_vec))
        val_vec.append(val)

with open(out_file,'w') as fp:
    for val in val_vec:
        tmp = fp.write(str(np.pi/2-val)+'\n')


