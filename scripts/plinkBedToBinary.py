import sys
import numpy as np
import math

in_fname = sys.argv[1]
num_sample = int(sys.argv[2])
num_snp = int(sys.argv[3])
out_fname = sys.argv[4]

print("Called plinkBedToBinary.py:", in_fname, num_sample, num_snp, out_fname)

x = np.fromfile(in_fname, dtype=np.uint8)[3:] # Skip magic numbers

byte_per_snp = int(math.ceil(num_sample / 4.0))

assert(len(x) == num_snp * byte_per_snp)

masks = [3, 12, 48, 192]
y = np.zeros((4, len(x)), dtype=np.int8)
for i in range(len(masks)):
    z = np.right_shift(np.bitwise_and(x, masks[i]), 2 * i).astype(np.int8)
    z0 = z == 0
    z1 = z == 1
    z[z0] = 1
    z = 3 - z
    z[z1] = -1
    y[i] = z
y = y.transpose().reshape((num_snp,-1)).transpose()[:num_sample]

print('Exporting matrix.. ', end='')
outfile = open(out_fname, 'wb')
y.tofile(outfile)
outfile.close()
print('done.')