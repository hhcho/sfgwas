import sys
import numpy as np

in_fname = sys.argv[1]
nrows = int(sys.argv[2])
ncols = int(sys.argv[3])
out_fname = sys.argv[4]

print("Called transposeMatrix.py:", in_fname, nrows, ncols, out_fname)

print('Loading matrix.. ', end='')
arr = np.fromfile(in_fname, dtype=np.int8)
print('done.')

print('Input array length:', len(arr))

arr = np.reshape(arr, (nrows, ncols))

print('Input dimensions:', arr.shape)

print('Transposing.. ', end='')
arr = arr.transpose()
print('done.')

print('Output dimensions:', arr.shape)

print('Writing to disk.. ', end='')
outfile = open(out_fname, 'wb')
arr.tofile(outfile)
outfile.close()
print('done.')
