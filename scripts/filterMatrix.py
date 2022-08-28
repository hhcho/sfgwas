import sys
import numpy as np

in_fname = sys.argv[1]
nrows = int(sys.argv[2])
ncols = int(sys.argv[3])
row_filt_fname = sys.argv[4]
col_filt_fname = sys.argv[5]
out_fname = sys.argv[6]

print("Called filterMatrix.py:", in_fname, nrows, ncols, row_filt_fname, col_filt_fname, out_fname)

print('Loading filters.. ', end='')
row_filt = np.fromfile(row_filt_fname, dtype=np.uint8)
col_filt = np.fromfile(col_filt_fname, dtype=np.uint8)
print('done.')

nrows_out = row_filt.sum()
ncols_out = col_filt.sum()

assert nrows == len(row_filt)
assert ncols == len(col_filt)

print("Row filter: %d/%d" % (nrows_out, nrows))
print("Column filter: %d/%d" % (ncols_out, ncols))

infile = open(in_fname, 'rb')
outfile = open(out_fname, 'wb')

print('Streaming matrix.. ', end='')
for r in range(nrows):
    arr = np.fromfile(infile, count=len(col_filt), dtype=np.int8)
    if row_filt[r]:
        arr = np.extract(col_filt, arr)
        arr.tofile(outfile)
print('done.')

infile.close()
outfile.close()

## For debugging
#arr = np.fromfile(in_fname, dtype=np.uint8)
#arr = np.reshape(arr, (nrows, ncols))
#print(arr)
#
#arr = np.fromfile(out_fname, dtype=np.uint8)
#arr = np.reshape(arr, (nrows_out, ncols_out))
#print(arr)
