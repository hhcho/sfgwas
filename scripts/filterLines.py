import sys
import numpy as np

in_fname = sys.argv[1]
bin_filt_fname = sys.argv[2]
shift = int(sys.argv[3])
out_fname = sys.argv[4]

print("Called filterLines.py:", in_fname, bin_filt_fname, shift, out_fname)

print('Loading filters.. ', end='')
filt = np.fromfile(bin_filt_fname, dtype=np.uint8)
print('done.')

print("Filter: %d/%d" % (np.sum(filt), len(filt)))
print("Starting index: %d" % shift)

infile = open(in_fname, 'rt')
outfile = open(out_fname, 'wt')

print('Streaming file.. ', end='')
lineno = -1
for line in infile:
    lineno += 1
    if lineno < shift:
        continue
    if lineno >= shift + len(filt):
        break
    if filt[lineno - shift]:
        outfile.write(line)    
print('done.')    

infile.close()
outfile.close()