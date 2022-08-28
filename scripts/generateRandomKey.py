import numpy as np
import sys

# Warning: only meant for toy example generation purposes
# Do not use these keys in real settings

outfile = sys.argv[1] # Output file
numBytes = int(sys.argv[2]) # Number of random bytes

b = np.random.bytes(numBytes)

fp = open(outfile, 'wb')
fp.write(b)
fp.close()