# compute Orion
import sys
import numpy as np
import pandas as pd

n = int(sys.argv[1])
n_orion = 1662
N_eff = 10000 * n / n_orion
window = 501

chrom = sys.argv[2]

snvs_fn = sys.argv[3]
snvs = pd.read_csv(snvs_fn, sep='\t', header=None, names=['pos', 'ac'])
snvs_pos_ac = snvs.groupby('pos')['ac'].sum().reset_index().to_numpy()
snvs_pos_ac[:,1] = np.minimum(snvs_pos_ac[:,1], 2*n - snvs_pos_ac[:,1]) # fold if necessary

n_snvs = len(snvs_pos_ac)

E_i = np.array([(1/i) + (1/(2*n-i))/(1+(i==2*n-i)) for i in np.arange(1, n+1)])
weights = np.array([(2*n)/(i+1) for i in range(n+1)])

sfs = np.zeros(n+1)

nr_L = 0
nr_H = 0
for line in sys.stdin:
	L, H = [int(x) for x in line.rstrip().split('\t')[1:3]]
	mu = float(line.rstrip().split('\t')[-1])

	while snvs_pos_ac[nr_L, 0] <= L:
		sfs[snvs_pos_ac[nr_L, 1]] -= 1
		nr_L += 1

	while (nr_H<n_snvs) and (snvs_pos_ac[nr_H, 0] <= H):
		sfs[snvs_pos_ac[nr_H, 1]] += 1
		nr_H += 1

	theta = 4 * N_eff * mu
	psi_i = theta * E_i
	psi_0 = window - np.sum(psi_i)
	psi = np.concatenate((np.array([psi_0]), psi_i))
	
	sfs[0] = window - np.sum(sfs[1:])

	Orion = np.sum(weights * (psi - sfs)) / (100000 * theta)

	print(chrom, L+window//2, H-window//2, Orion, sep='\t')




