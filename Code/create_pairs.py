import numpy as np
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Provide an input phenotype file (consisting only of traits as columns) to generate an xy_array of indices encoding all possible trait pairs.")

parser.add_argument('-p', '--pheno', type = str, help = 'Specify the phenotype file.')

args = parser.parse_args()
pheno = args.pheno
pheno = pd.read_csv(pheno)

p = pheno.shape[1]
xy_array = pd.DataFrame(np.tril_indices(p,k=-1)).T + 1
np.savetxt('xy_array.txt', xy_array, fmt='%i',delimiter=',')
print("Done!") 