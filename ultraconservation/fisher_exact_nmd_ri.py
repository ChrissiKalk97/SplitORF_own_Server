# NMD
import numpy as np
from scipy.stats import fisher_exact
table = table = np.array([[7, 271], [3, 1708]])
res = fisher_exact(table, alternative='greater')
res.pvalue
# np.float64(7.995869822999143e-05)


table = table = np.array([[1, 523], [5, 1551]])
res = fisher_exact(table, alternative='greater')
res.pvalue
