import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
import numpy as np
from scipy.stats import beta


fig, axs = plt.subplots(1, 1)

axs.set_xlabel("(a,b)")
axs.set_ylabel("AUC under binomial sampling")

aucs = []
N = int(1e5)
param_tuples = [(0.5,0.5),(1,1),(1000,1),(1,1000)]
colors = ["red", "green", "blue", "orange"]

for params in param_tuples:
    a = params[0]
    b = params[1]

    prob_samples = np.random.beta(a, b, N)
    y_samples = np.random.binomial(p=prob_samples, n=1)
    
    auc = roc_auc_score(y_samples, prob_samples)
    aucs.append(auc)
    
axs[1].bar(list(map(str, param_tuples)), aucs, color=colors)
fig.savefig("/home/jwt44/auc_by_beta_params.pdf")



fig, ax = plt.subplots(1, 1)
ax.set_ylim(0,3)

color_by_param = {(0.5,0.5): 'red', (1,1): 'green', (2,5): 'blue', (5,5): 'orange'}

for params, color in color_by_param.items():
    a, b = params
    x = np.linspace(0, 1, 100)
    ax.plot(x, beta.pdf(x, a, b),
           'r-', lw=5, alpha=1, label='beta pdf', color=color)

fig.savefig("/home/jwt44/beta_densities.pdf")



