import cvxpy as cp
import numpy as np
import matplotlib.pyplot as plt

# Constants - TEST
# mu1 = 8
# mu2 = 20
# sigma1 = 6
# sigma2 = 17.5
# rho = -0.25
# n = 100
# r_min = -30
# r_max = 70

# Constants - NEW
mu1 = 5
mu2 = 15
sigma1 = 3
sigma2 = 10
rho = -0.3
n = 100
r_min = -20
r_max = 50

# Discretized outcomes for R1 and R2
r = np.linspace(r_min, r_max, n)

# Marginal distributions for R1 and R2
p1 = np.exp(-(r - mu1) ** 2 / (2 * sigma1 ** 2))
p1 /= np.sum(p1)

p2 = np.exp(-(r - mu2) ** 2 / (2 * sigma2 ** 2))
p2 /= np.sum(p2)

plt.rcParams['figure.figsize'] = [6, 3]
fig = plt.figure()
plt.plot(r, p1, 'g', linewidth=1)
plt.plot(r, p2, 'r', linewidth=1)
plt.show()

# Create the joint probability matrix P
P = cp.Variable((n, n), nonneg=True)
# Objective: Maximize probability of R1 + R2 <= 0
r1p = r[:, np.newaxis]
r2p = r[np.newaxis, :]
loss_mask = (r1p + r2p <= 0)
objective = cp.Maximize(cp.sum(P[loss_mask]))

diff_r1 = r - mu1
diff_r2 = r - mu2

# Constraints
constraints = [
    cp.sum(P, axis=0) == p2,
    cp.sum(P, axis=1) == p1,
    cp.matmul((cp.matmul(diff_r1.T, P)), diff_r2) == rho * sigma1 * sigma2
]

# Problem
problem = cp.Problem(objective, constraints)
problem.solve(solver=cp.SCS, verbose=True)

# Results
p_max = problem.value
p_table = problem.var_dict['var1'].value
print(p_max)
print(diff_r1.T @ p_table @ diff_r2 - rho * sigma1 * sigma2)
print(np.sum(p_table))
print(np.sum(np.sum(p_table, axis=0)))
print(np.sum(np.sum(p_table, axis=1)))

# Smooth the results a bit because weird lines appear #
p_table[p_table < 5e-3] = 0
X, Y = np.meshgrid(r, r)
plt.axes().set_aspect('equal')
plt.contour(X, Y, p_table)
plt.plot(r, list(map(lambda x: max(x, -20), -r)), color="r", linewidth=0.3)  # R1+R2==0 red line
plt.show()
