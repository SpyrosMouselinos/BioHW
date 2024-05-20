import cvxpy as cp
import numpy as np

# Define the size of the matrix
n = 3

# Example PSD matrix P
P = np.array([[4, 1, 1],
              [1, 2, 1],
              [1, 1, 3]])

# Define the variable X
X = cp.Variable((n, n), symmetric=True)

# Define the constraints
constraints = [X >> 0, cp.bmat([[X, P], [P, X]]) >> 0]

# Define the problem
# Since we are only interested in finding X, we can minimize a constant (e.g., 0)
objective = cp.Minimize(0)
problem = cp.Problem(objective, constraints)

# Solve the problem
problem.solve()

# Print the result
print("Square root of P (X):")
print(X.value)
