# Import packages.
import cvxpy as cp
import numpy as np

import cvxpy as cp
import numpy as np

# Define the dimensions of the matrix
n_rows, n_cols = 3, 4

# Create a variable matrix of size n_rows x n_cols
matrix = cp.Variable((n_rows, n_cols))

# Known boolean mask (as a 2D array)
mask = np.array([
    [1, 0, 0, 1],
    [0, 1, 1, 0],
    [1, 0, 0, 1]
], dtype=bool)

# Objective function: minimize the sum of the masked elements
objective = cp.Minimize(cp.sum(cp.multiply(mask, matrix)))

# Constraints
constraints = [
    matrix >= 0,  # All elements must be non-negative
    cp.sum(matrix, axis=1) == [10, 15, 12]  # Sum of rows must match these values
]

# Define and solve the problem
problem = cp.Problem(objective, constraints)
problem.solve(solver=cp.GUROBI)

# Output the results
print("The optimal matrix is:")
print(matrix.value)
print("Minimum sum of selected elements:", problem.value)
