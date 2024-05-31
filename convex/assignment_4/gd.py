import numpy as np
import matplotlib.pyplot as plt


def objective_function_calculation(x, A):
    """
    Calculates f(x) parametrized by A matrix (which contains the m different a_i's) and current x value.
    f(x) = - ∑_{i=1}^m log(1 - a_i^T * x) - ∑_{i=1}^n log(1 - x_i²)
    """
    return -np.sum(np.log(1 - A @ x)) - np.sum(np.log(1 - x ** 2))


def gradient_calculation(x, A):
    """
    Calculates gradient of f(x) parametrized by A at point x.
    ∇f(x) =  ∑_{i=1}^m (a_i / (1 - a_i^T * x))  +  ∑_{i=1}^n (2x_i / (1 - x_i²))
    """
    return A.T @ (1 / (1 - A @ x)) + (2 * x) / (1 - x ** 2)


def hessian_calculation(x, A):
    """
    Computes the Hessian of f(x) parametrized by A at point x.
    ∇²f(x) = ∑_{i=1}^m (a_i * a_i^T) / (1 - a_i^T * x)²  +  diag(2 / (1 - x_i²)²)
    """
    denominator = (1 - A @ x) ** 2
    outer_products = A[:, :, None] * A[:, None, :]
    diagonal_term = np.diag(2 / (1 - x ** 2) ** 2)
    return (outer_products / denominator[:, None, None]).sum(axis=0) + diagonal_term


def backtracking_line_search(x, A, direction, gradient, alpha, beta):
    """
    Calculates BLS given x, A, the gradient of the objective function at x, as well as parameters a and b.
    Checks 3 conditions:
    Check if the condition for ai^T * x < 1 holds.
    Check if the condition |x| < 1 holds.
    THe main BLS check
    """
    t = 1
    while True:
        x_new = x + t * direction
        if np.all(A @ x_new <= 1) and \
                np.all(np.abs(x_new) <= 1) and \
                objective_function_calculation(x_new, A) <= \
                objective_function_calculation(x, A) + alpha * t * np.dot(gradient, direction):
            break
        t *= beta
    return t


def gradient_descent_bls(A, x0, alpha, beta, eta, max_iter):
    x = x0
    obj_values = [objective_function_calculation(x, A)]
    step_sizes = [1]
    iterations = 0

    while iterations < max_iter:
        gradient = gradient_calculation(x, A)
        if np.linalg.norm(gradient, 2) ** 2 < eta:
            break

        direction = -gradient  # Descent direction is the negative gradient
        t = backtracking_line_search(x, A, direction, gradient, alpha, beta)
        x = x + t * direction

        obj_values.append(objective_function_calculation(x, A))
        step_sizes.append(t)
        iterations += 1

    # Calculate f(x^k) - p*, assuming that the last value is p*
    obj_values = [f - objective_function_calculation(x, A) for f in obj_values]

    return x, obj_values, step_sizes, None, iterations


def newton_method_bls(A, x0, alpha, beta, eta, max_iter):
    """Finds the analytic center using Newton's method with backtracking line search."""
    x = x0
    obj_values = [objective_function_calculation(x, A)]
    step_sizes = [1]
    newton_decrements = [np.inf]  # Initialize with a large decrement
    iterations = 0

    while iterations < max_iter:
        gradient = gradient_calculation(x, A)
        hessian = hessian_calculation(x, A)

        # Solve for Newton step using a linear system solver instead of inverse as we said #
        direction = np.linalg.solve(hessian, -gradient)  # ∆x = - (∇²f(x))⁻¹ * ∇f(x)
        newton_lambda_square = -np.dot(gradient, direction)  # λ²(x) = ∇f(x) * (∇²f(x))⁻¹ * ∇f(x)

        if newton_lambda_square <= eta:
            break

        t = backtracking_line_search(x, A, direction, gradient, alpha, beta)
        x = x + t * direction

        obj_values.append(objective_function_calculation(x, A))
        step_sizes.append(t)
        newton_decrements.append(np.sqrt(newton_lambda_square))
        iterations += 1

    # Calculate f(x^k) - p*, assuming that the last value is p*
    obj_values = [f - objective_function_calculation(x, A) for f in obj_values]

    return x, obj_values, step_sizes, newton_decrements, iterations


def run_method(name, **kwargs):
    if name == 'Gradient Descent':
        return gradient_descent_bls(**kwargs)
    elif name == 'Newton Method':
        return newton_method_bls(**kwargs)
    else:
        raise ValueError("Choose Gradient Descend or Newton Method!")


# --- Experiment Setup --- #
np.random.seed(1995)
alpha_values = [0.05, 0.25, 0.5]
beta_values = [0.1, 0.5, 0.9]
eta = 1e-6
max_iter = 2000
methods = ['Newton Method', 'Gradient Descent']

A = np.random.rand(200, 100)
x0 = np.zeros(100)

for method_name in methods:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for alpha in alpha_values:
        for beta in beta_values:
            x_bls, obj_values, step_sizes, newton_decrements, iterations = run_method(name=method_name,
                                                                                      A=A,
                                                                                      x0=x0,
                                                                                      alpha=alpha,
                                                                                      beta=beta,
                                                                                      eta=eta,
                                                                                      max_iter=max_iter
                                                                                      )

            # Plotting
            axes[0].semilogy(obj_values, label=f"α={alpha}, β={beta}")
            axes[1].plot(step_sizes, 'o-', markersize=3, label=f"α={alpha}, β={beta}")

    axes[0].set_xlabel('Iteration (k)')
    axes[0].set_ylabel('f(x^(k))')
    axes[0].set_title(f'Objective Function Value - {method_name}')

    axes[1].set_xlabel('Iteration (k)')
    axes[1].set_ylabel('Step Size (t)')
    axes[1].set_title(f'Step Size vs. Iteration - {method_name} ')
    axes[1].legend()
    plt.tight_layout()
    plt.show()
