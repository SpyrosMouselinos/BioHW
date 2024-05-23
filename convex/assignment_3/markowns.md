## Solution
--------

### Theoretical Part

##### 1) Show that the objective function $f$ is convex - Note: The exersize says $f(S)$ but since our variable is P, I think it makes sense to refer to it as $f(P)$ instead

To show that the objective function $f(P)$ is convex, we first expand the expression for $d(x_i, y_i)$:
$d(x_i, y_i) = (x_i - y_i)^T P (x_i - y_i) ^{1/2}.$

Substitute this into the objective function $f(P)$:
$f(P) = \frac{1}{N} \sum_{i=1}^N ( d_i - \left( (x_i - y_i)^T P (x_i - y_i) \right)^{1/2})^2.$

Let $z_i = x_i - y_i$. Then, the objective function becomes:
$f(P) = \frac{1}{N} \sum_{i=1}^N (d_i - \left( z_i^T P z_i \right)^{1/2})^2.$

To prove convexity, we need to show that $f(P)$ is a convex function of $P$. By expanding the square we end up:
$f(P) = \frac{1}{N} \sum_{i=1}^N ( d_i^2 - 2 d_i (z_i^T P z_i)^{1/2} + z_i^T P z_i).$

We can ignore the $\frac{1}{N}$ for our analysis, end we end up with a sum of terms that look like:

$d_i^2 - 2 d_i (z_i^T P z_i)^{1/2} + z_i^T P z_i$

For each term:
- $d_i^2$ is constant with respect to P.
- $z_i^T P z_i$ is linear in $P$ because $P \in \mathbf{S}^n_+ $, and the product $z_i^T P z_i$ is a linear function of $P$.
- The square root of the affine function $z_i^T P z_i$ is concave and multiplying by a negative constant $-2 d_i$ makes it convex. $-2 d_i (z_i^T P z_i)^{1/2}$ is convex.

Therefore, each term of the summation is a sum of a constant, a linear term, and a convex term, thus a convex term. The sum of convex terms is convex, which makes $f(P)$ a convex function.

### Theoretical Part

##### 2) Show that the convex program can be expressed by an equivalent conic program with linear objective and a number of conic constraints using the $R^n_+$ (nonnegative orthant cone), $Q^n$ (second order cone), $Q_r^n$ (rotated second order cone), $S^n_+$ (positive semidefinite cone).

#### Conversion to conic problem with linear objective

We can start by defining  $t_i$ to represent $\sqrt{((x_i - y_i)^T P (x_i - y_i))}$:

- $t_i = \sqrt{(x_i - y_i)^T P (x_i - y_i)}$

Then we can rewrite the objective function using slack variables $t_i$ and $u_i$ to represent $((d_i - t_i)^2)$:
- $u_i = (d_i - t_i)^2 $

Our goal is to minimize:
$f(u) = \frac{1}{N} \sum_{i=1}^N u_i.$ or the equivalent $f(u) = \sum_{i=1}^N u_i.$ or even $f(u) = \quad 1^T \mathbf{u}.$ with $\mathbf{u}$ a N-dimensional vector with N the size of the datapoints.

#### Constraints

1. **Second-Order Cone (SOC) Constraint for $t_i$**:
   $$ t_i = \sqrt{(x_i - y_i)^T P (x_i - y_i)}$$
   $$ t_i^2 = (x_i - y_i)^T P (x_i - y_i)$$
   or as a constraint:
   $$ t_i^2 \geq (x_i - y_i)^T P (x_i - y_i), \quad \forall i $$
   $$
   \| \begin{pmatrix}P^{1/2} (x_i - y_i) \end{pmatrix}\|_2 \leq t_i
   $$
   $$
   (t_i, P^{1/2} (x_i - y_i)) \in Q^n
   $$
   where
   $$
   Q^n = \{(t, x) \mid \|x\|_2 \leq t\}
   $$
   Another way to express this would also be to directly employ the fact that P is a PSD matrix.
   Then we know that it can be written as $P = L^TL$, where $L$ is the Cholesky factorization (or any equivalent one)
   Thus:
    $$
    (x_i - y_i)^T P (x_i - y_i)=(x_i - y_i)^T L^T L (x_i - y_i)=\|L (x_i - y_i)\|_2^2
    $$
   and:
   $$
   t_i^2 \geq \|L (x_i - y_i)\|_2^2
   $$
   $$
   (t, L (x_i - y_i)) \in Q^n
   $$


2. **Rotated Second-Order Cone (RSOC) Constraint for $u_i$**:
   We can express $(d_i - t_i)^2 \leq u_i$ as the equivalent $\| (d_i - t_i) \|^2_2 \leq u_i$:
   $$
   \| \begin{pmatrix} 2\sqrt{u_i} \\ d_i - t_i \end{pmatrix} \|_2 \leq d_i + t_i, \quad \forall i
   $$

   $$
   (1/2, u_i, d_i - t_i) \in Q_r^n, u \geq 0, v \geq 0
   $$

3. **Positive Semidefinite Cone (PSD) Constraint**:
   $$
   P \succeq 0
   $$ 
   or 
   $$
   P \in S^n_+
   $$

4. **Non-Negativity Constraint**:
   $$ 
   u_i \geq 0, \quad t_i \geq 0, \quad \forall i
   $$
   or 
   $$
   u_i, t_i \in R^n_+
   $$

The final conic program is:

$$
\text{minimize} \sum \mathbf{u_i}
$$

subject to:

$$
(t_i, P^{1/2} (x_i - y_i)) \in Q^n \quad \text{(SOC)}
$$

$$
(1/2, u_i, d_i - t_i) \in Q_r^n \quad \text{(RSOC)}
$$

$$
P \in S^n_+
$$

$$
u_i \in R^n_+, t_i \in R^n_+, \forall i \quad \text{(Non-Negativity)}
$$



