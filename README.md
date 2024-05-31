# BioHW

Rewrite this with proper way that can render in the window of this chat
3. **Norm Interpretation**

   The supremum of a linear function over a unit ball occurs in the direction of the gradient. Thus:
   $\lambda(x) = \|\nabla^2 f(x)^{-1/2} \nabla f(x)\|$

   Since the Euclidean norm of a vector \( y \) is defined as \( \|y\| = \sqrt{y^T y} \), we have:
   \[ \lambda(x) = \left( \nabla f(x)^T \nabla^2 f(x)^{-1} \nabla f(x) \right)^{1/2}. \]

This confirms the first expression. 

### Transition to the Second Expression

The second expression:
\[ \lambda(x) = \sup_{v \neq 0} \frac{-v^T \nabla f(x)}{(v^T \nabla^2 f(x) v)^{1/2}}, \]
follows from the first expression by relaxing the constraint \( v^T \nabla^2 f(x) v = 1 \) to \( v \neq 0 \).

1. **Relaxing the Constraint**

   For any non-zero \( v \), consider the normalized vector:
   \[ u = \frac{v}{(v^T \nabla^2 f(x) v)^{1/2}}. \]

   This normalization ensures \( u^T \nabla^2 f(x) u = 1 \), and therefore:
   \[ \lambda(x) = \sup_{u^T \nabla^2 f(x) u = 1} (-u^T \nabla f(x)). \]

2. **Generalizing to Non-Zero \( v \)**

   The expression can be rewritten using the general \( v \neq 0 \) as:
   \[ \lambda(x) = \sup_{v \neq 0} \frac{-v^T \nabla f(x)}{(v^T \nabla^2 f(x) v)^{1/2}}. \]

This shows that the second expression follows immediately from the first by considering the appropriate normalization of \( v \).