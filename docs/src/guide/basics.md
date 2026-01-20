# Mathematical background

Given a domain `$\Omega\subset\mathbb{R}^2$`, a triangular mesh `$\mathcal{T}$` of `$\Omega$` is constructed. In this mesh each edge has a degree assigned to it. Hence any given triangle `$T\in \mathcal{T}$` has three degrees `$p_1\le p_2 \le p_3$`, one per edge. We consider the space `$\mathcal{P}_{p_1 p_2 p_3}(T)$` formed by all bivariate polynomials over `$T$` of degree `$p_3$` such that restricted to the edge with degree `$p_i$` have degree `$p_i$` (`$i=1,2,3$`). Over the whole triangulation, we work with the space `$\mathcal{P}_\mathcal{T}$` formed by piecewise polynomials `$q$` such that `$q|_T \in \mathcal{P}_{p_1 p_2 p_3}$` where `$p_1$`, `$p_2$` and `$p_3$` are the degrees corresponding to the edges of `$T$`. For this space to be well defined the degrees ought to verify the _`$p$` confirmity_ condition:

```math
p_3\le p_1+p_2.
```

The method needs a problem given in weak form and an a posteriori error estimator `$\eta_T(u_h)$` to be computed from a numerical solution `$u_h$` over each triangle `$T$`. The procedure, then, is similar to a standard a posteriori method:

- A solution `$u_h$` is computed over `$\mathcal{T}$`.
- The estimator `$\eta_T(u_h)$` is computed for every `$T\in\mathcal{T}$`.
- Every `$T'\in\mathcal{T}$` such that `$\eta_{T'}(u_h)^2>\rho \frac{1}{|\mathcal{T}|}\sum_{T\in\mathcal{T}}\eta_T(u_h)$` is marked for refinement.
- An heuristic criterion is used to determine how `$T'$` will be refined. There are two possibilities:
  - `$T'$` is partitioned into 4 new triangles, by bisecting each edge (`$h$`-refinement).
  - The polynomial space associated to `$T'$` is enlarged by augmenting `$p_1\to p_1+1$` (`$p$`-refinement).
- The refinement of a triangle typically implies some refinement of neighbouring triangles, to mantain conformity. The `$h$`-conformity is mantained using a classical newest vertex bisection algorithm (red-blue-green variation). The `$p$` conformity is mantained by a recursive routine that augments some degrees when the `$p$`-conformity condition is not verified.
- The solver goes back to the first step, finding a new solution `$u_h$` over the new mesh (with its updated degrees). This is repeated until no triangle is marked for refinement or a maximal number of iterations is reached.

All parameters can be set by an `Options` structure.

