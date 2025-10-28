# Cut-Based Piecewise-Affine (PWA) Approximation 
This repo contains codes related to cut-baed PWA approximation approach in my PhD research
These MATLAB functions construct and evaluate **piecewise-affine (PWA)** models that approximate nonlinear systems. 
---
## Requirements
- **MATLAB R2022b** or later
- **Optimization Toolbox**
- **Gurobi Optimizer** (for the parametric test script)
---
## Main Functions
### nldyn.m
Defines the **true nonlinear system** to approximate. 
Example implementation:
```
y = x(1) * cos(x(2));
```
You can modify this function to represent any nonlinear dynamics of interest.
---
### pwaapprox.m
Evaluates the **piecewise-affine model output** at a given input `x`.
1. Determines which region `x` belongs to, using the region constraints (`Ccal`).
2. Applies the corresponding affine law:
	f=Jpx+Kpf = J_p x + K_pf=Jp​x+Kp​
3. Returns the predicted output `f`.
---
### approxerr.m
Computes the **L2-norm of the total approximation error** between:
- The true nonlinear function (`nldyn`), and
- The PWA model (`pwaapprox`) over a set of test inputs `G`.
```
errn = norm(f_i - F_i, 2)
```
---
### continuity.m
Defines **nonlinear continuity constraints** for the optimizer `fmincon`.
These constraints ensure that the PWA model is **continuous across region boundaries**, preventing output jumps between adjacent affine pieces.
---
### hyperplanes.m
Generates the set of **hyperplanes** (`Hcal`) that divide the input space.
---
### chambers.m
Finds all **feasible chambers (regions)** formed by the intersection of hyperplanes that lie within the input domain.
- Enumerates all possible sign patterns (sides of each hyperplane).
- Uses **linear programming (linprog)** to test whether each combination corresponds to a nonempty region.
- Returns a matrix `Sigma`, whose columns represent the valid sign patterns.
---
### regions.m
Builds region-level data structures from the feasible chambers:
- `A`: the **adjacency matrix** indicating which regions share a boundary and through which hyperplane.
- `Ccal`: the **constraint matrix** defining the inequalities for each region.
This structure allows `pwaapprox` to determine region membership and ensures consistent geometry across the model.
---
### regapprox.m
The **main fitting routine** that learns the PWA approximation.
Steps:
1. Generates hyperplanes (`hyperplanes.m`),
2. Identifies feasible regions (`chambers.m`),
3. Builds adjacency and constraints (`regions.m`),
4. Solves a constrained optimization problem via `fmincon` to minimize the approximation error (`approxerr`),
5. Enforces continuity using `continuity.m`.
Returns either:
- The optimal affine parameters (`Jcal`, `Kcal`), or
- The minimized total error.
