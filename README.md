# Piecewise-Affine (PWA) Approximation Toolbox

This MATLAB toolbox constructs and evaluates **piecewise-affine (PWA)** models that approximate nonlinear systems.  
It partitions the input space into polyhedral regions using hyperplanes, fits local affine models, and ensures continuity between them.

---

## 🔧 Core Functions

### `nldyn.m`
Defines the **true nonlinear system** to approximate.  
Example implementation:
```matlab
y = x(1) * cos(x(2));