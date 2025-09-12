# fast-poly
Fast arbitrary polynomial math

# Support
- Addition and Subtraction
- Multiplication (Fast Fourier Transform)
- Root finding (Halley's method with loss accumulator to detect rootlessness)
- Polynomial positional shift and translation across X axis
- Intersection

- All functions are chainable (with the exception of root finding as it returns the root, not the polynomial)

# Additional support
- Root finding of arbitrary single variable function using Secant Approximation

# TODO

- Powell's method for blackbox multivariate functions
- Bayesian optimization for expensive blackbox multivariate functions
