[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10675749.svg)](https://doi.org/10.5281/zenodo.10675749)
[![Testing suite](https://github.com/irukoa/Extrapolation_Integration/actions/workflows/CI.yml/badge.svg)](https://github.com/irukoa/Extrapolation_Integration/actions/workflows/CI.yml)
# Extrapolation Integration

This is a tiny modern Fortran library to implement the [Richardson extrapolation](https://en.wikipedia.org/wiki/Richardson_extrapolation) method for integration. For more information, see Ref. [[1]](#ref1).

# API

The library adds the `extrapolation(array)` array reduction function. This is used as
```fortran
result = extrapolation(array)
```
where the type of `array(:)` is any of `real(sp)`, `complex(sp)`, `real(dp)` or `complex(dp)` and `result` is an scalar of the corresponding type.

If `array(:)`, of size $N$, holds a sampling of $f(x)$ in the range $x\in[a, b]$ discretized according to `array(i)` $=f(x_i)$, where

$$
x_i = a + (b-a) \frac{i-1}{N-1},
$$

then, on output, `result` represents the unnormalized integral of $f(x)$, i.e.,

$$
\frac{1}{b - a}\int_{a}^{b}f(x)dx.
$$

Three cases are distinguished,
- $N$ can be expressed as $2^M + 1$ for some $M\in\mathbb{N}$: the extrapolation method is employed, and the expected accuracy is $\mathcal{O}(h^{2M})$, where $h = (b-a)/2$.
- $N$ cannot be expressed as $2^M + 1$: the [trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule) is employed and the expected accuracy is $\mathcal{O}(h^{2})$, where $h = (b-a)/(N-1)$.
- $N = 1$: no integration and `result = array(1)`.

## Interface
```fortran
function extrapolation(array)
  ${type}$, intent(in) :: array(:)
  ${type}$ :: extrapolation
end function extrapolation
```

# Working principle

The extrapolation method computes an approximation to

$$
\frac{1}{b - a}\int_{a}^{b}f(x)dx,
$$

by employing relatively inaccurate trapezoidal rule approximations. The method works by averaging these approximations such that error elimination is achieved iteratively [[1]](#ref1). This requires, that for a final result accurate up to $\mathcal{O}(h^{2M})$, $M$ trapezoidal rule approximations be computed, each with step size $h_i = (b-a)/2^i$.

The implementation takes advantage of the fact that if the number of sampling points can be expressed as $2^M+1$, then, only a single sampling of $f(x)$ with step size $h_M = (b-a)/2^M$ is required to compute all $M$ trapezoidal rule approximations. The library reorders the data of the discretization of $f(x)$ and computes iteratively an extrapolation reduction.

# Build

An automated build is available for [Fortran Package Manager](https://fpm.fortran-lang.org/) users. This is the recommended way to build and use extrapolation integration in your projects. You can add extrapolation integration to your project dependencies by including

```
[dependencies]
Extrapolation_Integration = { git="https://github.com/irukoa/Extrapolation_Integration.git" }
```
to the `fpm.toml` file.

<a id="ref1"></a>
[1] R. L. Burden y J. D. Faires, Numerical analysis, 9. ed., International ed. Belmont, Calif.: Brooks/Cole, 2011, pp. 213-220.
