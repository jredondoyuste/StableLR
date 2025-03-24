[![arXiv:2502.18643](https://img.shields.io/badge/arXiv-2502.18643-b31b1b.svg)](https://arxiv.org/abs/2502.18643) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/jredondoyuste/StableLR/blob/main/LICENSE)
![GitHub last commit](https://img.shields.io/github/last-commit/jredondoyuste/StableLR)

# StableLR
Accompanying routines to study perturbative dynamics of nonlinear waves on stable lightrings. These codes were used in obtaining the results of https://arxiv.org/abs/2502.18643

![Snapshots from the evolution of the cubic defocusing wave equation on a spacetime with stable trapping](https://github.com/jredondoyuste/StableLR/blob/main/stableLR_illustration.png)

In particular, it includes:

- `linear_perturbations.nb`: Code to compute the QNM frequencies and radial functions for the potential with stable trapping, including direct integration, Breit-Wigner resonance method, and WKB approximation.
- `1d_nlkg.jl`: Julia language code that solves the cubic defocusing wave equation on a circle (1+1 dimensions), using a pseudospectral method, based on [`ApproxFun.jl`](https://github.com/JuliaApproximation/ApproxFun.jl).
- `2d_nlkg.jl` and `2d_nlkg_diss.jl`: Same as the above, but for the wave equation for axisymmetric modes on the sphere, and including artificial dissipation terms.


## Citing

If you find this code useful for your research, please consider citing the following paper:

```bibtex
@article{Redondo-Yuste:2025hlv,
    author = "Redondo-Yuste, Jaime and C\'ardenas-Avenda\~no, Alejandro",
    title = "{Perturbative and non-linear analyses of gravitational turbulence in spacetimes with stable light rings}",
    eprint = "2502.18643",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    month = "2",
    year = "2025"
}
```
