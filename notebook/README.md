To run Python script with uv:
```bash
uv run python ODESLV_tutorial.py
```
To lauch Jupyter Lab:
```bash
uv run --with jupyter jupyter lab
```

The tutorials:

| notebook | what it shows |
|---|---|
| `FFModel_tutorial` | declaring a model: domains, states, inputs, equations, outputs, and its analysis |
| `ODESLV_tutorial` | integrating an ODE model with SUNDIALS CVODES, with forward and adjoint sensitivities |
| `DAESLV_tutorial` | an index-2 DAE with SUNDIALS IDAS: index reduction and consistent initialisation |
| `OCFESLV_tutorial` | a nonlinear PDE by space-time collocation, monolithic and marching |
| `MOL_tutorial` | the same PDE by the method of lines on `DAESLV`, compared with `OCFESLV` |
| `PSA_tutorial` | breakthrough in a non-isothermal two-component adsorption column, by either solver |

`MOL_tutorial` and `PSA_tutorial` use CRONOS features newer than the 5.1.0 release (sensitivities by
name on `OCFESLV`, the collocation options at `options.INTERFACE` on both solvers, nonlinear outputs
under the method of lines): run them against CRONOS built from this tree until the next release.
