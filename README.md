# VADERS

**Vertical Atmospheric Dispersion with dEposition, Resuspension, and Surface Inactivation**

This repository contains the source code used to generate the simulations presented in the manuscript describing a vertically resolved mechanistic model for atmospheric dispersion with a reactive air-surface exchange boundary. The model includes deposition, resuspension, and surface inactivation processes.

## Requirements

- GNU Octave (tested with version 11.3.0)

## Repository Structure

### vaders_main

Input:
- Ratio $\Lambda^*_{\mathrm{res}}/v^*_{\mathrm{dep}}$ (dimensionless resuspension-rate coefficient divided by dimensionless deposition velocity).

Output:
- Air-clearance times for selected Peclet numbers and surface-inactivation rates.

Description:
- Simulations assuming a constant turbulent diffusivity.

### vaders_var_diff

Input:
- Ratio $\Lambda^*_{\mathrm{res}}/v^*_{\mathrm{dep}}$.

Output:
- Air-clearance times for selected Peclet numbers and surface-inactivation rates.

Description:
- Simulations using a linear parameterization of the turbulent diffusivity.

### vaders_var_lambda

Input:
- Number of hours before the resuspension rate transitions from a low value to a high value.

Description:
- The transition is represented using a logistic decay function.
- The specified time corresponds to the inflection point of the logistic function.

Output:
- Air-clearance times and concentration profiles for time-dependent resuspension scenarios.

### verification_vaders

Contains the verification framework for the numerical implementation.

Verification includes:

- Comparison against analytical solutions for two reference cases:
  - Case (i): perfect deposition
  - Case (ii): zero deposition
 
Input:
- `D` : turbulent diffusivity
- `N` : number of finite-volume cells

Output:
- Relative $L_2$ errors

- Mass-conservation verification of the full model.

Input:
- Ratio $\Lambda^*_{\mathrm{res}}/v^*_{\mathrm{dep}}$.

Output:
- total, and initial mass in system with respect to time

## Running the code

All the scripts are run as functions which can be included in ensemble runs for creating figures for different input arguments.

