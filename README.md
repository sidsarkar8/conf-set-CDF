# Unified Inference Framework for Functionals

This repository contains code related to the development of a **unified inference framework** for statistical functionals.

## Overview

This project aims to provide **assumption-flexible and post-selection-valid** inference methods for a broad class of statistical functionals (e.g., the mean, median). The core idea is to:

- Construct confidence sets for cumulative distribution functions (CDFs),
- Derive confidence intervals for functionals based on these CDF sets,
- Support user-specified assumptions (e.g., finite variance, bounded support, tail decay),
- Provide rigorous guarantees even under **post-selection inference settings**.

## Key Features

- **Assumption-flexible**: Choose assumptions best suited to your application (e.g., bounded moments, quantile constraints).
- **Functional-agnostic**: Works for a wide range of functionals beyond just the mean.
- **Post-selection validity**: Maintains coverage guarantees even when functionals are chosen based on data exploration.

## Applications

- Robust and reliable inference in exploratory data analysis.
- Data-driven uncertainty quantification with user-specified structural or distributional knowledge.
- Methodologically principled alternative to bootstrap and asymptotic approximations in sensitive or complex settings.
