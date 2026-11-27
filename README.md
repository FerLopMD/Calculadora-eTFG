# Calculadora eTFG / eGFR Calculator

Educational web application for adults implementing CKD-EPI 2021 creatinine, CKD-EPI 2012 cystatin C and CKD-EPI 2021 combined equations, with KDIGO 2024 categories and transparent intermediate calculations.

**Authors:** Abel Suárez-Castro, Carmen I. Bárcenas-López, Fernando Rosales-López, Genaro F. Sánchez-Mejorada.

## Run locally / Ejecutar localmente

Open the folder in VS Code and run in the terminal:

```powershell
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

For mathematical tests: `python -m unittest -v test_app.py`.

## Scope

Preliminary educational software, not clinically validated. No diagnosis, individual prognostic probabilities, or medication dosing. No database or application-level clinical-value logging. Do not enter identifying information. G categories use the rounded laboratory result, with a warning if the unrounded value crosses a threshold. ACR categories use the published unit-specific thresholds. Reference intervals of filtration markers are laboratory/method-specific. Cystatin C 2012 uses the original published age factor 0.996. No software DOI.

## References

- Inker et al. 2012. https://doi.org/10.1056/NEJMoa1114248
- Inker et al. 2021. https://doi.org/10.1056/NEJMoa2102953
- Miller et al. 2022. https://doi.org/10.1093/clinchem/hvab278
- KDIGO 2024. https://doi.org/10.1016/j.kint.2023.10.018
- NIDDK adult equations: https://www.niddk.nih.gov/research-funding/research-programs/kidney-clinical-research-epidemiology/laboratory/glomerular-filtration-rate-equations/adults

Version 1.0 (preliminary) — 2026-11-26.
