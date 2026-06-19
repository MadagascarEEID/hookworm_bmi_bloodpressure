# hookworm_bmi_bloodpressure

Analyses and figures associated with ms, "Hookworm parasitism shapes energy status and cardiovascular disease risk in northeast Madagascar"

**Authors**: Lev Kolinski, Georgia Titcomb, Kayla Kauffman, Jean Yves Rabezara, Prisca Rahary, Michelle Pender, Alma Solis, Camille M.M. DeSisto, Voahangy Soarimalala, Randall A. Kramer, James M. Moody, Charles L. Nunn

**hookworm_analysis_combined.Rmd** contains all code to:
- wrangle data
- generate descriptive statistics
- run averaged GLMs of BMI and blood pressure
- run SEMs to investigate direct and indirect associations of BMI, blood pressure, and hookworm infection
- generate figures

**hookworm_analysis_combined_R1.Rmd** contains all code to:
- wrangle data
- generate descriptive statistics
- run averaged GLMs of BMI and blood pressure
- run SEMs to investigate direct and indirect associations of BMI, blood pressure, and hookworm infection
- generate figures
- Updated for first revisions for AJHB (hence the R1 suffix)

**marginal_effects.R** contains all code to:
- generate marginal effects plots based on the averaged model of BMI
- population-averaged predictions per model, across wealth grid × infection group
- also plots observed values on top of marginal effects

**supp_anthropometrics_corr_test.R** contains all code to:
- run supplementary correlation tests for anthropometrics calculated in 2025
- generate correlation plot for supplementary materials

