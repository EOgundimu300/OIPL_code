# Observed-information penalization in binomial regression

R code and data for the article

Ogundimu, E. O. and Nadhrah, D. A. Finiteness and bias under observed-information penalization in binomial-response generalized linear models.

Table, figure and section numbers below refer to the article and its Supplementary Material.

## Files

| File | Contents |
| --- | --- |
| `obsinfo_core.R` | Likelihoods, penalties and adjusted scores, model fitting, separation checks, profile sets and exact enumeration for intercept-only models |
| `Endometrial_example.R` | Endometrial cancer analysis (Table 6) |
| `Endometrial_profile.R` | Profile likelihoods and intervals (Figure 2) |
| `two_group.R` | Two-group calculation (Table S13) |
| `final_checks.R` | Short deterministic checks of selected results |
| `obsinfo_simulation.R` | Fixed-design calculations and the simulation study |
| `profile.csv` | Profile-set summaries reported in Tables 5 and S11 |
| `Endometrial.dat` | Endometrial cancer data |

## Requirements

The code uses only packages included with R and was run with R 4.5.2. `Endometrial_profile.R` requires an R installation with Cairo support. Run all commands from the directory that contains the R files and `Endometrial.dat`.

## Data

`Endometrial.dat` contains the 79 endometrial cancer cases analysed by Heinze and Schemper (2002), as distributed with Agresti (2015), *Foundations of Linear and Generalized Linear Models*, at https://www.alanagresti.com/glm/data/Endometrial.dat. The columns are NV (neovascularization, 0 or 1), PI (pulsatility index), EH (endometrial height) and HG (histology grade, 1 for high grade). The same data are available as `endometrial` in the R package brglm2.

## Fitting a model

```r
source("obsinfo_core.R")
dat <- read.table("Endometrial.dat", header = TRUE)
fit <- obsinfo_fit(HG ~ NV + PI + EH, data = dat, link = "probit", method = "O")
summary(fit)
```

The argument `method` selects the observed-information penalty (`"O"`), the expected-information penalty (`"E"`), Firth's adjusted score (`"A"`) or maximum likelihood (`"ML"`). The available links are logit, probit, cloglog, loglog and cauchit.

## Endometrial cancer data (Section 5.4)

```sh
Rscript Endometrial_example.R Endometrial.dat
Rscript Endometrial_profile.R Endometrial.dat
```

The first command prints Table 6. The second prints the profile intervals and writes Figure 2 to `Endometrial_profile.pdf`.

## Fixed-design calculations

```sh
Rscript obsinfo_simulation.R tables
Rscript two_group.R
Rscript final_checks.R
```

The first command prints Tables 2, 3, S1 to S4, S6 and S7, the second prints Table S13, and the third runs the deterministic checks.

## Simulation study (Section 5.3 and Supplementary Section E)

```sh
Rscript obsinfo_simulation.R pilot
Rscript obsinfo_simulation.R production
Rscript obsinfo_simulation.R summary
```

Run the stages in this order; `Rscript obsinfo_simulation.R --help` lists the options. A second argument selects scenarios, for example `production 1:8`. Scenarios 1 to 8 use the probit link and 9 to 16 the complementary log-log link; within each link they run through n = 25, 50, 100 and 200 for setting S1 and then S2. Results, including the summaries for Tables 4, 5 and S8 to S12 and Figures 1 and S1, are written to `sim_out`. The number of cores, the output directory and the random-number seed are set in `cfg` at the top of the script.

## Citation

If you use this code, please cite the article above.
