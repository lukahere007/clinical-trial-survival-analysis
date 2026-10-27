# Statistical Analysis Plan

**Study:** Reanalysis of the NCCTG / Intergroup adjuvant chemotherapy trial in resected stage III colon cancer (`survival::colon`)
**Version:** 1.0 (fixed before any outcome analysis was run in this repository; see the note in the README about what "pre-specified" can and cannot mean for a public dataset)
**Author:** Luke Wamalwa

---

## 1. Objectives

| | Objective |
|---|---|
| **Primary** | Estimate the effect of levamisole + 5-fluorouracil (Lev+5FU) versus observation (Obs) on overall survival. |
| **Key secondary** | Estimate the effect of Lev+5FU versus Obs on recurrence-free survival. |
| **Secondary** | Estimate the effect of levamisole alone (Lev) versus Obs on overall and recurrence-free survival. |
| **Supportive** | Describe time to recurrence with death without recurrence treated as a competing event. |

## 2. Endpoints

Time origin is randomisation. The dataset records time in days; results are reported in years (days / 365.25).

| Endpoint | Event | Censoring |
|---|---|---|
| **Overall survival (OS)** – primary | Death from any cause | Alive at last follow-up |
| **Recurrence-free survival (RFS)** – key secondary | Recurrence or death from any cause, whichever comes first | Alive and recurrence-free at last follow-up |
| **Time to recurrence (TTR)** – supportive | Recurrence | Last follow-up; death without recurrence is a competing event (Section 6.6) |

## 3. Estimands (plain words)

- **Population:** patients with resected stage III colon cancer eligible for the trial.
- **Treatment condition:** being *assigned* to Lev+5FU (or Lev) versus being assigned to observation, regardless of how much treatment was actually received (treatment-policy strategy). Dose modifications and discontinuation are part of the assigned policy.
- **Variable:** time from randomisation to the event defined in Section 2.
- **Intercurrent events:** treatment discontinuation, handled by treatment policy (ignored). For RFS, death before recurrence is part of the endpoint (composite strategy). For TTR, death before recurrence is a competing event (Section 6.6).
- **Population-level summaries:**
  1. Hazard ratio (Cox model), interpreted as an average hazard ratio over follow-up if proportional hazards do not hold.
  2. Difference in restricted mean survival time (RMST) up to 5 years: the extra event-free time, on average, gained in the first 5 years after randomisation. This does not rely on proportional hazards.

## 4. Analysis population

**Intention-to-treat (ITT):** all randomised patients, analysed in the arm they were assigned to. No per-protocol population can be defined because the dataset contains no adherence information.

## 5. Multiplicity

Four confirmatory hypotheses are tested in a fixed sequence, each at two-sided α = 0.05. A hypothesis is formally tested only if every hypothesis before it was rejected; once one fails, all later ones are reported descriptively (nominal p-values only). This controls the family-wise error rate at 0.05.

1. **H1:** OS, Lev+5FU vs Obs
2. **H2:** RFS, Lev+5FU vs Obs
3. **H3:** OS, Lev vs Obs
4. **H4:** RFS, Lev vs Obs

All other analyses (adjusted models, sensitivity analyses, subgroups of covariates) are supportive and carry no multiplicity adjustment.

## 6. Methods

### 6.1 Baseline characteristics
Summarised by arm: n (%) for categorical variables, median (IQR) for continuous ones, with missing counts shown. No hypothesis tests on baseline differences (randomisation makes them uninformative).

### 6.2 Primary test
Two-sided log-rank test comparing the two arms in each hypothesis (patients in the third arm excluded from that comparison). Unstratified.

### 6.3 Effect estimates
- **Kaplan–Meier** curves by arm with numbers at risk; median survival with 95% CI (Brookmeyer–Crowley, log transform as in `survival::survfit`); 5-year event-free probability with 95% CI.
- **Unadjusted Cox model** with arm as a three-level factor (reference Obs). HR with Wald 95% CI.
- **Adjusted Cox model** (supportive), covariates fixed in advance on clinical grounds:
  age (linear), sex, more than 4 positive nodes (`node4`), extent of local spread (submucosa/muscle, serosa, contiguous structures), differentiation (well, moderate, poor), obstruction, perforation, adherence to nearby organs, time from surgery to registration (`surg`).
  `node4` is used rather than the raw node count because the raw count is missing for some patients and the two variables are known to disagree in a few records.
- **RMST** up to τ = 5 years for each arm, computed by integrating the Kaplan–Meier curve, with a Greenwood-type standard error. Difference between arms with 95% CI and a z-test.

### 6.4 Proportional-hazards assessment
Scaled Schoenfeld residuals (`cox.zph`) for each term in the unadjusted and adjusted models, plus a plot of the residuals against time for the treatment term. A term is flagged if p < 0.05 *and* the plot shows a visible trend.

Pre-specified response:
- **Covariate violates PH:** refit the adjusted model stratified on that covariate (it keeps its own baseline hazard). Report the treatment HR from the stratified model as a sensitivity analysis.
- **Treatment violates PH:** the HR is reported as an average over follow-up and the RMST difference becomes the main effect summary. A piecewise HR (0–2 years, beyond 2 years) is reported descriptively.

The piecewise HR is shown in all cases as a descriptive check.

### 6.5 Missing data
Outcomes and arm are complete. Baseline covariates may be missing.
- The unadjusted (primary) analyses use all ITT patients and are unaffected.
- The adjusted Cox model is fitted to complete cases (primary adjusted analysis).
- **Sensitivity:** a missing-indicator category for each partially missing categorical covariate. In a randomised trial this keeps all patients and gives an unbiased treatment effect, because baseline covariates are independent of arm by design.
Missing counts are reported by variable and arm.

### 6.6 Competing risks (supportive)
For TTR, cumulative incidence of recurrence and of death without recurrence is estimated with the Aalen–Johansen estimator by arm. The cause-specific HR for recurrence (Cox model censoring at death without recurrence) is reported. No subdistribution (Fine–Gray) model is planned.

### 6.7 Sensitivity analyses (summary)

| # | Analysis | Purpose |
|---|---|---|
| S1 | Adjusted Cox, complete cases | Precision gain; check for chance imbalance |
| S2 | Adjusted Cox, missing-indicator | Missing covariate data |
| S3 | Adjusted Cox stratified on PH-violating covariates | Non-proportional covariate effects |
| S4 | Piecewise treatment HR (0–2 y, > 2 y) | Non-proportional treatment effect |
| S5 | RMST at τ = 7 years | Dependence of RMST on τ |

## 7. Software
R (≥ 4.3) with `survival` for all estimation. Code in `R/` is unit-tested with `testthat`.
