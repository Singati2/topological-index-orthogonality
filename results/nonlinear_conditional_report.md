# Conditional non-linear association test (scripts/73_nonlinear_conditional.py)

Question: does candidate z carry target information beyond the 30-index baseline X under NON-LINEAR dependence?

## Method

- Cross-fitted (double-ML) residualization: 5-fold KFold; within each fold RandomForestRegressor(300 trees, min_samples_leaf=5, seed=20261005) of y on X and of z on X fitted on training folds, predicted on the held-out fold; ry = y - yhat, rz = z - zhat (z standardized first). BBBP: y in {0,1}, RF regressor gives a probability residual.
- rho_nl = Pearson(rz, ry), Fisher-z 95% CI with se = 1/sqrt(n-3), all molecules.
- dcor = distance correlation(rz, ry) (V-statistic, from scratch, O(n^2) memory), permutation p-value with 199 permutations (minimum attainable p = 0.005). **For the dCor part only, Lipophilicity and BBBP (n > 2000) are subsampled to 2000 molecules with a fixed seed**; rho_nl uses all molecules.
- linear_pcor = Pearson of OLS residuals of y and z on [1, X] (in-sample).
- BH adjustment across all real candidates x datasets (reference rows excluded from the family).
- Reference rows: REF_invNirmala (sum 1/sqrt(du+dv); exactly in linear span of X), REF_noise (N(0,1) column; calibration), REF_leak (z = y + N(0, sd(y)); positive control).

Caveat: a candidate's non-linear residual rz is exactly zero only if the RF reproduces z perfectly; RF cannot extrapolate/represent exact linear identities, so an in-span index has rz != 0 (nuisance error). Cross-fitting makes this error independent of the held-out y noise, so rho_nl stays centred near 0 unless the nuisance error itself tracks ry.

## esol (n = 1127, n_dcor = 1127)

```
     candidate linear_pcor rho_nl  ci_lo  ci_hi   dcor perm_p  bh_q            note
     InfoH_deg      +0.154 +0.142 +0.085 +0.199 +0.120  0.005 0.018                
    InfoH_dist      -0.049 +0.009 -0.049 +0.068 +0.053  0.745 0.745                
    InfoH_spec      -0.009 +0.027 -0.031 +0.086 +0.056  0.560 0.605                
     InfoH_ecc      -0.039 +0.038 -0.020 +0.096 +0.073  0.090 0.176                
    Bonchev_Id      -0.113 +0.015 -0.043 +0.074 +0.066  0.130 0.236                
       Btw_sum      -0.073 -0.018 -0.076 +0.041 +0.063  0.430 0.518                
       Cls_sum      -0.002 -0.035 -0.093 +0.024 +0.070  0.145 0.242                
       Eig_sum      +0.001 -0.058 -0.116 +0.000 +0.096  0.015 0.048                
     Harm_cent      +0.000 +0.052 -0.007 +0.110 +0.073  0.070 0.147 in_span(linear)
       Btw_var      +0.026 +0.007 -0.051 +0.066 +0.077  0.050 0.114                
       AlgConn      +0.037 -0.001 -0.059 +0.058 +0.059  0.620 0.644                
     Rho_x_Dia      -0.077 -0.004 -0.063 +0.054 +0.061  0.370 0.469                
        EE_x_W      -0.056 +0.039 -0.019 +0.097 +0.071  0.185 0.279                
     Triangles      -0.070 +0.003 -0.055 +0.062 +0.048  0.595 0.626                
     MeanClust      +0.039 +0.015 -0.043 +0.074 +0.050  0.705 0.719                
       FourCyc      +0.112 +0.001 -0.058 +0.059 +0.049  0.555 0.605                
      CutVerts      -0.118 -0.015 -0.074 +0.043 +0.069  0.170 0.267                
        TotEcc      -0.073 +0.005 -0.053 +0.063 +0.073  0.120 0.223                
        Radius      -0.070 -0.049 -0.107 +0.009 +0.066  0.200 0.296                
       MeanEcc      -0.109 -0.028 -0.087 +0.030 +0.060  0.490 0.568                
REF_invNirmala      +0.000 +0.043 -0.015 +0.101 +0.070  0.185     - in_span(linear)
     REF_noise      -0.038 -0.009 -0.067 +0.049 +0.050  0.560     -                
      REF_leak      +0.623 +0.576 +0.535 +0.613 +0.518  0.005     -                
```

## freesolv (n = 639, n_dcor = 639)

```
     candidate linear_pcor rho_nl  ci_lo  ci_hi   dcor perm_p  bh_q            note
     InfoH_deg      -0.061 -0.104 -0.181 -0.027 +0.108  0.025 0.067                
    InfoH_dist      +0.086 -0.042 -0.120 +0.035 +0.091  0.145 0.242                
    InfoH_spec      +0.080 -0.030 -0.108 +0.047 +0.105  0.060 0.133                
     InfoH_ecc      +0.026 -0.060 -0.137 +0.018 +0.087  0.230 0.307                
    Bonchev_Id      -0.039 +0.006 -0.071 +0.084 +0.077  0.535 0.594                
       Btw_sum      +0.088 +0.003 -0.074 +0.081 +0.087  0.225 0.305                
       Cls_sum      -0.065 -0.029 -0.107 +0.048 +0.081  0.375 0.469                
       Eig_sum      +0.059 +0.029 -0.049 +0.106 +0.086  0.220 0.305                
     Harm_cent      +0.000 +0.001 -0.077 +0.078 +0.082  0.440 0.518 in_span(linear)
       Btw_var      -0.003 +0.006 -0.072 +0.083 +0.110  0.020 0.055                
       AlgConn      -0.094 -0.003 -0.081 +0.075 +0.109  0.020 0.055                
     Rho_x_Dia      -0.002 -0.095 -0.172 -0.018 +0.090  0.160 0.256                
        EE_x_W      +0.129 +0.027 -0.050 +0.105 +0.081  0.435 0.518                
     Triangles      +0.107 +0.011 -0.067 +0.089 +0.065  0.520 0.586                
     MeanClust      +0.051 +0.020 -0.058 +0.097 +0.085  0.180 0.277                
       FourCyc      -0.074 +0.033 -0.044 +0.111 +0.061  0.585 0.624                
      CutVerts      +0.044 +0.014 -0.064 +0.091 +0.086  0.215 0.305                
        TotEcc      -0.018 -0.029 -0.106 +0.049 +0.093  0.140 0.242                
        Radius      +0.021 -0.001 -0.078 +0.077 +0.068  0.710 0.719                
       MeanEcc      +0.061 -0.027 -0.104 +0.051 +0.083  0.340 0.446                
REF_invNirmala      +0.000 -0.004 -0.082 +0.073 +0.079  0.580     - in_span(linear)
     REF_noise      -0.020 -0.030 -0.107 +0.048 +0.068  0.500     -                
      REF_leak      +0.668 +0.661 +0.615 +0.703 +0.585  0.005     -                
```

## lipophilicity (n = 4200, n_dcor = 2000)

```
     candidate linear_pcor rho_nl  ci_lo  ci_hi   dcor perm_p  bh_q            note
     InfoH_deg      -0.045 -0.013 -0.043 +0.018 +0.062  0.020 0.055                
    InfoH_dist      +0.042 -0.038 -0.068 -0.007 +0.057  0.090 0.176                
    InfoH_spec      +0.077 -0.026 -0.057 +0.004 +0.086  0.005 0.018                
     InfoH_ecc      +0.022 -0.026 -0.056 +0.005 +0.063  0.020 0.055                
    Bonchev_Id      -0.011 -0.049 -0.079 -0.019 +0.078  0.005 0.018                
       Btw_sum      -0.016 -0.035 -0.065 -0.004 +0.064  0.035 0.085                
       Cls_sum      -0.026 +0.006 -0.024 +0.036 +0.070  0.005 0.018                
       Eig_sum      +0.014 -0.020 -0.050 +0.011 +0.062  0.030 0.075                
     Harm_cent      +0.000 -0.047 -0.077 -0.017 +0.077  0.005 0.018 in_span(linear)
       Btw_var      -0.002 +0.022 -0.008 +0.052 +0.075  0.005 0.018                
       AlgConn      -0.055 -0.031 -0.061 -0.000 +0.085  0.005 0.018                
     Rho_x_Dia      +0.001 -0.032 -0.062 -0.002 +0.051  0.160 0.256                
        EE_x_W      +0.068 -0.038 -0.068 -0.007 +0.078  0.005 0.018                
     Triangles      +0.013 +0.003 -0.027 +0.033 +0.043  0.360 0.465                
     MeanClust      -0.006 -0.001 -0.031 +0.029 +0.043  0.400 0.492                
       FourCyc      +0.016 -0.034 -0.065 -0.004 +0.038  0.500 0.571                
      CutVerts      -0.008 -0.026 -0.056 +0.004 +0.056  0.070 0.147                
        TotEcc      -0.018 -0.046 -0.076 -0.016 +0.067  0.015 0.048                
        Radius      -0.014 -0.037 -0.067 -0.007 +0.047  0.225 0.305                
       MeanEcc      +0.012 -0.028 -0.058 +0.002 +0.048  0.225 0.305                
REF_invNirmala      +0.000 -0.044 -0.074 -0.014 +0.078  0.005     - in_span(linear)
     REF_noise      +0.040 +0.039 +0.009 +0.069 +0.052  0.085     -                
      REF_leak      +0.683 +0.636 +0.618 +0.654 +0.587  0.005     -                
```

## bbbp (n = 2039, n_dcor = 2000)

```
     candidate linear_pcor rho_nl  ci_lo  ci_hi   dcor perm_p  bh_q            note
     InfoH_deg      +0.006 +0.040 -0.003 +0.084 +0.071  0.010 0.035                
    InfoH_dist      -0.016 -0.024 -0.067 +0.020 +0.062  0.030 0.075                
    InfoH_spec      +0.014 -0.004 -0.048 +0.039 +0.107  0.005 0.018                
     InfoH_ecc      -0.033 -0.000 -0.044 +0.043 +0.052  0.135 0.240                
    Bonchev_Id      +0.038 -0.025 -0.068 +0.019 +0.085  0.005 0.018                
       Btw_sum      +0.010 -0.018 -0.061 +0.026 +0.086  0.005 0.018                
       Cls_sum      +0.004 -0.010 -0.053 +0.034 +0.109  0.005 0.018                
       Eig_sum      +0.017 +0.000 -0.043 +0.044 +0.107  0.005 0.018                
     Harm_cent      +0.000 -0.029 -0.072 +0.014 +0.082  0.005 0.018 in_span(linear)
       Btw_var      +0.033 +0.016 -0.028 +0.059 +0.053  0.090 0.176                
       AlgConn      +0.006 -0.001 -0.044 +0.043 +0.090  0.005 0.018                
     Rho_x_Dia      -0.006 -0.002 -0.045 +0.041 +0.097  0.005 0.018                
        EE_x_W      +0.012 -0.026 -0.070 +0.017 +0.088  0.005 0.018                
     Triangles      +0.225 +0.013 -0.031 +0.056 +0.051  0.100 0.190                
     MeanClust      +0.078 +0.008 -0.035 +0.051 +0.060  0.040 0.094                
       FourCyc      -0.220 -0.171 -0.212 -0.128 +0.218  0.005 0.018                
      CutVerts      +0.002 -0.015 -0.058 +0.029 +0.102  0.005 0.018                
        TotEcc      +0.035 -0.019 -0.063 +0.024 +0.120  0.005 0.018                
        Radius      -0.012 -0.036 -0.079 +0.007 +0.084  0.005 0.018                
       MeanEcc      +0.001 -0.012 -0.055 +0.032 +0.087  0.005 0.018                
REF_invNirmala      +0.000 -0.029 -0.072 +0.015 +0.076  0.010     - in_span(linear)
     REF_noise      -0.008 -0.026 -0.069 +0.017 +0.039  0.485     -                
      REF_leak      +0.617 +0.567 +0.537 +0.596 +0.508  0.005     -                
```

## Flagged real candidates (bh_q < 0.05 or |rho_nl| >= 0.10)

```
      dataset  candidate  linear_pcor  rho_nl  ci_lo  ci_hi  dcor  perm_p  bh_q
         esol  InfoH_deg        0.154   0.142  0.085  0.199 0.120   0.005 0.018
         esol    Eig_sum        0.001  -0.058 -0.116  0.000 0.096   0.015 0.048
     freesolv  InfoH_deg       -0.061  -0.104 -0.181 -0.027 0.108   0.025 0.067
lipophilicity InfoH_spec        0.077  -0.026 -0.057  0.004 0.086   0.005 0.018
lipophilicity Bonchev_Id       -0.011  -0.049 -0.079 -0.019 0.078   0.005 0.018
lipophilicity    Cls_sum       -0.026   0.006 -0.024  0.036 0.070   0.005 0.018
lipophilicity  Harm_cent        0.000  -0.047 -0.077 -0.017 0.077   0.005 0.018
lipophilicity    Btw_var       -0.002   0.022 -0.008  0.052 0.075   0.005 0.018
lipophilicity    AlgConn       -0.055  -0.031 -0.061 -0.000 0.085   0.005 0.018
lipophilicity     EE_x_W        0.068  -0.038 -0.068 -0.007 0.078   0.005 0.018
lipophilicity     TotEcc       -0.018  -0.046 -0.076 -0.016 0.067   0.015 0.048
         bbbp  InfoH_deg        0.006   0.040 -0.003  0.084 0.071   0.010 0.035
         bbbp InfoH_spec        0.014  -0.004 -0.048  0.039 0.107   0.005 0.018
         bbbp Bonchev_Id        0.038  -0.025 -0.068  0.019 0.085   0.005 0.018
         bbbp    Btw_sum        0.010  -0.018 -0.061  0.026 0.086   0.005 0.018
         bbbp    Cls_sum        0.004  -0.010 -0.053  0.034 0.109   0.005 0.018
         bbbp    Eig_sum        0.017   0.000 -0.043  0.044 0.107   0.005 0.018
         bbbp  Harm_cent        0.000  -0.029 -0.072  0.014 0.082   0.005 0.018
         bbbp    AlgConn        0.006  -0.001 -0.044  0.043 0.090   0.005 0.018
         bbbp  Rho_x_Dia       -0.006  -0.002 -0.045  0.041 0.097   0.005 0.018
         bbbp     EE_x_W        0.012  -0.026 -0.070  0.017 0.088   0.005 0.018
         bbbp    FourCyc       -0.220  -0.171 -0.212 -0.128 0.218   0.005 0.018
         bbbp   CutVerts        0.002  -0.015 -0.058  0.029 0.102   0.005 0.018
         bbbp     TotEcc        0.035  -0.019 -0.063  0.024 0.120   0.005 0.018
         bbbp     Radius       -0.012  -0.036 -0.079  0.007 0.084   0.005 0.018
         bbbp    MeanEcc        0.001  -0.012 -0.055  0.032 0.087   0.005 0.018
```

## Linear vs non-linear screen agreement (|.| >= 0.10)

- both: 2; linear only: 7; non-linear only: 1; neither: 70; BH q < 0.05 (dCor perm): 25 of 80.

## Calibration warning (added after the run)

The in-span negative control REF_invNirmala is rejected by the dCor permutation test on
Lipophilicity (p = 0.005) and BBBP (p = 0.010), with rho_nl CI excluding 0 on Lipophilicity
(-0.044 [-0.074, -0.014]). The pure-noise control is not rejected anywhere. So the dCor
permutation test on cross-fitted RF residuals is anti-conservative for X-dependent z: RF
nuisance error in rz and ry are both functions of X, and permuting rz against ry breaks that
shared X-dependence. The 25/80 BH "hits" therefore must NOT be read as target information.
An effect-size floor of about |rho_nl| ~ 0.05 (and dCor ~ 0.08-0.10) is set by the in-span
control; only effects clearly above that floor are interpretable.
