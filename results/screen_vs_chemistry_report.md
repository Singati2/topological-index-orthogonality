# Target-aware screen against a chemistry-aware baseline

                                              0           1              2               3
dataset                                    esol    freesolv  lipophilicity            bbbp
task                                 regression  regression     regression  classification
n                                          1127         639           4200            2039
n_rdkit_cols                                161         149            158             162
n_both_cols                                 191         179            188             192
n_topo_dropped_collinear                      0           0              0               0
n_rdkit_dropped_collinear_with_topo           0           0              0               0
r2_rdkit                               0.929264    0.971889       0.591484        0.544732
r2_rdkit_topo                          0.939751    0.975553       0.606514        0.558965
r2_rdkit_topo_cands                    0.944433    0.977222       0.616637        0.566434
r2_topo                                0.489407    0.246408       0.145626        0.314388
cv_rdkit                                0.72879    1.050502       0.829467        0.891282
cv_rdkit_topo                          0.730811    1.043591       0.824318        0.891816
cv_rdkit_topo_cands                    0.704663    1.043112       0.819709        0.893801
cv_topo                                1.541339    3.493388        1.12298        0.801026
cand_pass_both                                0           0              0               0
cand_inconclusive_both                       13          18              3               1
cand_negligible_both                          6           1             16              18
cand_inspan_both                              1           1              1               1
topo_pass_given_rdkit                         1           0              0               0
rdkit_pass_given_topo                        72          64             25              18

## Candidates not negligible given X_both

      dataset  candidate  max_abs_r_both  pcor_topo  ci_lo_topo  ci_hi_topo  relres_topo verdict_topo  pcor_rdkit  ci_lo_rdkit  ci_hi_rdkit  relres_rdkit verdict_rdkit  pcor_both  ci_lo_both  ci_hi_both  relres_both verdict_both
         esol InfoH_dist        0.886888  -0.049020   -0.107682    0.009982     0.195688 INCONCLUSIVE   -0.254592    -0.312294    -0.195017      0.170621          PASS  -0.081018   -0.143587   -0.017804     0.098339 INCONCLUSIVE
         esol InfoH_spec        0.943353  -0.008813   -0.067751    0.050186     0.110200   NEGLIGIBLE   -0.170960    -0.231169    -0.109446      0.099155          PASS  -0.062786   -0.125592    0.000521     0.056775 INCONCLUSIVE
         esol  InfoH_ecc        0.840725  -0.038994   -0.097743    0.020025     0.348839   NEGLIGIBLE   -0.250371    -0.308221    -0.190677      0.315404          PASS  -0.133024   -0.194689   -0.070311     0.253695 INCONCLUSIVE
         esol Bonchev_Id        0.991588  -0.113053   -0.170887   -0.054443     0.011734 INCONCLUSIVE    0.069217     0.006559     0.131334      0.344315  INCONCLUSIVE   0.059781   -0.003537    0.122622     0.007731 INCONCLUSIVE
         esol    Btw_sum        0.918439  -0.072980   -0.131388   -0.014068     0.094481 INCONCLUSIVE   -0.171205    -0.231408    -0.109695      0.109734          PASS  -0.054268   -0.117171    0.009068     0.041521 INCONCLUSIVE
         esol    Cls_sum        0.938358  -0.001524   -0.060492    0.057454     0.087193   NEGLIGIBLE    0.135279     0.073213     0.196301      0.118147  INCONCLUSIVE   0.067686    0.004400    0.130432     0.048450 INCONCLUSIVE
         esol    Eig_sum        0.917407   0.001172   -0.057805    0.060141     0.148404   NEGLIGIBLE    0.021503    -0.041239     0.084076      0.152224    NEGLIGIBLE   0.037218   -0.026149    0.100286     0.091568 INCONCLUSIVE
         esol  Rho_x_Dia        0.962778  -0.076849   -0.135209   -0.017957     0.134030 INCONCLUSIVE   -0.169072    -0.229328    -0.107525      0.143132          PASS  -0.057363   -0.120232    0.005963     0.092803 INCONCLUSIVE
         esol     EE_x_W        0.979979  -0.055546   -0.114145    0.003438     0.028698 INCONCLUSIVE   -0.234215    -0.292606    -0.174084      0.058508          PASS  -0.112010   -0.174081   -0.049053     0.016387 INCONCLUSIVE
         esol  MeanClust        0.757765   0.038683   -0.020336    0.097434     0.780136   NEGLIGIBLE   -0.053271    -0.115571     0.009448      0.574933  INCONCLUSIVE  -0.047812   -0.110781    0.015540     0.506385 INCONCLUSIVE
         esol     TotEcc        0.957318  -0.072526   -0.130939   -0.013612     0.055762 INCONCLUSIVE    0.069913     0.007258     0.132021      0.157662  INCONCLUSIVE   0.074929    0.011679    0.137581     0.035017 INCONCLUSIVE
         esol     Radius        0.914312  -0.070126   -0.128568   -0.011200     0.199944 INCONCLUSIVE   -0.014507    -0.077123     0.048224      0.200743    NEGLIGIBLE   0.065169    0.001872    0.127946     0.161591 INCONCLUSIVE
         esol    MeanEcc        0.929576  -0.109008   -0.166908   -0.050358     0.117631 INCONCLUSIVE   -0.167490    -0.227785    -0.105915      0.134440          PASS  -0.042557   -0.105578    0.020804     0.068268 INCONCLUSIVE
     freesolv  InfoH_deg        0.649173  -0.061361   -0.139618    0.017658     0.469868 INCONCLUSIVE   -0.096875    -0.182502    -0.009789      0.440998  INCONCLUSIVE   0.017194   -0.071679    0.105795     0.331673 INCONCLUSIVE
     freesolv InfoH_spec        0.949212   0.079738    0.000810    0.157678     0.110438 INCONCLUSIVE   -0.062665    -0.149019     0.024639      0.096357  INCONCLUSIVE   0.106591    0.017998    0.193523     0.054428 INCONCLUSIVE
     freesolv  InfoH_ecc        0.758580   0.026329   -0.052714    0.105043     0.334748 INCONCLUSIVE   -0.163263    -0.246917    -0.077193      0.311758  INCONCLUSIVE  -0.074083   -0.161782    0.014778     0.243220 INCONCLUSIVE
     freesolv Bonchev_Id        0.998574  -0.039489   -0.118054    0.039567     0.023205 INCONCLUSIVE    0.074453    -0.012799     0.160579      0.136276  INCONCLUSIVE  -0.101523   -0.188586   -0.012875     0.014159 INCONCLUSIVE
     freesolv    Btw_sum        0.879048   0.088339    0.009472    0.166114     0.079809 INCONCLUSIVE   -0.094236    -0.179927    -0.007126      0.095296  INCONCLUSIVE   0.059582   -0.029337    0.147564     0.031806 INCONCLUSIVE
     freesolv    Cls_sum        0.949045  -0.065054   -0.143252    0.013950     0.059575 INCONCLUSIVE    0.065878    -0.021414     0.152173      0.068557  INCONCLUSIVE  -0.038151   -0.126486    0.050784     0.030280 INCONCLUSIVE
     freesolv    Eig_sum        0.922335   0.059234   -0.019791    0.137525     0.064804 INCONCLUSIVE    0.015626    -0.071641     0.102655      0.129197  INCONCLUSIVE   0.030999   -0.057924    0.119433     0.036035 INCONCLUSIVE
     freesolv    Btw_var        0.576320  -0.003067   -0.081980    0.075885     0.377937   NEGLIGIBLE    0.011556    -0.075689     0.098626      0.402021    NEGLIGIBLE  -0.069246   -0.157044    0.019637     0.179429 INCONCLUSIVE
     freesolv    AlgConn        0.720049  -0.094342   -0.171994   -0.015524     0.317338 INCONCLUSIVE    0.063888    -0.023411     0.150220      0.268433  INCONCLUSIVE  -0.016345   -0.104956    0.072523     0.184589 INCONCLUSIVE
     freesolv  Rho_x_Dia        0.955309  -0.001798   -0.080720    0.077146     0.143037   NEGLIGIBLE   -0.213510    -0.295185    -0.128738      0.138792          PASS  -0.099217   -0.186339   -0.010547     0.097768 INCONCLUSIVE
     freesolv     EE_x_W        0.985111   0.129118    0.050701    0.205952     0.019321 INCONCLUSIVE   -0.151269    -0.235335    -0.064957      0.042395  INCONCLUSIVE   0.046991   -0.041947    0.135191     0.009709 INCONCLUSIVE
     freesolv  Triangles        0.495461   0.106995    0.028301    0.184371     0.145816 INCONCLUSIVE   -0.024284    -0.111217     0.063018      0.511749  INCONCLUSIVE  -0.021444   -0.109997    0.067448     0.102642 INCONCLUSIVE
     freesolv  MeanClust        0.287956   0.051448   -0.027597    0.129854     0.615283 INCONCLUSIVE   -0.036856    -0.123627     0.050475      0.589594  INCONCLUSIVE   0.014105   -0.074752    0.102739     0.399635 INCONCLUSIVE
     freesolv    FourCyc        0.240893  -0.074049   -0.152093    0.004913     0.727829 INCONCLUSIVE   -0.001286    -0.088444     0.085893      0.704707    NEGLIGIBLE  -0.029726   -0.118178    0.059193     0.555746 INCONCLUSIVE
     freesolv   CutVerts        0.875777   0.044098   -0.034957    0.122604     0.177049 INCONCLUSIVE   -0.068380    -0.154627     0.018902      0.080385  INCONCLUSIVE  -0.015568   -0.104188    0.073296     0.061058 INCONCLUSIVE
     freesolv     TotEcc        0.979932  -0.017684   -0.096483    0.061335     0.069321   NEGLIGIBLE   -0.105636    -0.191046    -0.018639      0.068689  INCONCLUSIVE  -0.085112   -0.172572    0.003679     0.036743 INCONCLUSIVE
     freesolv     Radius        0.878064   0.020630   -0.058398    0.099402     0.256440   NEGLIGIBLE   -0.041428    -0.128134     0.045906      0.220916  INCONCLUSIVE   0.020841   -0.068048    0.109401     0.190813 INCONCLUSIVE
     freesolv    MeanEcc        0.912328   0.061272   -0.017747    0.139530     0.114579 INCONCLUSIVE   -0.172912    -0.256219    -0.087056      0.114931  INCONCLUSIVE  -0.025661   -0.114164    0.063246     0.061802 INCONCLUSIVE
lipophilicity InfoH_dist        0.858278   0.041993    0.011688    0.072222     0.165051   NEGLIGIBLE    0.047948     0.017190     0.078617      0.300095    NEGLIGIBLE   0.075718    0.044970    0.106323     0.128547 INCONCLUSIVE
lipophilicity    Cls_sum        0.755826  -0.026238   -0.056514    0.004085     0.181226   NEGLIGIBLE   -0.060425    -0.091040    -0.029696      0.405094    NEGLIGIBLE  -0.094925   -0.125411   -0.064260     0.137180 INCONCLUSIVE
lipophilicity     EE_x_W        0.981655   0.068324    0.038082    0.098440     0.032966   NEGLIGIBLE    0.004622    -0.026166     0.035401      0.129600    NEGLIGIBLE   0.070764    0.039998    0.101395     0.025858 INCONCLUSIVE
         bbbp Bonchev_Id        0.996609   0.038097   -0.005558    0.081607     0.013022   NEGLIGIBLE    0.102870     0.057991     0.147334      0.391769  INCONCLUSIVE   0.057870    0.012525    0.102977     0.011119 INCONCLUSIVE

## Topological indices with residual signal beyond RDKit 2D (PASS)

dataset  index  pcor_given_rdkit    ci_lo    ci_hi  relres_on_rdkit verdict
   esol Mostar          0.185234 0.123986 0.245075         0.153343    PASS
