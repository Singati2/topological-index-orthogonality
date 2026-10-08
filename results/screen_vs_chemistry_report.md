# Target-aware screen against a chemistry-aware baseline

                                              0           1              2               3
dataset                                    esol    freesolv  lipophilicity            bbbp
task                                 regression  regression     regression  classification
n                                          1127         639           4200            2039
n_rdkit_cols                                182         164            176             195
n_both_cols                                 201         182            194             214
n_topo_dropped_collinear                      8           8              9               8
n_rdkit_dropped_collinear_with_topo           3           4              3               3
r2_rdkit                               0.935621    0.976624       0.624779        0.576991
r2_rdkit_topo                          0.944775    0.979923       0.634301        0.586534
r2_rdkit_topo_cands                    0.948717    0.981208       0.641195        0.592645
r2_topo                                0.489407    0.246408       0.145626        0.314388
cv_rdkit                               0.695282    1.042058       0.811769        0.889988
cv_rdkit_topo                          0.699126    1.032332       0.812655        0.892021
cv_rdkit_topo_cands                    0.688302    1.033066       0.814165        0.895276
cv_topo                                1.541339    3.493388        1.12298        0.801026
cand_pass_both                                0           0              0               0
cand_inconclusive_both                       11          15              0               0
cand_negligible_both                          8           4             19              19
cand_inspan_both                              1           1              1               1
topo_pass_given_rdkit                         0           1              0               0
rdkit_pass_given_topo                        81          74             22              20

## Candidates not negligible given X_both

 dataset  candidate  max_abs_r_both  pcor_topo  ci_lo_topo  ci_hi_topo  relres_topo verdict_topo  pcor_rdkit  ci_lo_rdkit  ci_hi_rdkit  relres_rdkit verdict_rdkit  pcor_both  ci_lo_both  ci_hi_both  relres_both verdict_both
    esol InfoH_dist        0.886888  -0.049020   -0.107682    0.009982     0.195688 INCONCLUSIVE   -0.197064    -0.257599    -0.134988      0.130700          PASS  -0.086826   -0.150408   -0.022528     0.080536 INCONCLUSIVE
    esol InfoH_spec        0.952455  -0.008813   -0.067751    0.050186     0.110200   NEGLIGIBLE   -0.192712    -0.253371    -0.130544      0.059510          PASS  -0.073819   -0.137588   -0.009440     0.042625 INCONCLUSIVE
    esol  InfoH_ecc        0.840725  -0.038994   -0.097743    0.020025     0.348839   NEGLIGIBLE   -0.208265    -0.268471    -0.146437      0.274557          PASS  -0.134954   -0.197659   -0.071148     0.235063 INCONCLUSIVE
    esol Bonchev_Id        0.991588  -0.113053   -0.170887   -0.054443     0.011734 INCONCLUSIVE   -0.013097    -0.076805     0.050718      0.258101    NEGLIGIBLE   0.107385    0.043261    0.170628     0.007062 INCONCLUSIVE
    esol    Btw_sum        0.918439  -0.072980   -0.131388   -0.014068     0.094481 INCONCLUSIVE   -0.136952    -0.198986    -0.073824      0.084589  INCONCLUSIVE  -0.040639   -0.104789    0.023847     0.034145 INCONCLUSIVE
    esol    Eig_sum        0.918529   0.001172   -0.057805    0.060141     0.148404   NEGLIGIBLE    0.024771    -0.039063     0.088404      0.139019    NEGLIGIBLE   0.056349   -0.008104    0.120336     0.086637 INCONCLUSIVE
    esol    AlgConn        0.704928   0.036894   -0.022127    0.095659     0.402947   NEGLIGIBLE    0.104516     0.041017     0.167174      0.271360  INCONCLUSIVE   0.047167   -0.017309    0.111252     0.213788 INCONCLUSIVE
    esol  Rho_x_Dia        0.962778  -0.076849   -0.135209   -0.017957     0.134030 INCONCLUSIVE   -0.139674    -0.201650    -0.076584      0.120752  INCONCLUSIVE  -0.054910   -0.118913    0.009547     0.086746 INCONCLUSIVE
    esol     EE_x_W        0.979979  -0.055546   -0.114145    0.003438     0.028698 INCONCLUSIVE   -0.236300    -0.295617    -0.175167      0.039171          PASS  -0.101265   -0.164615   -0.037083     0.013133 INCONCLUSIVE
    esol     TotEcc        0.957318  -0.072526   -0.130939   -0.013612     0.055762 INCONCLUSIVE    0.019704    -0.044124     0.083372      0.123414    NEGLIGIBLE   0.109997    0.045899    0.173193     0.032098 INCONCLUSIVE
    esol     Radius        0.914312  -0.070126   -0.128568   -0.011200     0.199944 INCONCLUSIVE    0.015350    -0.048470     0.079045      0.187314    NEGLIGIBLE   0.056263   -0.008190    0.120251     0.153823 INCONCLUSIVE
freesolv InfoH_dist        0.846008   0.086287    0.007404    0.164102     0.150718 INCONCLUSIVE   -0.122652    -0.210302    -0.033046      0.095849  INCONCLUSIVE   0.020838   -0.071025    0.112350     0.055597 INCONCLUSIVE
freesolv InfoH_spec        0.959491   0.079738    0.000810    0.157678     0.110438 INCONCLUSIVE   -0.069254    -0.158238     0.020847      0.049447  INCONCLUSIVE   0.042922   -0.048998    0.134121     0.031968 INCONCLUSIVE
freesolv Bonchev_Id        0.998574  -0.039489   -0.118054    0.039567     0.023205 INCONCLUSIVE    0.037306    -0.052842     0.126851      0.108179  INCONCLUSIVE  -0.057027   -0.147980    0.034882     0.011926 INCONCLUSIVE
freesolv    Btw_sum        0.879048   0.088339    0.009472    0.166114     0.079809 INCONCLUSIVE   -0.086566    -0.175172     0.003432      0.075369  INCONCLUSIVE   0.120463    0.029057    0.209871     0.024180 INCONCLUSIVE
freesolv    Cls_sum        0.946927  -0.065054   -0.143252    0.013950     0.059575 INCONCLUSIVE    0.100541     0.010667     0.188804      0.041122  INCONCLUSIVE  -0.083801   -0.174190    0.007987     0.024602 INCONCLUSIVE
freesolv    Eig_sum        0.914768   0.059234   -0.019791    0.137525     0.064804 INCONCLUSIVE    0.042932    -0.047221     0.132391      0.110093  INCONCLUSIVE   0.060670   -0.031231    0.151554     0.030139 INCONCLUSIVE
freesolv    Btw_var        0.726891  -0.003067   -0.081980    0.075885     0.377937   NEGLIGIBLE    0.018479    -0.071611     0.108270      0.344120  INCONCLUSIVE  -0.015723   -0.107296    0.076114     0.149489 INCONCLUSIVE
freesolv  Rho_x_Dia        0.955309  -0.001798   -0.080720    0.077146     0.143037   NEGLIGIBLE   -0.179627    -0.265310    -0.091129      0.121023  INCONCLUSIVE  -0.037411   -0.128696    0.054503     0.088238 INCONCLUSIVE
freesolv     EE_x_W        0.985111   0.129118    0.050701    0.205952     0.019321 INCONCLUSIVE   -0.163146    -0.249455    -0.074266      0.029905  INCONCLUSIVE   0.078693   -0.013128    0.169199     0.008041 INCONCLUSIVE
freesolv  Triangles        0.495461   0.106995    0.028301    0.184371     0.145816 INCONCLUSIVE   -0.052587    -0.141886     0.037562      0.428216  INCONCLUSIVE  -0.048988   -0.140086    0.042932     0.059545 INCONCLUSIVE
freesolv  MeanClust        0.287956   0.051448   -0.027597    0.129854     0.615283 INCONCLUSIVE   -0.004122    -0.094058     0.085881      0.513715    NEGLIGIBLE   0.039409   -0.052508    0.130664     0.322918 INCONCLUSIVE
freesolv    FourCyc        0.561472  -0.074049   -0.152093    0.004913     0.727829 INCONCLUSIVE   -0.027724    -0.117402     0.062402      0.393882  INCONCLUSIVE   0.018963   -0.072891    0.110498     0.304131 INCONCLUSIVE
freesolv   CutVerts        0.875777   0.044098   -0.034957    0.122604     0.177049 INCONCLUSIVE   -0.055963    -0.145202     0.034180      0.059946  INCONCLUSIVE  -0.015605   -0.107178    0.076232     0.044012 INCONCLUSIVE
freesolv     TotEcc        0.979932  -0.017684   -0.096483    0.061335     0.069321   NEGLIGIBLE   -0.132324    -0.219680    -0.042864      0.053325  INCONCLUSIVE  -0.027040   -0.118473    0.064848     0.030960 INCONCLUSIVE
freesolv    MeanEcc        0.912328   0.061272   -0.017747    0.139530     0.114579 INCONCLUSIVE   -0.169456    -0.255531    -0.080716      0.091595  INCONCLUSIVE   0.022350   -0.069520    0.113844     0.054629 INCONCLUSIVE

## Topological indices with residual signal beyond RDKit 2D (PASS)

 dataset index  pcor_given_rdkit    ci_lo    ci_hi  relres_on_rdkit verdict
freesolv    HR          0.192974 0.104823 0.278116         0.026712    PASS
