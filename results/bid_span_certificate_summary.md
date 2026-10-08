Data-free certificate: rank of the 18 x 10 f-matrix of the BID baseline over P4 = 10 (sigma_min/sigma_max = 1.20e-04); hence every BID index is an exact linear combination of the baseline on every graph with max degree <= 4 and no isolated vertex.

# BID span certificate

      dataset    n  max_degree  degree_pairs_realised  rank_mij  rank_BID18  rank_BID18_plus_mij  max_rel_resid_mij_on_BID18  certificate_holds greedy_minimal_certifying_subset  size_minimal_subset
         esol 1127           4                     10        10          10                   10                4.522735e-13               True    M1 M2 mM1 mM2 F R SCI H GA AG                   10
     freesolv  639           4                     10        10          10                   10                2.071692e-12               True    M1 M2 mM1 mM2 F R SCI H GA AG                   10
lipophilicity 4200           4                      9         9           9                    9                2.585556e-12               True       M1 M2 mM1 mM2 F R SCI H GA                    9
         bbbp 2039           4                     10        10          10                   10                1.099688e-12               True    M1 M2 mM1 mM2 F R SCI H GA AG                   10

## Conventional pairwise screen vs certificate

      dataset          family   n  pairwise_pass  in_span  max_rel_resid  max_abs_naive_pcor
         bbbp     Loyola grid 100              0      100   1.281607e-12            0.002195
         bbbp non-BID control   2              1        0   7.526229e-02            0.068185
         bbbp   published BID  18              0       18   1.805983e-13            0.001851
         esol     Loyola grid 100              0      100   6.915227e-13            0.028225
         esol non-BID control   2              0        0   5.663453e-02            0.074655
         esol   published BID  18              0       18   1.890472e-13            0.000384
     freesolv     Loyola grid 100              1      100   1.254205e-13            0.036177
     freesolv non-BID control   2              0        0   6.779820e-02            0.068281
     freesolv   published BID  18              0       18   5.552912e-14            0.003636
lipophilicity     Loyola grid 100              0      100   1.929710e-12            0.009328
lipophilicity non-BID control   2              1        0   7.069375e-02            0.022359
lipophilicity   published BID  18              0       18   3.540437e-13            0.001663

## Published variants passing the pairwise screen (all are exactly in span)

      dataset                                   index          family  max_abs_r_baseline30 most_correlated  rel_resid_on_baseline30  in_span  pcor_certified
         esol                   mSO (modified Sombor)   published BID              0.999440               H             5.574251e-14     True        0.000000
         esol                   ESO (elliptic Sombor)   published BID              0.999716               F             6.326795e-15     True        0.000000
         esol                       EU (Euler Sombor)   published BID              0.999957              M1             2.994792e-14     True        0.000000
         esol                 DSO (diminished Sombor)   published BID              0.999920              AG             1.428853e-14     True        0.000000
         esol                 HSO (hyperbolic Sombor)   published BID              0.987333          SO_red             1.890472e-13     True        0.000000
         esol                  SO_p, p=1/2 (p-Sombor)   published BID              0.999871              M1             1.579413e-14     True        0.000000
         esol                    SO_p, p=3 (p-Sombor)   published BID              0.999850              SO             1.386631e-14     True        0.000000
         esol                       KG-SO (KG-Sombor)   published BID              0.999461              SO             2.821153e-14     True        0.000000
         esol                  ASO (augmented Sombor)   published BID              0.999851              AG             2.698970e-14     True        0.000000
         esol CoRSO_60 (cosine-rule Sombor, theta=60)   published BID              0.999169              SO             1.118971e-14     True        0.000000
         esol   EU_lambda=1/2 (variable Euler-Sombor)   published BID              0.999952              SO             2.949986e-14     True        0.000000
         esol                 ISI (inverse sum indeg)   published BID              0.998108              M1             4.253057e-14     True        0.000000
         esol                SDD (symm. division deg)   published BID              0.993356              SO             5.261682e-14     True        0.000000
         esol                             N (Nirmala)   published BID              0.999800             ABS             2.584465e-14     True        0.000000
         esol             IN1 (first inverse Nirmala)   published BID              0.999071              EE             1.621969e-14     True        0.000000
         esol            IN2 (second inverse Nirmala)   published BID              0.999521             ABS             4.016459e-14     True        0.000000
         esol                     GG1 (first Gourava)   published BID              0.999165              M2             1.688575e-14     True        0.000000
         esol                       HM (hyper-Zagreb)   published BID              0.998948               F             6.469995e-15     True        0.000000
         esol                              SO coindex non-BID control              0.998974              PI             7.143747e-03    False        0.039918
         esol               SO_ecc (eccentric Sombor) non-BID control              0.970924              PI             5.663453e-02    False       -0.074655
     freesolv                   mSO (modified Sombor)   published BID              0.999004               H             4.383112e-14     True        0.000000
     freesolv                   ESO (elliptic Sombor)   published BID              0.999403               F             4.101119e-14     True        0.000000
     freesolv                       EU (Euler Sombor)   published BID              0.999916              M1             4.497165e-14     True        0.000000
     freesolv                 DSO (diminished Sombor)   published BID              0.999830              AG             4.567699e-14     True        0.000000
     freesolv                 HSO (hyperbolic Sombor)   published BID              0.982834          SO_red             5.552912e-14     True        0.000000
     freesolv                  SO_p, p=1/2 (p-Sombor)   published BID              0.999745              M1             4.703327e-14     True        0.000000
     freesolv                    SO_p, p=3 (p-Sombor)   published BID              0.999728              SO             4.602774e-14     True        0.000000
     freesolv                       KG-SO (KG-Sombor)   published BID              0.999424              SO             4.943197e-14     True        0.000000
     freesolv                  ASO (augmented Sombor)   published BID              0.999498              AG             4.436818e-14     True        0.000000
     freesolv CoRSO_60 (cosine-rule Sombor, theta=60)   published BID              0.998448              SO             4.277072e-14     True        0.000000
     freesolv   EU_lambda=1/2 (variable Euler-Sombor)   published BID              0.999907              SO             4.527406e-14     True        0.000000
     freesolv                 ISI (inverse sum indeg)   published BID              0.997062             ABS             5.078807e-14     True        0.000000
     freesolv                SDD (symm. division deg)   published BID              0.991531          SO_red             5.518834e-14     True        0.000000
     freesolv                             N (Nirmala)   published BID              0.999635             ABS             4.453479e-14     True        0.000000
     freesolv             IN1 (first inverse Nirmala)   published BID              0.999001              EE             4.917605e-14     True        0.000000
     freesolv            IN2 (second inverse Nirmala)   published BID              0.998806             ABS             3.990634e-14     True        0.000000
     freesolv                     GG1 (first Gourava)   published BID              0.999174              M2             4.023749e-14     True        0.000000
     freesolv                       HM (hyper-Zagreb)   published BID              0.997808               F             4.078929e-14     True        0.000000
     freesolv                              SO coindex non-BID control              0.997741              PI             6.487413e-03    False        0.068281
     freesolv               SO_ecc (eccentric Sombor) non-BID control              0.975785               R             6.779820e-02    False       -0.015508
lipophilicity                   mSO (modified Sombor)   published BID              0.999468               H             8.574418e-14     True        0.000000
lipophilicity                   ESO (elliptic Sombor)   published BID              0.999491               F             9.295853e-15     True        0.000000
lipophilicity                       EU (Euler Sombor)   published BID              0.999939              M1             6.781518e-14     True        0.000000
lipophilicity                 DSO (diminished Sombor)   published BID              0.999886              AG             1.318161e-14     True        0.000000
lipophilicity                 HSO (hyperbolic Sombor)   published BID              0.984904          SO_red             3.540437e-13     True        0.000000
lipophilicity                  SO_p, p=1/2 (p-Sombor)   published BID              0.999814              M1             5.931628e-14     True        0.000000
lipophilicity                    SO_p, p=3 (p-Sombor)   published BID              0.999800              SO             2.145090e-14     True        0.000000
lipophilicity                       KG-SO (KG-Sombor)   published BID              0.999482              SO             4.879548e-14     True        0.000000
lipophilicity                  ASO (augmented Sombor)   published BID              0.999863              AG             4.410558e-14     True        0.000000
lipophilicity CoRSO_60 (cosine-rule Sombor, theta=60)   published BID              0.998843              SO             4.974751e-15     True        0.000000
lipophilicity   EU_lambda=1/2 (variable Euler-Sombor)   published BID              0.999932              SO             3.674057e-14     True        0.000000
lipophilicity                 ISI (inverse sum indeg)   published BID              0.997642             ABS             1.132222e-13     True        0.000000
lipophilicity                SDD (symm. division deg)   published BID              0.991805          SO_red             1.125384e-13     True        0.000000
lipophilicity                             N (Nirmala)   published BID              0.999874             ABS             4.723894e-14     True        0.000000
lipophilicity             IN1 (first inverse Nirmala)   published BID              0.999326             ABC             1.296261e-14     True        0.000000
lipophilicity            IN2 (second inverse Nirmala)   published BID              0.999340             ABS             7.823959e-14     True        0.000000
lipophilicity                     GG1 (first Gourava)   published BID              0.999174              M2             3.621650e-14     True        0.000000
lipophilicity                       HM (hyper-Zagreb)   published BID              0.998139               F             1.185377e-14     True        0.000000
lipophilicity                              SO coindex non-BID control              0.997353              PI             7.697771e-03    False        0.022359
lipophilicity               SO_ecc (eccentric Sombor) non-BID control              0.942664              PI             7.069375e-02    False       -0.021578
         bbbp                   mSO (modified Sombor)   published BID              0.999506               H             5.118482e-14     True        0.000000
         bbbp                   ESO (elliptic Sombor)   published BID              0.999748               F             6.254930e-15     True        0.000000
         bbbp                       EU (Euler Sombor)   published BID              0.999962              M1             6.121142e-14     True        0.000000
         bbbp                 DSO (diminished Sombor)   published BID              0.999930              AG             4.799616e-15     True        0.000000
         bbbp                 HSO (hyperbolic Sombor)   published BID              0.990095          SO_red             1.184437e-13     True        0.000000
         bbbp                  SO_p, p=1/2 (p-Sombor)   published BID              0.999885              M1             2.855953e-14     True        0.000000
         bbbp                    SO_p, p=3 (p-Sombor)   published BID              0.999873              SO             7.145404e-15     True        0.000000
         bbbp                       KG-SO (KG-Sombor)   published BID              0.999492              SO             1.122458e-13     True        0.000000
         bbbp                  ASO (augmented Sombor)   published BID              0.999875              AG             2.582534e-14     True        0.000000
         bbbp CoRSO_60 (cosine-rule Sombor, theta=60)   published BID              0.999287              SO             4.282126e-15     True        0.000000
         bbbp   EU_lambda=1/2 (variable Euler-Sombor)   published BID              0.999958              SO             1.303895e-14     True        0.000000
         bbbp                 ISI (inverse sum indeg)   published BID              0.998303              M1             1.233448e-13     True        0.000000
         bbbp                SDD (symm. division deg)   published BID              0.994351              SO             9.182504e-14     True        0.000000
         bbbp                             N (Nirmala)   published BID              0.999803             ABS             2.459725e-14     True        0.000000
         bbbp             IN1 (first inverse Nirmala)   published BID              0.999401             ABC             7.788908e-15     True        0.000000
         bbbp            IN2 (second inverse Nirmala)   published BID              0.999623             ABS             1.805983e-13     True        0.000000
         bbbp                     GG1 (first Gourava)   published BID              0.999159              M2             4.639666e-14     True        0.000000
         bbbp                       HM (hyper-Zagreb)   published BID              0.999069               F             1.032110e-14     True        0.000000
         bbbp                              SO coindex non-BID control              0.998933              PI             7.118233e-03    False       -0.068185
         bbbp               SO_ecc (eccentric Sombor) non-BID control              0.947485             mM2             7.526229e-02    False        0.037907
Empty DataFrame
Columns: [dataset, index, max_abs_r_baseline30, most_correlated, rel_resid_on_baseline30, in_span]
Index: []
