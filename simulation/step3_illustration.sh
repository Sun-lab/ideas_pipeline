
R CMD BATCH --no-save --no-restore '--args ncase=13 nctrl=10 ncell=360 r_mean=1.2 r_var=1.5' step3_illustration.R step3_illustration_Rout/step3_illustration_13_10_360_1.2_1.5.Rout &

R CMD BATCH --no-save --no-restore '--args ncase=13 nctrl=10 ncell=0 r_mean=1.2 r_var=1.5' step3_illustration.R step3_illustration_Rout/step3_illustration_13_10_unequal_n_cell_1.2_1.5.Rout &

R CMD BATCH --no-save --no-restore '--args ncase=13 nctrl=10 ncell=120 r_mean=1.2 r_var=1.5' step3_illustration.R step3_illustration_Rout/step3_illustration_13_10_120_1.2_1.5.Rout &



R CMD BATCH --no-save --no-restore '--args ncase=5 nctrl=5 ncell=360 r_mean=1.2 r_var=1.5' step3_illustration.R step3_illustration_Rout/step3_illustration_5_5_360_1.2_1.5.Rout &

R CMD BATCH --no-save --no-restore '--args ncase=5 nctrl=5 ncell=120 r_mean=1.2 r_var=1.5' step3_illustration.R step3_illustration_Rout/step3_illustration_5_5_120_1.2_1.5.Rout &

R CMD BATCH --no-save --no-restore '--args ncase=5 nctrl=5 ncell=1080 r_mean=1.2 r_var=1.5' step3_illustration.R step3_illustration_Rout/step3_illustration_5_5_1080_1.2_1.5.Rout &


R CMD BATCH --no-save --no-restore '--args ncase=20 nctrl=20 ncell=360 r_mean=1.2 r_var=1.5' step3_illustration.R step3_illustration_Rout/step3_illustration_20_20_360_1.2_1.5.Rout &

R CMD BATCH --no-save --no-restore '--args ncase=20 nctrl=20 ncell=120 r_mean=1.2 r_var=1.5' step3_illustration.R step3_illustration_Rout/step3_illustration_20_20_120_1.2_1.5.Rout &


