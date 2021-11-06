
R CMD BATCH --no-save --no-restore '--args ncase=13 nctrl=10 ncell=360 r_mean=1.2 r_var=1.5' step2_evaluate_methods_scDD.R step2_evaluate_methods_scDD_13_10_360_1.2_1.5.Rout &
R CMD BATCH --no-save --no-restore '--args ncase=13 nctrl=10 ncell=1080 r_mean=1.2 r_var=1.5' step2_evaluate_methods_scDD.R step2_evaluate_methods_scDD_13_10_1080_1.2_1.5.Rout &

R CMD BATCH --no-save --no-restore '--args ncase=5 nctrl=5 ncell=360 r_mean=1.2 r_var=1.5' step2_evaluate_methods_scDD.R step2_evaluate_methods_scDD_5_5_360_1.2_1.5.Rout &
R CMD BATCH --no-save --no-restore '--args ncase=5 nctrl=5 ncell=1080 r_mean=1.2 r_var=1.5' step2_evaluate_methods_scDD.R step2_evaluate_methods_scDD_5_5_1080_1.2_1.5.Rout &


R CMD BATCH --no-save --no-restore '--args ncase=13 nctrl=10 ncell=0 r_mean=1.2 r_var=1.5' step2_evaluate_methods_scDD.R step2_evaluate_methods_scDD_13_10_unequal_n_cell_1.2_1.5.Rout &


R CMD BATCH --no-save --no-restore '--args ncase=20 nctrl=20 ncell=360 r_mean=1.2 r_var=1.5' step2_evaluate_methods_scDD.R step2_evaluate_methods_scDD_20_20_360_1.2_1.5.Rout &
R CMD BATCH --no-save --no-restore '--args ncase=20 nctrl=20 ncell=1080 r_mean=1.2 r_var=1.5' step2_evaluate_methods_scDD.R step2_evaluate_methods_scDD_20_20_1080_1.2_1.5.Rout &

R CMD BATCH --no-save --no-restore '--args ncase=10 nctrl=10 ncell=360 r_mean=1.2 r_var=1.5' step2_evaluate_methods_scDD.R step2_evaluate_methods_scDD_10_10_360_1.2_1.5.Rout &
