library(data.table)
library(ggplot2)

# ==========================================================
# NOTE: Create a csv file with individual subject model fits 
#       for the Bayes factor analysis in Jasp 
# ==========================================================

rm( list = ls() )

dir_name <- 'D:/Studium/Auslandsstudium/TuitionWaver_Master/Masterthesis/Analysis/Exp_Variance/Modelling/sorted/fits'

f_names <- list.files(paste(dir_name, sep=''),
           pattern=paste('fit_2state_bcee*'),
           full.names=TRUE)
f_names <- f_names[1:48]

ldf <- lapply(seq_along(f_names), function(j) {
  z<-fread(f_names[j]);
  z<-t(z);
  z<-data.table(z);
  setnames(z, c('alpha_s', 'beta_s', 'g_sigma_s', 'alpha_f', 'beta_f', 'g_sigma_f'));
})
ldf <- lapply(seq_along(ldf), function(z) ldf[[z]][, sub := rep(z, .N)])
d <- rbindlist(ldf)

## NOTE: add rot_dir indicator column
subs_cnd <- c(1, 2, 2, 1, 2, 2, 0, 2, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 2,
              0, 1, 0, 0, 0, 1, 2, 2, 2, 0, 2, 1, 0, 0, 1, 0, 0, 2, 1, 2, 1, 2,
              2, 2, 2, 1)
condition <- rep(subs_cnd, each=d[, .N, .(sub)][, unique(N)])
d[, cnd := condition]

## NOTE: add rot_dir indicator column
subs_ccw = c(1, 3, 4, 6, 7, 10, 12, 14, 15, 16, 17, 19, 21, 24, 29, 31, 33, 34,
             35, 36, 37, 41, 43, 47)
d[, rot_dir := 'cw']
d[sub %in% subs_ccw, rot_dir := 'ccw']

fwrite(d, '../bayesfactor/fits_subs.csv')
