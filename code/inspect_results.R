library(data.table)
library(ggplot2)
library(ggpubr)

rm( list = ls() )

setwd('D:/Studium/Auslandsstudium/TuitionWaver_Master/Masterthesis/Analysis/Exp_Variance/Modelling/sorted/code')
# # ------------------------------------------------------------------------------------------------------------

# # Extracting csv files for SubjData
# ----------------------------------- 
# setwd("../../../MatlabAnalysis/AnalysisData/Exp_Variance_MissingsReplaced_CorrectRTs")
# 
# dir_name <- "../../../Modelling/sorted/data"
# 
# 
# # Extracting csv files for SubjData
# # ---------------------------------
# myfolders <- list.dirs(recursive = FALSE)
# myfolders <- myfolders[1:8]
# 
# for (h in 1:length(myfolders)){
#   subjdatafiles <- list.files(myfolders[h], pattern = "*SubjData.csv")
#   for (g in 1:length(subjdatafiles)){
#     file.copy(paste0(myfolders[h], "/", subjdatafiles[g]), dir_name)
#   }
# }
# 
# rm( list = ls() )
# setwd("../../../Modelling/sorted/code")
# # ------------------------------------------------------------------------------------------------------------

# Transform Dataframe Subject Main Data
# -------------------------------------
col_names <- c("Target", "Trial_Phase", "Appl_Perturb",
               "imv_X", "imv_Y",
               "Endpoint_X", "Endpoint_Y",
               "imv_Error", "Endpoint_Error",
               "imv_Error_Mean", "Endpoint_Error_Mean",
               "MT", "RT", "Max_Vel")

phase_names_f <- c('Prebaseline_SubjData', 'Familiarisation_SubjData',
                 'Baseline_NFB_SubjData', 'Baseline_SubjData',
                 'Postbaseline_SubjData', 'Training_SubjData',
                 'Generalisation_SubjData', 'Relearning_SubjData',
                 'Washout_SubjData')

phase_names <- c('Prebaseline', 'Familiarisation', 'Baseline_NFB', 'Baseline_S',
                 'Postbaseline', 'Training', 'Generalisation', 'Relearning',
                 'Washout')

phase_order <- 1:9

ldf_rec <- list()
for (i in 1:length(phase_names_f)) {
    f_names <- list.files(paste('../data/', sep=''),
                          pattern=paste(phase_names_f[i]),
                          full.names=TRUE)

    ldf <- lapply(f_names, function(z) {z<-fread(z); setnames(z,col_names)})
    ldf <- lapply(seq_along(ldf), function(z) ldf[[z]][, sub := rep(z, .N)])
    ldf <- lapply(seq_along(ldf), function(z) ldf[[z]][, phase := rep(phase_names[i], .N)])
    ldf <- lapply(seq_along(ldf), function(z) ldf[[z]][, phase_order := rep(phase_order[i], .N)])
    ldf_rec <- c(ldf_rec, ldf)
}

d <- rbindlist(ldf_rec)
d <- d[order(sub)]

## NOTE: Adapt the NFB targets (1 to 12)
d[Target %in% 13:24, Target := as.integer(Target - 12)]

## NOTE: add target_deg indicator column
d[, target_deg := -1]
d[Target == 1, target_deg := 0]
d[Target == 2, target_deg := 30]
d[Target == 3, target_deg := 60]
d[Target == 4, target_deg := 90]
d[Target == 5, target_deg := 120]
d[Target == 6, target_deg := 150]
d[Target == 7, target_deg := 180]
d[Target == 8, target_deg := -150]
d[Target == 9, target_deg := -120]
d[Target == 10, target_deg := -90]
d[Target == 11, target_deg := -60]
d[Target == 12, target_deg := -30]

## NOTE: add trial indicator column (fix problem with nfb and fb trials in baseline)
d <- d[order(sub, phase_order, Trial_Phase)]
d[, trial := 1:.N, .(sub)]
d[phase %in% c("Baseline_S", "Baseline_NFB"), trial := Trial_Phase+240]

## NOTE: add aware indicator column
subs_explicit = c(1, 2, 3, 5, 9, 10, 16, 18, 19, 20, 24, 28, 29, 30, 31, 33, 34,
                  40, 41, 44, 45, 46, 47, 48)
d[, aware := 'implicit']
d[sub %in% subs_explicit, aware := 'explicit']

## NOTE: add condition (0: 0°, 1: 4°, 2: 12°) indicator column
subs_rot <- c(1, 2, 2, 1, 2, 2, 0, 2, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 2,
              0, 1, 0, 0, 0, 1, 2, 2, 2, 0, 2, 1, 0, 0, 1, 0, 0, 2, 1, 2, 1, 2,
              2, 2, 2, 1)
condition <- rep(subs_rot, each=d[, .N, .(sub)][, unique(N)])
d[, cnd := condition]

## NOTE: add rot_dir indicator column
subs_ccw = c(1, 3, 4, 6, 7, 10, 12, 14, 15, 16, 17, 19, 21, 24, 29, 31, 33, 34,
             35, 36, 37, 41, 43, 47)
d[, rot_dir := 'cw']
d[sub %in% subs_ccw, rot_dir := 'ccw']

## NOTE: add phase means for baseline correction
d[, ee_mean := mean(Endpoint_Error, na.rm=T), .(sub, phase, Target)]
d[, ee_mean_correction_fb := ee_mean[which(phase=="Baseline_S")][1], .(Target, rot_dir)]
d[, ee_mean_correction_nfb := ee_mean[which(phase=="Baseline_NFB")][1], .(Target, rot_dir)]

## NOTE: flip rot_dir for bcee and plots
d[rot_dir == 'cw' &
    phase %in% c('Training', 'Generalisation', 'Relearning', 'Washout'),
  Endpoint_Error := -1 * Endpoint_Error]

d[rot_dir == 'cw' &
    phase %in% c('Training', 'Generalisation', 'Relearning', 'Washout'),
  imv_Error := -1 * imv_Error]

d[phase %in% c('Prebaseline', 'Baseline_NFB', 'Generalisation', 'Washout'), 
  Appl_Perturb := NaN]
d[rot_dir == 'cw' & 
    phase %in% c('Training', 'Generalisation', 'Relearning', 'Washout'),
  Appl_Perturb := -1 * Appl_Perturb]

## NOTE: Perform baseline correction
d[, bcee := Endpoint_Error]
d[phase %in% c("Generalisation", "Washout"), bcee := Endpoint_Error - ee_mean_correction_nfb]
d[phase %in% c("Postbaseline", "Training", "Relearning"), bcee := Endpoint_Error - ee_mean_correction_nfb]
# --------------------------------------------------------------------------------------------------------

# # Create csv file for python modelling
# # ------------------------------------
# fwrite(d, '../fit_input/master_data.csv')
# # =======================================================================================================

# # Create NA overview tables per subject and block
# # -----------------------------------------------
# ## NOTE: Missings per subject and Block
# d <- d[, phase := factor(phase, levels = c("Prebaseline", "Familiarisation", 
#                                         "Baseline", "Baseline_NFB", "Postbaseline",
#                                         "Training", "Generalisation", 
#                                         "Relearning", "Washout"))]
# df_MissBlock = as.data.frame(matrix(ncol = 0, nrow = 0))  
# 
# nsubs = c(1:48)
# for(i in nsubs) {
#   for(j in levels(d$phase)) {
#     # df_MissBlock[i, j] = sum(is.na(tbl[Block == j & Subj == i, Endpoint_X]))
#     abs <- sum(is.na(d[phase == j & sub == i, Endpoint_X]))
#     perc <- round(abs / max(d$Trial_Phase[d$phase == j &
#                                             d$sub == i], na.rm = TRUE), 2) * 100
#     df_MissBlock[i, j] = paste0(abs, " (", perc, "%)")
#   }
# }
# 
# ## NOTE: Missings per target, block and subject
# df_MissTarget = as.data.frame(matrix(ncol = 0, nrow = 0))
# 
# find_level <- function(phase) {
#   if(phase %in% c(6, 8, 9)) {
#     level <- "0"
#   } else {
#     level <- levels(as.factor(d$target_deg))
#   }
#   return(level)
# }
# 
# colname <- c()
# for(i in 1:48) {
#   m = 1
#   print(paste0("i = ", i))
#   for(j in 1:length(levels(d$phase))) {
#     print(paste0("j = ", j))
#     for(k in 1:length(find_level(j))) {
#       print(paste0("k = ", k))
#       l = m + k - 1
#       df_MissTarget[i, l] = sum(is.na(d[phase == levels(d$phase)[j] &
#                                           sub == i &
#                                           target_deg == find_level(j)[k],
#                                         Endpoint_X]))
#       # col <- paste(levels(d$phase)[j], "(", find_level(j)[k], "°)")
#       col <- paste(find_level(j)[k])
#       colname <- c(colname, col)
#     }
#     m = m + length(find_level(j))
#   }
# }
# colname <- colname[1:75]
# colnames(df_MissTarget) <- colname
# 
# ## NOTE: Missings in Percent and over Target for Generalisation
# subs_idx <- factor(d$sub)
# Sub <- levels(subs_idx)
# df_MissPhase <- cbind(Sub, df_MissBlock[,1:6], df_MissBlock[,8:9])
# df_MissGen <- cbind(Sub, df_MissTarget[, 62:73])
# df_MissTotal <- cbind(Sub, df_MissBlock[,1:6], df_MissTarget[, 62:73], df_MissBlock[,8:9])
# 
# fwrite(df_MissPhase, '../missings/phase_missings_per_sub.csv')
# fwrite(df_MissGen, '../missings/general_missings_per_sub.csv')
# fwrite(df_MissTotal, '../missings/total_missings_per_sub.csv')
# # =======================================================================================================

# Plots
# -----

# Illustration experimental design
# --------------------------------
## NOTE: plot experimental procedure with applied rotation for all groups
dd <- d[target_deg == 0, .(mean(Appl_Perturb, na.rm=T), sd(Appl_Perturb, na.rm=T)/sqrt(16)), .(cnd, trial)]
ggplot(dd, aes(trial, V1, color=as.factor(cnd))) +
  geom_line(alpha=0.4) +
  labs(title = 'Illustration Experimantal Design',
       x = 'Trial',
       y = 'Applied Roation Perturbation (in °)',
       color = 'Variance \n Group') +
  scale_x_continuous(breaks=c(120, 240, 480, 600, 1000, 1072, 1172, 1272), 
                     limits = c(0, 1272), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-30, 50), expand = c(0, 0)) +
  geom_segment(aes(x=120,xend=120,y=-30,yend=15), color = 'darkgrey', linetype = 'dashed') +
  geom_segment(aes(x=240,xend=240,y=-30,yend=15), color = 'darkgrey', linetype = 'dashed') +
  geom_segment(aes(x=480,xend=480,y=-30,yend=15), color = 'darkgrey', linetype = 'dashed') +
  geom_segment(aes(x=600,xend=600,y=-30,yend=15), color = 'darkgrey') +
  geom_segment(aes(x=1000,xend=1000,y=-30,yend=15), color = 'darkgrey') +
  geom_segment(aes(x=1072,xend=1072,y=-30,yend=15), color = 'darkgrey') +
  geom_segment(aes(x=1172,xend=1172,y=-30,yend=15), color = 'darkgrey') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = c(60, 180, 360, 540, 800, 1036, 1122, 1222), y = -25, colour = 'grey', size = 2.5,
           label = c('Prebaseline', 'Familiarisation', 'Baseline', 'Postbaseline', 
                     'Training', 'General.', 'Relearning', 'Washout'))
# ggsave('../figures/illustrate_exp_design/exp_design_applperturb.png', width=10, height=4)


# Results experiment: Illustrate Design Overall performance
# ---------------------------------------------------------
## NOTE: plot all trials collapsed across CW / CCW
dd <- d[, .(mean(bcee, na.rm=T), sd(bcee, na.rm=T)/sqrt(16)), .(cnd, trial)]
ggplot(dd, aes(trial, V1, color=as.factor(cnd))) +
  geom_line(alpha=0.4) +
  labs(title = 'All Trials Collapsed across Variance Condition',
       x = 'Trial', 
       y = 'Baseline corrected EE (in °)',
       color = 'Variance \n Group') +
  scale_x_continuous(breaks=c(120, 240, 480, 600, 1000, 1072, 1172, 1272), limits = c(0, 1272), expand = c(0, 0)) +
  scale_y_continuous(breaks=c(-10, -5, 0, 5, 10, 15, 20), expand = c(0, 0)) +
  geom_segment(aes(x=120,xend=240,y=0,yend=0), color = 'azure3') +
  geom_segment(aes(x=240,xend=480,y=0,yend=0), color = 'azure3', linetype = 'dashed') +
  geom_segment(aes(x=480,xend=600,y=0,yend=0), color = 'azure3') +
  geom_segment(aes(x=600,xend=1000,y=15,yend=15), color = 'azure3') +
  geom_segment(aes(x=1072,xend=1172,y=15,yend=15), color = 'azure3') +
  geom_segment(aes(x=120,xend=120,y=-10,yend=15), color = 'azure3', linetype = 'dashed') +
  geom_segment(aes(x=240,xend=240,y=-10,yend=15), color = 'azure3', linetype = 'dashed') +
  geom_segment(aes(x=480,xend=480,y=-10,yend=15), color = 'azure3', linetype = 'dashed') +
  geom_segment(aes(x=600,xend=600,y=-10,yend=15), color = 'azure3') +
  geom_segment(aes(x=1000,xend=1000,y=-10,yend=15), color = 'azure3') +
  geom_segment(aes(x=1072,xend=1072,y=-10,yend=15), color = 'azure3') +
  geom_segment(aes(x=1172,xend=1172,y=-10,yend=15), color = 'azure3') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = c(60, 180, 360, 540, 800, 1036, 1122, 1222), y = -9, colour = 'grey', size = 2.5,
           label = c('Prebaseline', 'Familiarisation', 'Baseline', 'Postbaseline', 
                     'Training', 'General.', 'Relearning', 'Washout'))
# ggsave('../figures/illustrate_exp_design/exp_design_performance.png', width=10, height=4)


# Results experiment: all trials collapsed across cnd
# -------------------------------------------------
## Note: Performance accross all trials bcee
dd <- d[phase %in% c('Training', 'Generalisation', 'Relearning', 'Washout'),
        .(mean(bcee, na.rm=T), sd(bcee, na.rm=T)/sqrt(16)), .(rot_dir, trial)]

plt1 <- ggplot(dd, aes(trial, V1, color=as.factor(rot_dir))) +
  geom_line(alpha=0.4) +
  labs(title = 'Overall Performance',
       x = 'Trial', 
       y = 'Baseline corrected EE (in °)',
       color = 'Rotation \n Direction') +
  # scale_x_continuous(breaks=seq(0, 1200, 200), limits = c(0, 1272), expand = c(0, 0)) +
  scale_y_continuous(breaks=c(-10, -5, 0, 5, 10, 15, 20), expand = c(0, 0)) +
  geom_segment(aes(x=600,xend=1000,y=15,yend=15), color = 'azure3') +
  geom_segment(aes(x=1072,xend=1172,y=15,yend=15), color = 'azure3') +
  geom_segment(aes(x=600,xend=600,y=0,yend=15), color = 'azure3') +
  geom_segment(aes(x=1000,xend=1000,y=0,yend=15), color = 'azure3') +
  geom_segment(aes(x=1072,xend=1072,y=0,yend=15), color = 'azure3') +
  geom_segment(aes(x=1172,xend=1172,y=0,yend=15), color = 'azure3') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

## NOTE: generalisation curve bcee
dd <- d[phase == "Generalisation",
        .(mean(bcee, na.rm = T), sd(bcee, na.rm = T)/sqrt(.N)),
        .(rot_dir, target_deg)]

plt2 <- ggplot(dd, aes(x = target_deg, y = V1, colour = as.factor(rot_dir))) +
  geom_line(alpha = 0.6) +
  labs(title = 'Generalisation across Targets',
       x = 'Target Position (in °)', 
       y = 'Baseline corrected EE (in °)',
       color = 'Rotation \n Direction') +
  geom_errorbar(aes(ymin=V1-V2, ymax=V1+V2), width=10) +
  scale_x_continuous(breaks=c(-120, -60, 0, 60, 120), 
                     limits = c(-150, 150)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_blank())

ggarrange(plt1, plt2, widths = c(2, 1), ncol = 2, legend = 'right', common.legend = T)
# ggsave('../figures/inspect_difference_rotation/comparison_cw_ccw.png', width=10, height=2)


# Results experiment: reaches against applied rotation
# ----------------------------------------------------
## NOTE: plot applied rotation against reaches seperately for each variance group
dd <- d[phase %in% c('Training', 'Generalisation', 'Relearning', 'Washout')]
dd <- dd[, .(mean(Appl_Perturb, na.rm=T), mean(bcee, na.rm = T), 
             sd(Appl_Perturb, na.rm=T)/sqrt(16)), .(cnd, rot_dir, trial)]
ggplot(dd, aes(trial, V2, color=as.factor(cnd))) +
  geom_line(alpha=0.6) +
  geom_line(aes(trial, V1), alpha=0.4, colour = 'grey') +
  labs(x = 'Trial',
       y = 'Baseline corrected EE (in °)',
       color = 'Variance \n Group') +
  scale_y_continuous(limits = c(-10, 40)) +
  geom_segment(aes(x=600,xend=1000,y=15,yend=15), color = 'azure3') +
  geom_segment(aes(x=1072,xend=1172,y=15,yend=15), color = 'azure3') +
  geom_segment(aes(x=600,xend=600,y=-10,yend=15), color = 'azure3') +
  geom_segment(aes(x=1000,xend=1000,y=-10,yend=15), color = 'azure3') +
  geom_segment(aes(x=1072,xend=1072,y=-10,yend=15), color = 'azure3') +
  geom_segment(aes(x=1172,xend=1172,y=-10,yend=15), color = 'azure3') +
  theme_classic() +
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~rot_dir*cnd, scales = 'free',
             labeller = label_wrap_gen(multi_line=FALSE))
# ggsave('../figures/inspect_overall_performance/rot_against_reaches_pergroup.png',
#        width=10, height=4)


# Results experiment: Generalisation overview
# -------------------------------------------
## NOTE: Plot generalisation function bcee
dd <- d[phase == "Generalisation",
        .(mean(bcee, na.rm = T), sd(bcee, na.rm = T)/sqrt(16)),
        .(cnd, target_deg, rot_dir)]

ggplot(dd, aes(x = target_deg, y = V1, colour = as.factor(cnd))) +
  geom_line(alpha = 0.6) +
  geom_errorbar(aes(ymin=V1-V2, ymax=V1+V2), width=10, alpha = 0.4) +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150, -150, 
                                -120, -90, -60, -30), 
                     limits = c(-150, 150)) +
  scale_y_continuous(sec.axis = sec_axis(~ . / 15 * 100, 
                                         name = "Adaptation (in %)")) +
  facet_wrap(~rot_dir) +
  theme_classic() +
  labs(title = 'Generalisation across Targets',
       y = 'Baseline corrected EE (in °)', 
       x = 'Target Position (in °)',
       colour = 'Variance \n Group') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(strip.background = element_rect(colour = "white", fill = 'white')) +
  geom_hline(yintercept = 0, color='azure3', linetype = 'dashed')
# ggsave('../figures/inspect_overall_performance/generalisation_bcee_percadapt.png',
#        width=10, height=4)


# Methodological Questions: Difference between imv vs. EE
# -------------------------------------------------------
## NOTE: overall performance all phases EE
dd <- d[, .(mean(Endpoint_Error, na.rm=T), 
            sd(Endpoint_Error, na.rm=T)/sqrt(16)), .(cnd, trial, rot_dir)]
plt1 <- ggplot(dd, aes(trial, V1, colour=as.factor(cnd))) +
  geom_line(alpha=0.4) +
  scale_y_continuous(limits = c(-15, 30)) +
  labs(title = 'Overall Performance across all Phases',
       x = 'Trial',
       y = 'EE (in °)',
       color = 'Variance \n Group') +
  scale_x_continuous(breaks=seq(0, 1200, 200), expand = c(0, 0)) +
  geom_segment(aes(x=0,xend=600,y=0,yend=0), color = 'azure3', linetype = 'dashed') +
  geom_segment(aes(x=600,xend=1000,y=15,yend=15), color = 'azure3') +
  geom_segment(aes(x=1072,xend=1172,y=15,yend=15), color = 'azure3') +
  geom_segment(aes(x=600,xend=600,y=0,yend=15), color = 'azure3') +
  geom_segment(aes(x=1000,xend=1000,y=0,yend=15), color = 'azure3') +
  geom_segment(aes(x=1072,xend=1072,y=0,yend=15), color = 'azure3') +
  geom_segment(aes(x=1172,xend=1172,y=0,yend=15), color = 'azure3') +
  geom_segment(aes(x=1000,xend=1072,y=0,yend=0), color = 'azure3', linetype = 'dashed') +
  geom_segment(aes(x=1172,xend=1272,y=0,yend=0), color = 'azure3', linetype = 'dashed') +
  theme_classic() +
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~rot_dir) +
  theme(axis.title.x = element_blank())

## NOTE: generlisation curves EE
dd <- d[phase == "Generalisation",
        .(mean(Endpoint_Error, na.rm = T), sd(Endpoint_Error, na.rm = T)/sqrt(.N)),
        .(cnd, target_deg, rot_dir)]

plt2 <- ggplot(dd, aes(x = target_deg, y = V1, colour = as.factor(cnd))) +
  geom_line(alpha = 0.6) +
  geom_errorbar(aes(ymin=V1-V2, ymax=V1+V2), width=10, alpha=0.6) +
  scale_y_continuous(limits = c(-10, 20)) +
  scale_x_continuous(breaks = c(-120, -60, 0, 60, 120), 
                     limits = c(-150, 150)) +
  facet_wrap(~rot_dir) +
  theme_classic() +
  labs(title = 'Generalisation across Targets',
       y = 'EE (in °)', 
       x = 'Target Position (in °)',
       colour = 'Variance \n Group') +
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())

## NOTE: overall performance all phases IMVE
dd <- d[, .(mean(imv_Error, na.rm=T), 
            sd(imv_Error, na.rm=T)/sqrt(16)), .(cnd, trial, rot_dir)]
plt3 <- ggplot(dd, aes(trial, V1, colour=as.factor(cnd))) +
  geom_line(alpha=0.4) +
  scale_y_continuous(limits = c(-15, 30)) +
  labs(x = 'Trial',
       y = 'IMVE (in °)',
       color = 'Variance \n Group') +
  scale_x_continuous(breaks=seq(0, 1200, 200), expand = c(0, 0)) +
  geom_segment(aes(x=0,xend=600,y=0,yend=0), color = 'azure3', linetype = 'dashed') +
  geom_segment(aes(x=600,xend=1000,y=15,yend=15), color = 'azure3') +
  geom_segment(aes(x=1072,xend=1172,y=15,yend=15), color = 'azure3') +
  geom_segment(aes(x=600,xend=600,y=0,yend=15), color = 'azure3') +
  geom_segment(aes(x=1000,xend=1000,y=0,yend=15), color = 'azure3') +
  geom_segment(aes(x=1072,xend=1072,y=0,yend=15), color = 'azure3') +
  geom_segment(aes(x=1172,xend=1172,y=0,yend=15), color = 'azure3') +
  geom_segment(aes(x=1000,xend=1072,y=0,yend=0), color = 'azure3', linetype = 'dashed') +
  geom_segment(aes(x=1172,xend=1272,y=0,yend=0), color = 'azure3', linetype = 'dashed') +
  facet_wrap(~rot_dir) +
  theme_classic() +
  theme(strip.background = element_rect(colour = "white", fill = "white"))

## NOTE: generlisation curves IMVE
dd <- d[phase == "Generalisation",
        .(mean(imv_Error, na.rm = T), sd(imv_Error, na.rm = T)/sqrt(.N)),
        .(cnd, target_deg, rot_dir)]

plt4 <- ggplot(dd, aes(x = target_deg, y = V1, colour = as.factor(cnd))) +
  geom_line(alpha = 0.6) +
  geom_errorbar(aes(ymin=V1-V2, ymax=V1+V2), width=10, alpha=0.6) +
  scale_y_continuous(limits = c(-10, 20)) +
  scale_x_continuous(breaks = c(-120, -60, 0, 60, 120), 
                     limits = c(-150, 150)) +
  facet_wrap(~rot_dir) +
  theme_classic() +
  labs(y = 'EE (in °)', 
       x = 'Target Position (in °)',
       colour = 'Variance \n Group') +
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
  theme(axis.title.y = element_blank())

ggarrange(plt1, plt2, plt3, plt4, nrow=2, ncol = 2, widths = c(2, 1), 
          legend = 'right', common.legend = T)
# ggsave('../figures/inspect_methodilogical_questions/comparison_imve_ee.png',
#        width=10, height=4)


## NOTE: Compute correlation between imve and ee
cor_total <- cor(d$Endpoint_Error, d$imv_Error, use = "complete.obs")
cor_subs <- d[, cor(Endpoint_Error, imv_Error, use = "complete.obs"), sub]
cor_subs

cor_subs[V1 %in% c(min(V1), max(V1))]

mean(d$Endpoint_Error, na.rm = T)
mean(d$imv_Error, na.rm = T)

cor_phase <- d[, .(mean_EE = mean(Endpoint_Error, na.rm = T), 
                   mean_IMVE = mean(imv_Error, na.rm = T),
                   sd_EE = sd(Endpoint_Error, na.rm = T), 
                   sd_IMVE = sd(imv_Error, na.rm = T)),
               .(phase)]
cor_phase


# Methodological Questions: Difference between bcee and ee
# --------------------------------------------------------
## NOTE: overall performance bcee
dd <- d[phase %in% c('Postbaseline', 'Training', 'Generalisation', 'Relearning', 'Washout')]
dd <- dd[, .(mean(bcee, na.rm=T), sd(bcee, na.rm=T)/sqrt(16)), .(cnd, trial, rot_dir)]
plt1 <- ggplot(dd, aes(trial, V1, color=as.factor(cnd))) +
  geom_line(alpha=0.4) +
  labs(title = 'Overall Performance across all Baseline Corrected Phases',
       x = 'Trial', 
       y = 'bcEE (in °)',
       color = 'Variance \n Group') +
  scale_y_continuous(limits = c(-10, 25)) +
  scale_x_continuous(breaks=seq(0, 1200, 200), expand = c(0, 0)) +
  geom_segment(aes(x=480,xend=600,y=0,yend=0), color = 'azure3', linetype = 'dashed') +
  geom_segment(aes(x=600,xend=1000,y=15,yend=15), color = 'azure3') +
  geom_segment(aes(x=1072,xend=1172,y=15,yend=15), color = 'azure3') +
  geom_segment(aes(x=600,xend=600,y=0,yend=15), color = 'azure3') +
  geom_segment(aes(x=1000,xend=1000,y=0,yend=15), color = 'azure3') +
  geom_segment(aes(x=1072,xend=1072,y=0,yend=15), color = 'azure3') +
  geom_segment(aes(x=1172,xend=1172,y=0,yend=15), color = 'azure3') +
  geom_segment(aes(x=1000,xend=1072,y=0,yend=0), color = 'azure3', linetype = 'dashed') +
  geom_segment(aes(x=1172,xend=1272,y=0,yend=0), color = 'azure3', linetype = 'dashed') +
  theme_classic() +
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~rot_dir) +
  theme(axis.title.x = element_blank())

## NOTE: generalisation curves bcee
dd <- d[phase == "Generalisation",
        .(mean(bcee, na.rm = T), sd(bcee, na.rm = T)/sqrt(16)),
        .(cnd, target_deg, rot_dir)]

plt2 <- ggplot(dd, aes(x = target_deg, y = V1, colour = as.factor(cnd))) +
  geom_line(alpha=0.6) +
  geom_errorbar(aes(ymin=V1-V2, ymax=V1+V2), width=10, alpha=0.6) +
  scale_y_continuous(limits = c(-5, 15)) +
  scale_x_continuous(breaks = c(-120, -60, 0, 60, 120), 
                     limits = c(-150, 150)) +
  facet_wrap(~rot_dir) +
  theme_classic() +
  labs(title = 'Generalisation across Targets',
       y = 'bcEE (in °)', 
       x = 'Target Position (in °)',
       colour = 'Variance \n Group') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

## NOTE: overall performance ee
dd <- d[phase %in% c('Postbaseline', 'Training', 'Generalisation', 'Relearning', 'Washout')]
dd <- dd[, .(mean(Endpoint_Error, na.rm=T), sd(Endpoint_Error, na.rm=T)/sqrt(16)), .(cnd, trial, rot_dir)]
plt3 <- ggplot(dd, aes(trial, V1, color=as.factor(cnd))) +
  geom_line(alpha=0.4) +
  labs(x = 'Trial', 
       y = 'EE (in °)',
       color = 'Variance \n Group') +
  scale_y_continuous(limits = c(-10, 25)) +
  scale_x_continuous(breaks=seq(0, 1200, 200), expand = c(0, 0)) +
  geom_segment(aes(x=480,xend=600,y=0,yend=0), color = 'azure3', linetype = 'dashed') +
  geom_segment(aes(x=600,xend=1000,y=15,yend=15), color = 'azure3') +
  geom_segment(aes(x=1072,xend=1172,y=15,yend=15), color = 'azure3') +
  geom_segment(aes(x=600,xend=600,y=0,yend=15), color = 'azure3') +
  geom_segment(aes(x=1000,xend=1000,y=0,yend=15), color = 'azure3') +
  geom_segment(aes(x=1072,xend=1072,y=0,yend=15), color = 'azure3') +
  geom_segment(aes(x=1172,xend=1172,y=0,yend=15), color = 'azure3') +
  geom_segment(aes(x=1000,xend=1072,y=0,yend=0), color = 'azure3', linetype = 'dashed') +
  geom_segment(aes(x=1172,xend=1272,y=0,yend=0), color = 'azure3', linetype = 'dashed') +
  theme_classic() +
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~rot_dir)

## NOTE: generlisation curves ee
dd <- d[phase == "Generalisation",
        .(mean(Endpoint_Error, na.rm = T), sd(Endpoint_Error, na.rm = T)/sqrt(.N)),
        .(cnd, target_deg, rot_dir)]

plt4 <- ggplot(dd, aes(x = target_deg, y = V1, colour = as.factor(cnd))) +
  geom_line(alpha = 0.6) +
  geom_errorbar(aes(ymin=V1-V2, ymax=V1+V2), width=10, alpha=0.6) +
  scale_y_continuous(limits = c(-5, 15)) +
  scale_x_continuous(breaks = c(-120, -60, 0, 60, 120), 
                     limits = c(-150, 150)) +
  facet_wrap(~rot_dir) +
  theme_classic() +
  labs(y = 'EE (in °)', 
       x = 'Target Position (in °)',
       colour = 'Variance \n Group') +
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
  theme(axis.title.y = element_blank())

ggarrange(plt1, plt2, plt3, plt4, nrow=2, ncol = 2, widths = c(2, 1), 
          legend = 'right', common.legend = T)
# ggsave('../figures/inspect_methodilogical_questions/comparison_bcee_ee.png',
#        width=10, height=4)


# Inspect: Reaction time as explanation for the difference between imv vs. EE?
# ----------------------------------------------------------------------------
## NOTE: Performance across all trials
dd <- d[RT > 0.6, RT := NaN]
dd <- dd[phase %in% c("Training", "Generalisation", "Relearning", "Washout"), 
         .(mean(RT, na.rm=T), sd(RT, na.rm=T)/sqrt(16)), .(cnd, trial, rot_dir)]
plt1 <- ggplot(dd, aes(trial, V1, colour=as.factor(cnd))) +
  geom_line(alpha=0.4) +
  labs(title = 'RTs across all Baseline Corrected Phases',
       x = 'Trial', 
       y = 'RT (in s)',
       color = 'Variance \n Group') +
  facet_wrap(~rot_dir, nrow=2, scales = 'free') +
  theme_classic() +
  scale_x_continuous(breaks=seq(0, 1200, 200)) +
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(colour = 'azure2'))

## NOTE: RT in generalisation
dd <- d[RT > 0.6, RT := NaN]
dd <- dd[phase == "Generalisation",
        .(mean(RT, na.rm = T), sd(RT, na.rm = T)/sqrt(16)),
        .(cnd, target_deg, rot_dir)]

plt2 <- ggplot(dd, aes(x = target_deg, y = V1, colour = as.factor(cnd))) +
  geom_line(alpha = 0.6) +
  geom_errorbar(aes(ymin=V1-V2, ymax=V1+V2), width=10, alpha = 0.4) +
  scale_y_continuous(limits = c(0.2, 0.4)) +
  scale_x_continuous(breaks = c(-120, -60, 0, 60, 120), 
                     limits = c(-150, 150)) +
  facet_wrap(~rot_dir, nrow=2, scales = 'free') +
  theme_classic() +
  labs(title = 'RTs in Generalisation',
       y = 'RT (in s)', 
       x = 'Target Position (in °)',
       colour = 'Variance \n Group') +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        panel.grid.major = element_line(colour = 'azure2')) +
  theme(strip.background = element_rect(colour = "white", fill = 'white'))

ggarrange(plt1, plt2, ncol = 2, 
          legend = 'right', common.legend = T)
# ggsave('../figures/inspect_difference_rotation/comparison_RTs_cw_ccw.png',
#        width=10, height=4)


## NOTE: Testing RT difference between cw and ccw 
## NOTE: for training
dd <- d[phase == 'Training']
m.pre <- lm(RT ~ as.factor(cnd), data = dd)
# and store the residual HE values. 
dd$RT.r <- residuals(m.pre)

m.RT <- lm(RT.r ~ as.factor(rot_dir), data = dd) 
summary(m.RT)

## NOTE: for generalisation
dd <- d[phase == 'Generalisation']
m.pre <- lm(RT ~ as.factor(cnd) + as.factor(target_deg), data = dd)
# and store the residual HE values. 
dd$RT.r <- residuals(m.pre)

m.RT <- lm(RT.r ~ as.factor(rot_dir), data = dd) 
summary(m.RT)

## NOTE: for washout
dd <- d[phase == 'Washout']
m.pre <- lm(RT ~ as.factor(cnd), data = dd)
# and store the residual HE values. 
dd$RT.r <- residuals(m.pre)

m.RT <- lm(RT.r ~ as.factor(rot_dir), data = dd) 
summary(m.RT)


# Inspect: Difference between implicit and explicit
# -------------------------------------------------
## NOTE: overall performance bcee
dd <- d[phase %in% c('Training', 'Generalisation', 'Relearning', 'Washout')]
dd <- dd[, .(mean(bcee, na.rm=T), sd(bcee, na.rm=T)/sqrt(16)), .(cnd, trial, rot_dir, aware)]

plt1 <- ggplot(dd, aes(trial, V1, color=as.factor(cnd))) +
  geom_line(alpha=0.4) +
  labs(title = 'Overall Performance across all Baseline Corrected Phases',
       x = 'Trial', 
       y = 'Baseline corrected EE (in °)',
       color = 'Variance \n Group') +
  scale_y_continuous(limits = c(-15, 40)) +
  geom_segment(aes(x=600,xend=1000,y=15,yend=15), color = 'azure3') +
  geom_segment(aes(x=1072,xend=1172,y=15,yend=15), color = 'azure3') +
  geom_segment(aes(x=600,xend=600,y=-10,yend=15), color = 'azure3') +
  geom_segment(aes(x=1000,xend=1000,y=-10,yend=15), color = 'azure3') +
  geom_segment(aes(x=1072,xend=1072,y=-10,yend=15), color = 'azure3') +
  geom_segment(aes(x=1172,xend=1172,y=-10,yend=15), color = 'azure3') +
  theme_classic() +
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~aware*rot_dir, scales = 'free',
             labeller = label_wrap_gen(multi_line=FALSE))

## NOTE: generalisation curve bcee
dd <- d[phase == "Generalisation",
        .(mean(bcee, na.rm = T), sd(bcee, na.rm = T)/sqrt(.N)),
        .(cnd, target_deg, rot_dir, aware)]

plt2 <- ggplot(dd, aes(x = target_deg, y = V1, colour = as.factor(cnd))) +
  geom_line(alpha = 0.6) +
  geom_errorbar(aes(ymin=V1-V2, ymax=V1+V2), width=10, alpha = 0.6) +
  labs(title = 'Generalisation across Targets',
       x = 'Target Position (in °)', 
       y = 'Baseline corrected EE (in °)',
       color = 'Variance \n Group') +
  scale_x_continuous(breaks = c(-120, -60, 0, 60, 120), limits = c(-150, 150)) +
  facet_wrap(~aware*rot_dir, scales = 'free',
             labeller = label_wrap_gen(multi_line=FALSE)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
  theme(axis.title.y = element_blank())

ggarrange(plt1, plt2, ncol = 2, widths = c(2, 1), 
          legend = 'right', common.legend = T)
# ggsave('../figures/inspect_difference_rotation/comparison_bcee_aware.png',
#        width=10, height=4)


# ## NOTE: Test for significance of additional awereness factor in training 
# dd <- d[phase %in% c('Training')]
# ml1 <- lm(bcee ~ as.factor(rot_dir) + as.factor(cnd),
#           data = dd)
# ml2 <- lm(bcee ~ as.factor(rot_dir) + as.factor(cnd) + as.factor(aware),
#           data = dd)
# 
# summary(ml1)
# summary(ml2)
# anova(ml1, ml2) 
# 
# ## NOTE: Test for significance of additional awereness factor in washout
# dd <- d[phase %in% c('Generalisation')]
# ml1 <- lm(bcee ~ as.factor(rot_dir) + as.factor(cnd),
#           data = dd)
# ml2 <- lm(bcee ~ as.factor(rot_dir) + as.factor(cnd) + as.factor(aware),
#           data = dd)
# 
# summary(ml1)
# summary(ml2)
# anova(ml1, ml2)
# 
# ## NOTE: Test for significance of additional awereness factor in washout
# dd <- d[phase %in% c('Washout')]
# ml1 <- lm(bcee ~ as.factor(rot_dir) + as.factor(cnd),
#           data = dd)
# ml2 <- lm(bcee ~ as.factor(rot_dir) + as.factor(cnd) + as.factor(aware),
#           data = dd)
# 
# summary(ml1)
# summary(ml2)
# anova(ml1, ml2)


# Inspect: Analysing Prebaseline
# ------------------------------
## NOTE: take a look at prebaseline histograms
dd <- d[phase == 'Prebaseline', mean(Endpoint_Error, na.rm=T), .(cnd, trial, rot_dir)]
ggplot(dd, aes(V1, fill=as.factor(cnd))) +
  geom_histogram(alpha=0.75) +
  theme_classic() +
  labs(y = 'Number of reaches', 
       x = 'EE (in °)',
       fill = 'Variance \n Group') +
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
  scale_x_continuous(limits = c(-12, 4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~rot_dir*cnd, scales = 'free',
             labeller = label_wrap_gen(multi_line = FALSE))
# ggsave('../figures/inspect_prebaseline/prebaseline_hist.png',
#        width=10, height=4)


## NOTE: look at prebaseline bias per target
dd <- d[phase == 'Prebaseline',
        .(mean(bcee, na.rm=T), sd(Endpoint_Error, na.rm=T)/sqrt(16)),
        .(cnd, rot_dir, target_deg)]
ggplot(dd, aes(target_deg, V1, colour=as.factor(cnd))) +
  geom_point(alpha=0.4) +
  geom_line(alpha=0.4) +
  geom_errorbar(aes(ymin=V1-V2, ymax=V1+V2), alpha=0.4) +
  theme_classic() +
  labs(y = 'EE (in °)', 
       x = 'Target Position (in °)',
       colour = 'Variance \n Group') +
  scale_x_continuous(breaks = c(-150, -120, -90, -60, -30, 0, 
                                30, 60, 90, 120, 150), 
                     limits = c(-150, 150)) +
  scale_y_continuous(limits = c(-10, 7)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
  facet_wrap(~rot_dir*cnd, scales = 'free',
             labeller = label_wrap_gen(multi_line = FALSE))
# ggsave('../figures/inspect_prebaseline/prebaseline_bias.png',
#        width=10, height=4)
