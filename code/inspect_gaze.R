library(data.table)
library(ggplot2)
library(tidyverse)
library(ggpubr)

rm( list = ls() )
setwd('D:/Studium/Auslandsstudium/TuitionWaver_Master/Masterthesis/Analysis/Exp_Variance/Modelling/sorted/code')
# --------------------------------------------------------------------------------------------------------------

# # Extracting csv files for GazeData
# # ---------------------------------
# setwd("../../../MatlabAnalysis/AnalysisData/Exp_Variance_MissingsReplaced_CorrectRTs")
# dir_name <- "../../../Modelling/sorted/data"
# 
# myfolders <- list.dirs(recursive = FALSE)
# myfolders <- myfolders[1:8]
# 
# for (h in 1:length(myfolders)){
#     subjdatafiles <- list.files(myfolders[h], pattern = "*GazeData.csv")
#     for (g in 1:length(subjdatafiles)){
#         file.copy(paste0(myfolders[h], "/", subjdatafiles[g]), dir_name)
#     }
# }
# 
# setwd("../../../Modelling/sorted/code")
# # ---------------------------------------------------------------------------------------------------------------

# Transform Dataframe Gaze Data
# -----------------------------
phase_names_f <- c('Prebaseline_GazeData', 'Familiarisation_GazeData',
                   'Baseline_NFB_GazeData', 'Baseline_GazeData',
                   'Postbaseline_GazeData', 'Training_GazeData',
                   'Generalisation_GazeData', 'Relearning_GazeData',
                   'Washout_GazeData')

phase_names <- c('Prebaseline', 'Familiarisation', 'Baseline_NFB', 'Baseline_FB',
                 'Postbaseline', 'Training', 'Generalisation', 'Relearning',
                 'Washout')

phase_order <- 1:9

ldf_rec <- list()
for (i in 1:length(phase_names_f)) {
    f_names <- list.files(paste('../data', sep=''),
                          pattern=paste(phase_names_f[i]),
                          full.names=TRUE)
    ldf <- lapply(seq_along(f_names),
                  function(j) {
                      z<-fread(f_names[j]);
                      setnames(z, c('trial', 'target', 'x', 'y', 't'));
                      })
    ldf <- lapply(seq_along(ldf), function(z) ldf[[z]][, sub := rep(z, .N)])
    ldf <- lapply(seq_along(ldf), function(z) ldf[[z]][, phase := rep(phase_names[i], .N)])
    ldf <- lapply(seq_along(ldf), function(z) ldf[[z]][, phase_order := rep(phase_order[i], .N)])
    ldf_rec <- c(ldf_rec, ldf)
}

d.Ga <- rbindlist(ldf_rec)
# -------------------------------------------------------------------------------------------------------------------

# # Extracting csv files for TrajData
# # ---------------------------------
# 
# setwd("../../../MatlabAnalysis/AnalysisData/Exp_Variance_MissingsReplaced_CorrectRTs")
# dir_name <- "../../../Modelling/sorted/data"
# 
# myfolders <- list.dirs(recursive = FALSE)
# myfolders <- myfolders[1:8]
# 
# for (h in 1:length(myfolders)){
#     subjdatafiles <- list.files(myfolders[h], pattern = "*TrajData.csv")
#     for (g in 1:length(subjdatafiles)){
#         file.copy(paste0(myfolders[h], "/", subjdatafiles[g]), dir_name)
#     }
# }
# 
# setwd("../../../Modelling/sorted/code")
# # --------------------------------------------------------------------------------------------------------------------

# Transform Dataframe Trajectory Data
# -----------------------------------
phase_names_f <- c('Prebaseline_TrajData', 'Familiarisation_TrajData',
                   'Baseline_NFB_TrajData', 'Baseline_TrajData',
                   'Postbaseline_TrajData', 'Training_TrajData',
                   'Generalisation_TrajData', 'Relearning_TrajData',
                   'Washout_TrajData')

phase_names <- c('Prebaseline', 'Familiarisation', 'Baseline_NFB', 'Baseline_FB',
                 'Postbaseline', 'Training', 'Generalisation', 'Relearning',
                 'Washout')

phase_order <- 1:9

ldf_rec <- list()
for (i in 1:length(phase_names_f)) {
    f_names <- list.files(paste('../data', sep=''),
                          pattern=paste(phase_names_f[i]),
                          full.names=TRUE)
    ldf <- lapply(seq_along(f_names),
                  function(j) {
                      z<-fread(f_names[j]);
                      setnames(z, c('trial', 'target', 'x', 'y'));
                  })
    ldf <- lapply(seq_along(ldf), function(z) ldf[[z]][, sub := rep(z, .N)])
    ldf <- lapply(seq_along(ldf), function(z) ldf[[z]][, phase := rep(phase_names[i], .N)])
    ldf <- lapply(seq_along(ldf), function(z) ldf[[z]][, phase_order := rep(phase_order[i], .N)])
    ldf_rec <- c(ldf_rec, ldf)
}

d.Tr <- rbindlist(ldf_rec)

d <- d.Ga[, `:=`(x.ga = x, y.ga = y, x.tr = d.Tr$x, y.tr = d.Tr$y)]
d <- d[, `:=`(x = NULL, y = NULL)]

sub_ccw <- c(1, 3, 4, 6, 7, 10, 12, 14, 15, 16, 17, 19, 21, 24, 29, 31, 33, 34,
             35, 36, 37, 41, 43, 47)
sub_expl <- c(1, 2, 3, 5, 9, 10, 16, 18, 19, 20, 24, 28, 29, 30, 31, 33, 34,
              40, 41, 44, 45, 46, 47, 48)

d <- d[,phase := factor(phase, levels = c('Prebaseline', 'Familiarisation', 'Baseline_FB',
                                          'Baseline_NFB', 'Postbaseline', 'Training',
                                          'Generalisation', 'Relearning', 'Washout'))]
# ========================================================================================================================

# Plots: Trajectories only
# ------------------------
# NOTE: Compte trial time
d[, t_trial := t - min(t), .(sub, trial, target, phase)]

# NOTE: Define Target Locations
targets = data.table(
    X = c(0, 0, 0.0500000000000000, 0.0866025403784000, 0.100000000000000, 0.0866025403784000, 
          0.0500000000000000, 1.22464679915000e-17, -0.0500000000000000, -0.0866025403784000, 
          -0.100000000000000, -0.0866025403784000, -0.0500000000000000),
    Y = c(0.2, 0.300000000000000, 0.286602540378000, 0.250000000000000, 0.200000000000000,
          0.150000000000000, 0.113397459622000, 0.100000000000000, 0.113397459622000,
          0.150000000000000, 0.200000000000000, 0.250000000000000, 0.286602540378000),
    target_num = 0:12
    )

targets_train = data.table(
    X = c(0, 0),
    Y = c(0.2, 0.300000000000000),
    target_num = 0:1
)


# Comparison IMV error vs. EE
# ---------------------------
# NOTE: Plot Traj and Gaze for sub with highest corr between imv and ee (46) and lowest (6)
ddd <- d[sub %in% c(6, 46) & phase == 'Baseline_NFB', 
         .(mean(x.tr, na.rm = T), mean(y.tr, na.rm = T)), 
         .(sub, target, t_trial)]
ggplot(ddd[t_trial < 0.6], aes(V1, V2)) + 
    geom_point(colour=alpha('blue', 0.15)) +
    geom_point(data = targets, mapping = aes(x = X, y = Y), colour = 'red') +
    xlim(-0.1, 0.1) +
    ylim(0.1, 0.3) +
    facet_wrap(~sub) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5)) +
    labs(x = 'X coordinates (in m)',
         y = 'Y coordinates (in m)')
# ggsave('../figures/inspect_methodilogical_questions/trajdata_comp_imve_ee.png', 
#        width=10, height=5.5)


# NOTE: Plot Traj and Gaze for 4 subs with highest corr between imv and ee and lowest
dd <- d[sub %in% c(6, 17, 21, 26, 30, 5, 45, 46) & phase == 'Baseline_NFB', 
         .(mean(x.tr, na.rm = T), mean(y.tr, na.rm = T)), 
         .(factor(sub, levels = c(6, 17, 21, 26, 30, 5, 45, 46)), target, t_trial)]
# dd <- dd[,sub := factor(sub, levels = c(6, 17, 21, 26, 30, 5, 45, 46))]
ggplot(dd[t_trial < 0.6], aes(V1, V2)) + 
    geom_point(colour=alpha('blue', 0.15)) +
    geom_point(data = targets, mapping = aes(x = X, y = Y), colour = 'red') +
    xlim(-0.1, 0.1) +
    ylim(0.1, 0.3) +
    facet_wrap(~factor, nrow = 2) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5)) +
    labs(x = 'X coordinates (in m)',
         y = 'Y coordinates (in m)')
# ggsave('../figures/inspect_methodilogical_questions/trajdata_comp_cor-imve-ee_subs.png', 
#        width=10, height=5.5)


# Prebaseline Trajectories
# ------------------------
# NOTE: Example trajectories for first 4 subs
ddd <- d[sub %in% 1:4 & phase == 'Prebaseline', 
         .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), .(sub, target, t_trial)]
ggplot(aes(V3, V4), data = ddd) +
    geom_point(colour=alpha('blue', 0.05)) +
    geom_point(data = targets, mapping = aes(x = X, y = Y), colour = 'red') +
    xlim(-0.1, 0.1) +
    ylim(0.1, 0.3) +
    facet_wrap(~sub) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5)) +
    labs(x = 'X coordinates (in m)',
         y = 'Y coordinates (in m)')
# ggsave('../figures/inspect_prebaseline/trajdata_prebase_examplsubs.png', 
#        width=5, height=5)

# NOTE: Mean trajcetoreis across all subs
ddd <- d[phase == 'Prebaseline', .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), .(target, t_trial)]
ggplot(aes(V3, V4), data = ddd) +
    geom_point(colour=alpha('blue', 0.05)) +
    geom_point(data = targets, mapping = aes(x = X, y = Y), colour = 'red') +
    xlim(-0.1, 0.1) +
    ylim(0.1, 0.3) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5)) +
    labs(x = 'X coordinates (in m)',
         y = 'Y coordinates (in m)')
# ggsave('../figures/inspect_prebaseline/trajdata_prebase_mean_subs.png', 
#        width=5, height=5)
# ------------------------------------------------------------------------------------------------------

# Plots: Comparison Gaze and Trajectory Data
# ------------------------------------------
# NOTE: Define Colours for the legend
colours <- c('Trajectories' = 'blue', 'Gaze' = 'orange', 'Targets' = 'red')


# Traj and Gaze Development over time
# -----------------------------------
## NOTE: ccw training 50ms exlicit
dd <- d[sub %in% sub_ccw & phase == 'Training' & sub %in% sub_expl, 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), 
        .(sub, target, t_trial)]

plt1 <- ggplot(dd[t_trial < 0.05], aes(V1, V2)) + 
    geom_point(aes(x = V1, y = V2, colour='Gaze'), alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4, colour='Trajectories'), alpha = 0.05) +
    geom_point(data = targets_train, mapping = aes(x = X, y = Y, colour = 'Targets')) +
    xlim(-0.05, 0.05) +
    ylim(0.18, 0.3) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(x = '100 ms, X coordinates (in m)',
         y = 'Y coordinates (in m)') +
    scale_color_manual(name = '', values = colours)

## NOTE: ccw training 200ms exlicit
dd <- d[sub %in% sub_ccw & phase == 'Training' & sub %in% sub_expl, 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), 
        .(sub, target, t_trial)]

plt2 <- ggplot(dd[t_trial < 0.2], aes(V1, V2)) + 
    geom_point(aes(x = V1, y = V2, colour='Gaze'), alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4, colour='Trajectories'), alpha = 0.05) +
    geom_point(data = targets_train, mapping = aes(x = X, y = Y, colour = 'Targets')) +
    xlim(-0.05, 0.05) +
    ylim(0.18, 0.3) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(x = '200 ms, X coordinates (in m)',
         y = 'Y coordinates (in m)') +
    scale_color_manual(name = '', values = colours)

## NOTE: ccw training 400ms exlicit
dd <- d[sub %in% sub_ccw & phase == 'Training' & sub %in% sub_expl, 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), 
        .(sub, target, t_trial)]

plt3 <- ggplot(dd[t_trial < 0.4], aes(V1, V2)) + 
    geom_point(aes(x = V1, y = V2, colour='Gaze'), alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4, colour='Trajectories'), alpha = 0.05) +
    geom_point(data = targets_train, mapping = aes(x = X, y = Y, colour = 'Targets')) +
    xlim(-0.05, 0.05) +
    ylim(0.18, 0.3) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(x = '400 ms, X coordinates (in m)',
         y = 'Y coordinates (in m)') +
    scale_color_manual(name = '', values = colours)

## NOTE: ccw training 600ms exlicit
dd <- d[sub %in% sub_ccw & phase == 'Training' & sub %in% sub_expl, 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), 
        .(sub, target, t_trial)]

plt4 <- ggplot(dd[t_trial < 0.6], aes(V1, V2)) + 
    geom_point(aes(x = V1, y = V2, colour='Gaze'), alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4, colour='Trajectories'), alpha = 0.05) +
    geom_point(data = targets_train, mapping = aes(x = X, y = Y, colour = 'Targets')) +
    xlim(-0.05, 0.05) +
    ylim(0.18, 0.3) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(x = '600ms, X coordinates (in m)',
         y = 'Y coordinates (in m)') +
    scale_color_manual(name = '', values = colours)

row1 <- ggarrange(plt1, plt2, plt3, plt4, nrow = 1, ncol = 4, 
          legend = 'top', common.legend = T)
row1 <- annotate_figure(row1, top = 'Training, Explicit, CCW')


## NOTE: cw training 50ms exlicit
dd <- d[!(sub %in% sub_ccw) & phase == 'Training' & sub %in% sub_expl, 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), 
        .(sub, target, t_trial)]

plt5 <- ggplot(dd[t_trial < 0.05], aes(V1, V2)) + 
    geom_point(aes(x = V1, y = V2), colour='orange', alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4), colour='blue', alpha = 0.05) +
    geom_point(data = targets_train, mapping = aes(x = X, y = Y), colour = 'red') +
    xlim(-0.05, 0.05) +
    ylim(0.18, 0.3) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(x = '100ms, X coordinates (in m)',
         y = 'Y coordinates (in m)')

## NOTE: cw training 200ms exlicit
dd <- d[!(sub %in% sub_ccw) & phase == 'Training' & sub %in% sub_expl, 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), 
        .(sub, target, t_trial)]

plt6 <- ggplot(dd[t_trial < 0.2], aes(V1, V2)) + 
    geom_point(aes(x = V1, y = V2), colour='orange', alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4), colour='blue', alpha = 0.05) +
    geom_point(data = targets_train, mapping = aes(x = X, y = Y), colour = 'red') +
    xlim(-0.05, 0.05) +
    ylim(0.18, 0.3) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(x = '200 ms, X coordinates (in m)',
         y = 'Y coordinates (in m)')

## NOTE: cw training 400ms exlicit
dd <- d[!(sub %in% sub_ccw) & phase == 'Training' & sub %in% sub_expl, 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), 
        .(sub, target, t_trial)]

plt7 <- ggplot(dd[t_trial < 0.4], aes(V1, V2)) + 
    geom_point(aes(x = V1, y = V2), colour='orange', alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4), colour='blue', alpha = 0.05) +
    geom_point(data = targets_train, mapping = aes(x = X, y = Y), colour = 'red') +
    xlim(-0.05, 0.05) +
    ylim(0.18, 0.3) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(x = '400 ms, X coordinates (in m)',
         y = 'Y coordinates (in m)')

## NOTE: cw training 600ms exlicit
dd <- d[!(sub %in% sub_ccw) & phase == 'Training' & sub %in% sub_expl, 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), 
        .(sub, target, t_trial)]

plt8 <- ggplot(dd[t_trial < 0.6], aes(V1, V2)) + 
    geom_point(aes(x = V1, y = V2), colour='orange', alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4), colour='blue', alpha = 0.05) +
    geom_point(data = targets_train, mapping = aes(x = X, y = Y), colour = 'red') +
    xlim(-0.05, 0.05) +
    ylim(0.18, 0.3) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(x = '600 ms, X coordinates (in m)',
         y = 'Y coordinates (in m)')

row2 <- ggarrange(plt5, plt6, plt7, plt8, nrow = 1, ncol = 4, 
          legend = 'right', common.legend = T)
row2 <- annotate_figure(row2, top = 'Explicit, CW')


## NOTE: ccw training 50ms implicit
dd <- d[sub %in% sub_ccw & phase == 'Training' & !(sub %in% sub_expl), 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), 
        .(sub, target, t_trial)]

plt9 <- ggplot(dd[t_trial < 0.05], aes(V1, V2)) + 
    geom_point(aes(x = V1, y = V2), colour='orange', alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4), colour='blue', alpha = 0.05) +
    geom_point(data = targets_train, mapping = aes(x = X, y = Y), colour = 'red') +
    xlim(-0.05, 0.05) +
    ylim(0.18, 0.3) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(x = '100 ms, X coordinates (in m)',
         y = 'Y coordinates (in m)')

## NOTE: ccw training 200ms implicit
dd <- d[sub %in% sub_ccw & phase == 'Training' & !(sub %in% sub_expl), 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), 
        .(sub, target, t_trial)]

plt10 <- ggplot(dd[t_trial < 0.2], aes(V1, V2)) + 
    geom_point(aes(x = V1, y = V2), colour='orange', alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4), colour='blue', alpha = 0.05) +
    geom_point(data = targets_train, mapping = aes(x = X, y = Y), colour = 'red') +
    xlim(-0.05, 0.05) +
    ylim(0.18, 0.3) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(x = '200 ms, X coordinates (in m)',
         y = 'Y coordinates (in m)')

## NOTE: ccw training 400ms implicit
dd <- d[sub %in% sub_ccw & phase == 'Training' & !(sub %in% sub_expl), 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), 
        .(sub, target, t_trial)]

plt11 <- ggplot(dd[t_trial < 0.4], aes(V1, V2)) + 
    geom_point(aes(x = V1, y = V2), colour='orange', alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4), colour='blue', alpha = 0.05) +
    geom_point(data = targets_train, mapping = aes(x = X, y = Y), colour = 'red') +
    xlim(-0.05, 0.05) +
    ylim(0.18, 0.3) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(x = '400 ms, X coordinates (in m)',
         y = 'Y coordinates (in m)')

## NOTE: ccw training 600ms implicit
dd <- d[sub %in% sub_ccw & phase == 'Training' & !(sub %in% sub_expl), 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), 
        .(sub, target, t_trial)]

plt12 <- ggplot(dd[t_trial < 0.6], aes(V1, V2)) + 
    geom_point(aes(x = V1, y = V2), colour='orange', alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4), colour='blue', alpha = 0.05) +
    geom_point(data = targets_train, mapping = aes(x = X, y = Y), colour = 'red') +
    xlim(-0.05, 0.05) +
    ylim(0.18, 0.3) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(x = '600 ms, X coordinates (in m)',
         y = 'Y coordinates (in m)')

row3 <- ggarrange(plt9, plt10, plt11, plt12, nrow = 1, ncol = 4, 
          legend = 'right', common.legend = T)
row3 <- annotate_figure(row3, top = 'Implicit, CCW')


## NOTE: cw training 50ms implicit
dd <- d[!(sub %in% sub_ccw) & phase == 'Training' & !(sub %in% sub_expl), 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), 
        .(sub, target, t_trial)]

plt13 <- ggplot(dd[t_trial < 0.05], aes(V1, V2)) + 
    geom_point(aes(x = V1, y = V2), colour='orange', alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4), colour='blue', alpha = 0.05) +
    geom_point(data = targets_train, mapping = aes(x = X, y = Y), colour = 'red') +
    xlim(-0.05, 0.05) +
    ylim(0.18, 0.3) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    labs(x = '50 ms, X coordinates (in m)',
         y = 'Y coordinates (in m)')

## NOTE: cw training 200ms implicit
dd <- d[!(sub %in% sub_ccw) & phase == 'Training' & !(sub %in% sub_expl), 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), 
        .(sub, target, t_trial)]

plt14 <- ggplot(dd[t_trial < 0.2], aes(V1, V2)) + 
    geom_point(aes(x = V1, y = V2), colour='orange', alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4), colour='blue', alpha = 0.05) +
    geom_point(data = targets_train, mapping = aes(x = X, y = Y), colour = 'red') +
    xlim(-0.05, 0.05) +
    ylim(0.18, 0.3) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank()) +
    labs(x = '200 ms, X coordinates (in m)',
         y = 'Y coordinates (in m)')

## NOTE: cw training 400ms implicit
dd <- d[!(sub %in% sub_ccw) & phase == 'Training' & !(sub %in% sub_expl), 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), 
        .(sub, target, t_trial)]

plt15 <- ggplot(dd[t_trial < 0.4], aes(V1, V2)) + 
    geom_point(aes(x = V1, y = V2), colour='orange', alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4), colour='blue', alpha = 0.05) +
    geom_point(data = targets_train, mapping = aes(x = X, y = Y), colour = 'red') +
    xlim(-0.05, 0.05) +
    ylim(0.18, 0.3) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank()) +
    labs(x = '400 ms, X coordinates (in m)',
         y = 'Y coordinates (in m)')

## NOTE: cw training 600ms implicit
dd <- d[!(sub %in% sub_ccw) & phase == 'Training' & !(sub %in% sub_expl), 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), 
        .(sub, target, t_trial)]

plt16 <- ggplot(dd[t_trial < 0.6], aes(V1, V2)) + 
    geom_point(aes(x = V1, y = V2), colour='orange', alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4), colour='blue', alpha = 0.05) +
    geom_point(data = targets_train, mapping = aes(x = X, y = Y), colour = 'red') +
    xlim(-0.05, 0.05) +
    ylim(0.18, 0.3) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank()) +
    labs(x = '600 ms, X coordinates (in m)',
         y = 'Y coordinates (in m)')

row4 <- ggarrange(plt13, plt14, plt15, plt16, nrow = 1, ncol = 4, 
                  legend = 'right', common.legend = T)
row4 <- annotate_figure(row4, top = 'Implicit, CW')

final <- ggarrange(row1, row2, row3, row4, nrow = 4, ncol = 1, 
          legend = 'right', common.legend = T)
annotate_figure(final, left = 'Y coordinates (in m)')
# ggsave('../figures/inspect_difference_rotation/trajgaze_data_over_time.png', 
#        width=10, height=10)


# Inspect Familiarisation
# -----------------------
# NOTE: Select only subjects with good eye tracking quality
sub_bad <- c(1, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16, 17, 18, 20, 21, 22, 23, 25, 
             26, 27, 28, 31, 32, 33, 34, 36, 38, 40, 42, 43, 44, 45, 48)
d[sub %in% sub_bad, `:=`(x.ga = NA, y.ga = NA, x.tr = NA, y.tr = NA)]


# NOTE: Plot Familiarisation
dd <- d[!(sub %in% sub_bad) & phase == 'Familiarisation', 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), 
        .(target, t_trial)]
ggplot(dd[t_trial < 0.6], aes(V1, V2)) +
    geom_point(aes(x = V1, y = V2), colour=alpha('orange', 0.1)) + 
    geom_point(aes(x = V3, y = V4), colour=alpha('blue', 0.1) ) +
    geom_point(data = targets, mapping = aes(x = X, y = Y), colour = 'red') +
    xlim(-0.1, 0.1) +
    ylim(0.1, 0.3) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5)) +
    labs(x = 'X coordinates (in m)',
         y = 'Y coordinates (in m)')
# ggsave('../figures/inspect_difference_rotation/gazetrajdata_all_subs_famil.png', 
#        width=5, height=5)


# Comparison CW vs. CCW
# ---------------------
# NOTE: In rotation relevant phases only
colours <- c('Trajectories' = 'blue', 'Gaze' = 'orange', 'Targets' = 'red')

# NOTE: CCW
dd <- d[sub %in% sub_ccw & !(sub %in% sub_bad) & 
            phase %in% c('Training', 'Generalisation', 'Relearning', 'Washout'), 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), .(phase, target, t_trial)]

plt1 <- ggplot(dd[t_trial < 0.6], aes(V1, V2)) +
    geom_point(aes(x = V1, y = V2, colour='Gaze'), alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4, colour='Trajectories'), alpha = 0.05) +
    geom_point(data = targets, mapping = aes(x = X, y = Y, colour = 'Targets')) +
    xlim(-0.1, 0.1) +
    ylim(0.1, 0.3) +
    facet_wrap(~phase) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = 'Subjects with CCW Rotation',
         x = 'X coordinates (in m)',
         y = 'Y coordinates (in m)') +
    scale_color_manual(name = '', values = colours)

# NOTE: CW
dd <- d[!(sub %in% sub_ccw) & !(sub %in% sub_bad) & 
            phase %in% c('Training', 'Generalisation', 'Relearning', 'Washout'), 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), .(phase, target, t_trial)]
plt2 <- ggplot(dd[t_trial < 0.6], aes(V1, V2)) +
    geom_point(aes(x = V1, y = V2, colour='Gaze'), alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4, colour='Trajectories'), alpha = 0.05) +
    geom_point(data = targets, mapping = aes(x = X, y = Y, colour = 'Targets')) +
    xlim(-0.1, 0.1) +
    ylim(0.1, 0.3) +
    facet_wrap(~phase) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(title = 'Subjects with CW Rotation',
         x = 'X coordinates (in m)',
         y = 'Y coordinates (in m)') +
    scale_color_manual(name = '', values = colours)

final <- ggarrange(plt1, plt2, ncol = 2, legend = 'top', common.legend = T)
annotate_figure(final, bottom = 'X coordinates (in m)')
# ggsave('../figures/inspect_difference_rotation/gazetrajdata_all_subs_cw_ccw_relevphases.png', 
#        width=12, height=6.5)


# NOTE: In all phases
colours <- c('Trajectories' = 'blue', 'Gaze' = 'orange', 'Targets' = 'red')

# NOTE: CCW
dd <- d[sub %in% sub_ccw & !(sub %in% sub_bad), 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), .(phase, target, t_trial)]
plt1 <- ggplot(dd[t_trial < 0.6], aes(V1, V2)) +
    geom_point(aes(x = V1, y = V2, colour='Gaze'), alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4, colour='Trajectories'), alpha = 0.05) +
    geom_point(data = targets, mapping = aes(x = X, y = Y, colour = 'Targets')) +
    xlim(-0.1, 0.1) +
    ylim(0.1, 0.3) +
    facet_wrap(~phase) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank()) +
    labs(title = 'Subjects with CCW Rotation',
         x = 'X coordinates (in m)',
         y = 'Y coordinates (in m)') +
    scale_color_manual(name = '', values = colours)

# NOTE: CW
dd <- d[!(sub %in% sub_ccw) & !(sub %in% sub_bad), 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), .(phase, target, t_trial)]
plt2 <- ggplot(dd[t_trial < 0.6], aes(V1, V2)) +
    geom_point(aes(x = V1, y = V2, colour='Gaze'), alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4, colour='Trajectories'), alpha = 0.05) +
    geom_point(data = targets, mapping = aes(x = X, y = Y, colour = 'Targets')) +
    xlim(-0.1, 0.1) +
    ylim(0.1, 0.3) +
    facet_wrap(~phase) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(title = 'Subjects with CW Rotation',
         x = 'X coordinates (in m)',
         y = 'Y coordinates (in m)') +
    scale_color_manual(name = '', values = colours)

final <- ggarrange(plt1, plt2, ncol = 2, legend = 'top', common.legend = T)
annotate_figure(final, bottom = 'X coordinates (in m)')
# ggsave('../figures/inspect_difference_rotation/gazetrajdata_all_subs_cw_ccw.png', 
#        width=12, height=6.5)


# Comparison explicit vs. implicit
# --------------------------------
# NOTE: In relevant phases only
colours <- c('Trajectories' = 'blue', 'Gaze' = 'orange', 'Targets' = 'red')

# NOTE: explicit, ccw
dd <- d[sub %in% sub_expl & sub %in% sub_ccw & !(sub %in% sub_bad) & 
            phase %in% c('Training', 'Generalisation', 'Relearning', 'Washout'), 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), .(phase, target, t_trial)]

plt1 <- ggplot(dd[t_trial < 0.6], aes(V1, V2)) +
    geom_point(aes(x = V1, y = V2, colour='Gaze'), alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4, colour='Trajectories'), alpha = 0.05) +
    geom_point(data = targets, mapping = aes(x = X, y = Y, colour = 'Targets')) +
    xlim(-0.1, 0.1) +
    ylim(0.1, 0.3) +
    facet_wrap(~phase) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(title = 'Explicit, CCW',
         x = 'X coordinates (in m)',
         y = 'Y coordinates (in m)') +
    scale_color_manual(name = '', values = colours)

# NOTE: explicit, cw
dd <- d[sub %in% sub_expl & !(sub %in% sub_ccw) & !(sub %in% sub_bad) & 
            phase %in% c('Training', 'Generalisation', 'Relearning', 'Washout'), 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), .(phase, target, t_trial)]

plt2 <- ggplot(dd[t_trial < 0.6], aes(V1, V2)) +
    geom_point(aes(x = V1, y = V2, colour='Gaze'), alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4, colour='Trajectories'), alpha = 0.05) +
    geom_point(data = targets, mapping = aes(x = X, y = Y, colour = 'Targets')) +
    xlim(-0.1, 0.1) +
    ylim(0.1, 0.3) +
    facet_wrap(~phase) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(title = 'Explicit, CW',
         x = 'X coordinates (in m)',
         y = 'Y coordinates (in m)') +
    scale_color_manual(name = '', values = colours)

# NOTE: implicit, ccw
dd <- d[!(sub %in% sub_expl) & sub %in% sub_ccw & !(sub %in% sub_bad) & 
            phase %in% c('Training', 'Generalisation', 'Relearning', 'Washout'), 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), .(phase, target, t_trial)]
plt3 <- ggplot(dd[t_trial < 0.6], aes(V1, V2)) +
    geom_point(aes(x = V1, y = V2, colour='Gaze'), alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4, colour='Trajectories'), alpha = 0.05) +
    geom_point(data = targets, mapping = aes(x = X, y = Y, colour = 'Targets')) +
    xlim(-0.1, 0.1) +
    ylim(0.1, 0.3) +
    facet_wrap(~phase) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = 'Implicit, CCW',
         x = 'X coordinates (in m)',
         y = 'Y coordinates (in m)') +
    scale_color_manual(name = '', values = colours)

# NOTE: implicit, cw
dd <- d[!(sub %in% sub_expl) & !(sub %in% sub_ccw) & !(sub %in% sub_bad) & 
            phase %in% c('Training', 'Generalisation', 'Relearning', 'Washout'), 
        .(mean(x.ga), mean(y.ga), mean(x.tr), mean(y.tr)), .(phase, target, t_trial)]
plt4 <- ggplot(dd[t_trial < 0.6], aes(V1, V2)) +
    geom_point(aes(x = V1, y = V2, colour='Gaze'), alpha = 0.05) + 
    geom_point(aes(x = V3, y = V4, colour='Trajectories'), alpha = 0.05) +
    geom_point(data = targets, mapping = aes(x = X, y = Y, colour = 'Targets')) +
    xlim(-0.1, 0.1) +
    ylim(0.1, 0.3) +
    facet_wrap(~phase) +
    theme_classic() +
    theme(strip.background = element_rect(colour = 'white', fill = 'white'),
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(title = 'Implicit, CW',
         x = 'X coordinates (in m)',
         y = 'Y coordinates (in m)') +
    scale_color_manual(name = '', values = colours)

final <- ggarrange(plt1, plt2, plt3, plt4, ncol = 2, nrow = 2, legend = 'top', common.legend = T)
annotate_figure(final, bottom = 'X coordinates (in m)', left = 'Y coordinates (in m)')
# ggsave('../figures/inspect_difference_rotation/gazetrajdata_all_subs_aware_relevphases.png', 
#        width=10, height=10)

