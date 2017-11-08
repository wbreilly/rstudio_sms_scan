# 5x5 mtx stuff
library(afex)
library(tidyverse)
library(viridis)
# read in 
rsa =read.table("/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/singletrial_4_rsatoolbox/RSAmeans_DVSS_firstgo_11_4_17.txt", header = TRUE, stringsAsFactors = TRUE)

# filter out bad subs
# s016 has one bad run, s008 has 3. 
rsa = rsa %>% filter(sub != "s008")
# make mtx fam a factor 

rsa$mtx = as.factor(rsa$mtx)
# code each mt square to its condition
rsa$mtx.group = ifelse(rsa$mtx == 1, "diag", ifelse(rsa$mtx == 6, "diag", ifelse(rsa$mtx == 10, "diag", ifelse(rsa$mtx == 13, "diag", ifelse(rsa$mtx == 15, "diag","off")))))
rsa$mtx.group = as.factor(rsa$mtx.group)

# all plot
# all.mtx = ggplot(data = rsa, aes(x = condition, y = similarity)) + geom_boxplot() + facet_wrap( ~roi, ncol=6)
# all.mtx
# not sure what that it supposed to show

# the rois
PM_rois = c('near-phc-ant-L', 'near-phc-ant-R', 'lh-ANG', 'rh-ANG', 'lh-PCC', 'rh-PCC',	'lh-Prec',	'rh-Prec', 'lh-RSC', 'rh-RSC')
AT_rois = c('near-prc-L',	'near-prc-R', 	'lh-TPole', 	'rh-TPole'	)
HIPP_rois = c('near-hipp-body-L',	'near-hipp-body-R',	'near-hipp-head-L',	'near-hipp-head-R')
PM_HIPP_rois = c('near-hipp-body-L',	'near-hipp-body-R',	'near-hipp-head-L',	'near-hipp-head-R', 
                 'near-phc-ant-L', 'near-phc-ant-R', 'lh-ANG', 'rh-ANG', 'lh-PCC', 'rh-PCC',	'lh-Prec',	'rh-Prec', 'lh-RSC', 'rh-RSC')
vmpfc = c('near-VMPFC') 

# make networks factor
# first make logical vectors for member regions
rsa$PM = rsa$roi %in% PM_rois
rsa$AT = rsa$roi %in% AT_rois
rsa$HIPP = rsa$roi %in% HIPP_rois
rsa$PM_HIPP = rsa$roi %in% PM_HIPP_rois

# now assign group label in single column
rsa$roi_group = ifelse(rsa$PM == TRUE, "PM", ifelse(rsa$AT == TRUE, "AT", ifelse(rsa$HIPP == TRUE,"HIPP", "Other" )))
# delete the logical vectors
rsa[,7:10] = NULL
# make group a factor
rsa$roi_group = as.factor(rsa$roi_group)

# wish I did this earlier but hindsight yada yada
rsa$roi = gsub(".*body-L$", "lh-hipp-body", rsa$roi)
rsa$roi = gsub(".*body-R$", "rh-hipp-body", rsa$roi)
rsa$roi = gsub(".*head-L$", "lh-hipp-head", rsa$roi)
rsa$roi = gsub(".*head-R$", "rh-hipp-head", rsa$roi)
rsa$roi = gsub(".*ant-L$", "lh-phc-ant", rsa$roi)
rsa$roi = gsub(".*ant-R$", "rh-phc-ant", rsa$roi)
rsa$roi = gsub(".*VMPFC", "VMPFC", rsa$roi)
rsa$roi = gsub(".*prc-L$", "lh-prc", rsa$roi)
rsa$roi = gsub(".*prc-R$", "rh-prc", rsa$roi)
########################
# PMAT AT and HIPP only
rsa.PMAT = rsa %>% filter(roi_group != "Other")

# plot by diag # + geom_point(aes(colour = factor(sub))) 
clusters.rois = ggplot(data = rsa.PMAT, aes(x = factor(condition), y = similarity))  + geom_boxplot() +  facet_wrap( mtx.group~roi_group) + scale_y_continuous(limits = c(-.2, .6))
clusters.rois

############################################

# first look
m1 = aov_ez("sub", "similarity", rsa.PMAT, within = c("condition", "roi_group", "mtx.group"))
m1

# don't forget about VMPFC
# rsa.VMPFC = rsa %>% filter(roi == "VMPFC")
# m1 = aov_ez("sub", "similarity", rsa.VMPFC, within = c("condition"))
# m1
# nada with condition
rsa.occ = rsa %>% filter(roi == "rh-occ-pole")

# looking for main effect so only intact and scrambled
rsa.main = rsa %>% filter(condition != "random")
# # split df's by roi_group and only look at intact and scrambled
# rsa.PM = rsa.main %>% filter(roi_group == "PM")
# rsa.AT = rsa.main %>% filter(roi_group == "AT")
# rsa.HIPP = rsa.main %>% filter(roi_group == "HIPP")
# look at all conds
rsa.PM = rsa %>% filter(roi_group == "PM")
rsa.AT = rsa %>% filter(roi_group == "AT")
rsa.HIPP = rsa %>% filter(roi_group == "HIPP")
#
# look at all conditions in HIPP
RSA.HIPP.allconds = rsa.PMAT %>% filter(roi_group == "HIPP")
# pairwise ttests
# sample code
# pairwise.t.test(rsa.AT$similarity,rsa.PM$condition, p.adjust.method = "bonferroni")

#  anova.. but this is doing same as paired ttest 
m1.PM = aov_ez("sub", "similarity", rsa.PM, within = c("condition", "mtx.group", "roi"))
m1.PM
m1.AT = aov_ez("sub", "similarity", rsa.AT, within = c("condition","mtx.group","roi"))# "roi"))
m1.AT
m1.HIPP = aov_ez("sub", "similarity", rsa.HIPP, within = c("condition","mtx.group" ,"roi"))#"roi"))
m1.HIPP

###############
# make some 5x5s
# split by condition and pm group
#PM intact
rsa.PM.intact = rsa.PM %>%  filter(condition == "intact")
mtx.stats.pm = rsa.PM.intact %>% group_by(condition,mtx, roi) %>% summarise(meansim = mean(similarity))
mtx.stats.pm = data.frame(mtx.stats.pm)
# mtx.stats.pm$x = c(1,1,1,1,1,2,2,2,2,3,3,3,4,4,5)
# mtx.stats.pm$y = c(5,4,3,2,1,4,3,2,1,3,2,1,2,1,1)
mtx.stats.pm$x = c(rep(1,10),rep(1,10),rep(1,10),rep(1,10),rep(1,10),rep(2,10),rep(2,10),rep(2,10),rep(2,10),rep(3,10),rep(3,10),rep(3,10),rep(4,10),rep(4,10),rep(5,10))
mtx.stats.pm$y = c(rep(5,10),rep(4,10),rep(3,10),rep(2,10),rep(1,10),rep(4,10),rep(3,10),rep(2,10),rep(1,10),rep(3,10),rep(2,10),rep(1,10),rep(2,10),rep(1,10),rep(1,10))

# mtx.stats.pm = mtx.stats.pm %>%  filter(roi == "lh-Prec" | "rh-Prec")

mtx.PM.intact = ggplot(mtx.stats.pm, aes(x = x, y = y, fill = meansim)) + 
  geom_tile() + scale_fill_viridis(na.value = "transparent",limits = (c(0,.35))) + facet_wrap(~ roi,ncol = 5) + theme(axis.text.y = element_blank(), axis.text.x = element_blank())
mtx.PM.intact 


# random
rsa.PM.random = rsa.PM %>%  filter(condition == "random")
mtx.stats.pm = rsa.PM.random %>% group_by(condition,mtx,roi) %>% summarise(meansim = mean(similarity))
mtx.stats.pm = data.frame(mtx.stats.pm)
# mtx.stats.pm$x = c(1,1,1,1,1,2,2,2,2,3,3,3,4,4,5)
# mtx.stats.pm$y = c(5,4,3,2,1,4,3,2,1,3,2,1,2,1,1)
# 
# mtx.PM.intact = ggplot(mtx.stats.pm, aes(x = x, y = y, fill = meansim)) + 
#   geom_tile() + scale_fill_viridis(na.value = "transparent",limits = (c(0,.2))) #+ facet_wrap(~ day) 
# mtx.PM.intact  
mtx.stats.pm$x = c(rep(1,10),rep(1,10),rep(1,10),rep(1,10),rep(1,10),rep(2,10),rep(2,10),rep(2,10),rep(2,10),rep(3,10),rep(3,10),rep(3,10),rep(4,10),rep(4,10),rep(5,10))
mtx.stats.pm$y = c(rep(5,10),rep(4,10),rep(3,10),rep(2,10),rep(1,10),rep(4,10),rep(3,10),rep(2,10),rep(1,10),rep(3,10),rep(2,10),rep(1,10),rep(2,10),rep(1,10),rep(1,10))

mtx.PM.intact = ggplot(mtx.stats.pm, aes(x = x, y = y, fill = meansim)) + 
  geom_tile() + scale_fill_viridis(na.value = "transparent",limits = (c(0,.35))) + facet_wrap(~ roi,ncol = 5) + theme(axis.text.y = element_blank(), axis.text.x = element_blank())
mtx.PM.intact 

# scrambled
rsa.PM.scrambled = rsa.PM %>%  filter(condition == "scrambled")
mtx.stats.pm = rsa.PM.scrambled %>% group_by(condition,mtx,roi) %>% summarise(meansim = mean(similarity))
mtx.stats.pm = data.frame(mtx.stats.pm)
# mtx.stats.pm$x = c(1,1,1,1,1,2,2,2,2,3,3,3,4,4,5)
# mtx.stats.pm$y = c(5,4,3,2,1,4,3,2,1,3,2,1,2,1,1)
# 
# mtx.PM.intact = ggplot(mtx.stats.pm, aes(x = x, y = y, fill = meansim)) + 
#   geom_tile() + scale_fill_viridis(na.value = "transparent",limits = (c(0,.4))) #+ facet_wrap(~ day) 
# mtx.PM.intact  
mtx.stats.pm$x = c(rep(1,10),rep(1,10),rep(1,10),rep(1,10),rep(1,10),rep(2,10),rep(2,10),rep(2,10),rep(2,10),rep(3,10),rep(3,10),rep(3,10),rep(4,10),rep(4,10),rep(5,10))
mtx.stats.pm$y = c(rep(5,10),rep(4,10),rep(3,10),rep(2,10),rep(1,10),rep(4,10),rep(3,10),rep(2,10),rep(1,10),rep(3,10),rep(2,10),rep(1,10),rep(2,10),rep(1,10),rep(1,10))

mtx.PM.intact = ggplot(mtx.stats.pm, aes(x = x, y = y, fill = meansim)) + 
  geom_tile() + scale_fill_viridis(na.value = "transparent",limits = (c(0,.35))) + facet_wrap(~ roi,ncol = 5) + theme(axis.text.y = element_blank(), axis.text.x = element_blank())
mtx.PM.intact 

####################
# AT 
#AT intact
rsa.AT.intact = rsa.AT %>%  filter(condition == "intact")
mtx.stats.pm = rsa.AT.intact %>% group_by(condition,mtx,roi) %>% summarise(meansim = mean(similarity))
mtx.stats.pm = data.frame(mtx.stats.pm)
mtx.stats.pm$x = c(rep(1,4),rep(1,4),rep(1,4),rep(1,4),rep(1,4),rep(2,4),rep(2,4),rep(2,4),rep(2,4),rep(3,4),rep(3,4),rep(3,4),rep(4,4),rep(4,4),rep(5,4))
mtx.stats.pm$y = c(rep(5,4),rep(4,4),rep(3,4),rep(2,4),rep(1,4),rep(4,4),rep(3,4),rep(2,4),rep(1,4),rep(3,4),rep(2,4),rep(1,4),rep(2,4),rep(1,4),rep(1,4))

mtx.PM.intact = ggplot(mtx.stats.pm, aes(x = x, y = y, fill = meansim)) + 
  geom_tile() + scale_fill_viridis(na.value = "transparent",limits = (c(-.05,.15))) + facet_wrap(~ roi,ncol =2) + theme(axis.text.y = element_blank(), axis.text.x = element_blank())
mtx.PM.intact 
# 

# random
rsa.AT.random = rsa.AT %>%  filter(condition == "random")
mtx.stats.pm = rsa.AT.random %>% group_by(condition,mtx,roi) %>% summarise(meansim = mean(similarity))
mtx.stats.pm = data.frame(mtx.stats.pm)
# mtx.stats.pm$x = c(1,1,1,1,1,2,2,2,2,3,3,3,4,4,5)
# mtx.stats.pm$y = c(5,4,3,2,1,4,3,2,1,3,2,1,2,1,1)
# 
# mtx.PM.intact = ggplot(mtx.stats.pm, aes(x = x, y = y, fill = meansim)) + 
#   geom_tile() + scale_fill_viridis(na.value = "transparent",limits = (c(-.02,.1))) #+ facet_wrap(~ day) 
# mtx.PM.intact  
mtx.stats.pm$x = c(rep(1,4),rep(1,4),rep(1,4),rep(1,4),rep(1,4),rep(2,4),rep(2,4),rep(2,4),rep(2,4),rep(3,4),rep(3,4),rep(3,4),rep(4,4),rep(4,4),rep(5,4))
mtx.stats.pm$y = c(rep(5,4),rep(4,4),rep(3,4),rep(2,4),rep(1,4),rep(4,4),rep(3,4),rep(2,4),rep(1,4),rep(3,4),rep(2,4),rep(1,4),rep(2,4),rep(1,4),rep(1,4))

mtx.PM.intact = ggplot(mtx.stats.pm, aes(x = x, y = y, fill = meansim)) + 
  geom_tile() + scale_fill_viridis(na.value = "transparent",limits = (c(-.05,.15))) + facet_wrap(~ roi,ncol =2) + theme(axis.text.y = element_blank(), axis.text.x = element_blank())
mtx.PM.intact 

# scrambled
rsa.AT.scrambled = rsa.AT %>%  filter(condition == "scrambled")
mtx.stats.pm = rsa.AT.scrambled %>% group_by(condition,mtx,roi) %>% summarise(meansim = mean(similarity))
mtx.stats.pm = data.frame(mtx.stats.pm)
# mtx.stats.pm$x = c(1,1,1,1,1,2,2,2,2,3,3,3,4,4,5)
# mtx.stats.pm$y = c(5,4,3,2,1,4,3,2,1,3,2,1,2,1,1)
# 
# mtx.PM.intact = ggplot(mtx.stats.pm, aes(x = x, y = y, fill = meansim)) + 
#   geom_tile() + scale_fill_viridis(na.value = "transparent",limits = (c(-.02,.1))) #+ facet_wrap(~ day) 
# mtx.PM.intact  
mtx.stats.pm$x = c(rep(1,4),rep(1,4),rep(1,4),rep(1,4),rep(1,4),rep(2,4),rep(2,4),rep(2,4),rep(2,4),rep(3,4),rep(3,4),rep(3,4),rep(4,4),rep(4,4),rep(5,4))
mtx.stats.pm$y = c(rep(5,4),rep(4,4),rep(3,4),rep(2,4),rep(1,4),rep(4,4),rep(3,4),rep(2,4),rep(1,4),rep(3,4),rep(2,4),rep(1,4),rep(2,4),rep(1,4),rep(1,4))

mtx.PM.intact = ggplot(mtx.stats.pm, aes(x = x, y = y, fill = meansim)) + 
  geom_tile() + scale_fill_viridis(na.value = "transparent",limits = (c(-.05,.15))) + facet_wrap(~ roi,ncol =2) + theme(axis.text.y = element_blank(), axis.text.x = element_blank())
mtx.PM.intact 

####################
# HIPP
#PM intact
# rsa.HIPP.body = RSA.HIPP.allconds %>%  filter(roi != "lh-hipp-head", roi != "rh-hipp-head")
rsa.HIPP.intact = rsa.HIPP %>%  filter(condition == "intact")
mtx.stats.pm = rsa.HIPP.intact %>% group_by(condition,mtx,roi) %>% summarise(meansim = mean(similarity))
mtx.stats.pm = data.frame(mtx.stats.pm)
# mtx.stats.pm$x = c(1,1,1,1,1,2,2,2,2,3,3,3,4,4,5)
# mtx.stats.pm$y = c(5,4,3,2,1,4,3,2,1,3,2,1,2,1,1)
# 
# mtx.PM.intact = ggplot(mtx.stats.pm, aes(x = x, y = y, fill = meansim)) + 
#   geom_tile() + scale_fill_viridis(na.value = "transparent",limits = (c(-.02,.06))) #+ facet_wrap(~ day) 
# mtx.PM.intact 
mtx.stats.pm$x = c(rep(1,4),rep(1,4),rep(1,4),rep(1,4),rep(1,4),rep(2,4),rep(2,4),rep(2,4),rep(2,4),rep(3,4),rep(3,4),rep(3,4),rep(4,4),rep(4,4),rep(5,4))
mtx.stats.pm$y = c(rep(5,4),rep(4,4),rep(3,4),rep(2,4),rep(1,4),rep(4,4),rep(3,4),rep(2,4),rep(1,4),rep(3,4),rep(2,4),rep(1,4),rep(2,4),rep(1,4),rep(1,4))

mtx.PM.intact = ggplot(mtx.stats.pm, aes(x = x, y = y, fill = meansim)) + 
  geom_tile() + scale_fill_viridis(na.value = "transparent",limits = (c(-.04,.06))) + facet_wrap(~ roi,ncol =2) + theme(axis.text.y = element_blank(), axis.text.x = element_blank())
mtx.PM.intact 
# 

# random
# rsa.HIPP.random = rsa.HIPP.body %>%  filter(condition == "random")
rsa.HIPP.random = rsa.HIPP %>%  filter(condition == "random")
mtx.stats.pm = rsa.HIPP.random %>% group_by(condition,mtx,roi) %>% summarise(meansim = mean(similarity))
mtx.stats.pm = data.frame(mtx.stats.pm)
# mtx.stats.pm$x = c(1,1,1,1,1,2,2,2,2,3,3,3,4,4,5)
# mtx.stats.pm$y = c(5,4,3,2,1,4,3,2,1,3,2,1,2,1,1)
# 
# mtx.PM.intact = ggplot(mtx.stats.pm, aes(x = x, y = y, fill = meansim)) + 
#   geom_tile() + scale_fill_viridis(na.value = "transparent",limits = (c(-.02,.06))) #+ facet_wrap(~ day) 
# mtx.PM.intact  
mtx.stats.pm$x = c(rep(1,4),rep(1,4),rep(1,4),rep(1,4),rep(1,4),rep(2,4),rep(2,4),rep(2,4),rep(2,4),rep(3,4),rep(3,4),rep(3,4),rep(4,4),rep(4,4),rep(5,4))
mtx.stats.pm$y = c(rep(5,4),rep(4,4),rep(3,4),rep(2,4),rep(1,4),rep(4,4),rep(3,4),rep(2,4),rep(1,4),rep(3,4),rep(2,4),rep(1,4),rep(2,4),rep(1,4),rep(1,4))

mtx.PM.intact = ggplot(mtx.stats.pm, aes(x = x, y = y, fill = meansim)) + 
  geom_tile() + scale_fill_viridis(na.value = "transparent",limits = (c(-.04,.06))) + facet_wrap(~ roi,ncol =2) + theme(axis.text.y = element_blank(), axis.text.x = element_blank())
mtx.PM.intact 

# scrambled
# rsa.HIPP.scrambled = rsa.HIPP.body %>%  filter(condition == "scrambled")
rsa.HIPP.scrambled = rsa.HIPP %>%  filter(condition == "scrambled")
mtx.stats.pm = rsa.HIPP.scrambled %>% group_by(condition,mtx,roi) %>% summarise(meansim = mean(similarity))
mtx.stats.pm = data.frame(mtx.stats.pm)
# mtx.stats.pm$x = c(1,1,1,1,1,2,2,2,2,3,3,3,4,4,5)
# mtx.stats.pm$y = c(5,4,3,2,1,4,3,2,1,3,2,1,2,1,1)
# 
# mtx.PM.intact = ggplot(mtx.stats.pm, aes(x = x, y = y, fill = meansim)) + 
#   geom_tile() + scale_fill_viridis(na.value = "transparent",limits = (c(-.02,.05))) #+ facet_wrap(~ day) 
# mtx.PM.intact  
mtx.stats.pm$x = c(rep(1,4),rep(1,4),rep(1,4),rep(1,4),rep(1,4),rep(2,4),rep(2,4),rep(2,4),rep(2,4),rep(3,4),rep(3,4),rep(3,4),rep(4,4),rep(4,4),rep(5,4))
mtx.stats.pm$y = c(rep(5,4),rep(4,4),rep(3,4),rep(2,4),rep(1,4),rep(4,4),rep(3,4),rep(2,4),rep(1,4),rep(3,4),rep(2,4),rep(1,4),rep(2,4),rep(1,4),rep(1,4))

mtx.PM.intact = ggplot(mtx.stats.pm, aes(x = x, y = y, fill = meansim)) + 
  geom_tile() + scale_fill_viridis(na.value = "transparent",limits = (c(-.04,.06))) + facet_wrap(~ roi,ncol =2) + theme(axis.text.y = element_blank(), axis.text.x = element_blank())
mtx.PM.intact 

##############################################
#  VMPFC comin to da r
rsa.VMPFC= rsa %>%  filter(roi == "rh-occ-pole", condition == "intact")
mtx.stats.pm = rsa.VMPFC %>% group_by(condition,mtx) %>% summarise(meansim = mean(similarity))
mtx.stats.pm = data.frame(mtx.stats.pm)
mtx.stats.pm$x = c(1,1,1,1,1,2,2,2,2,3,3,3,4,4,5)
mtx.stats.pm$y = c(5,4,3,2,1,4,3,2,1,3,2,1,2,1,1)

mtx.PM.intact = ggplot(mtx.stats.pm, aes(x = x, y = y, fill = meansim)) +
  geom_tile() + scale_fill_viridis(na.value = "transparent",limits = (c(0,.3))) #+ fcet_wrap(~ day)
mtx.PM.intact

# random
rsa.VMPFC= rsa %>%  filter(roi == "rh-occ-pole", condition == "random")
mtx.stats.pm = rsa.VMPFC %>% group_by(condition,mtx) %>% summarise(meansim = mean(similarity))
mtx.stats.pm = data.frame(mtx.stats.pm)
mtx.stats.pm$x = c(1,1,1,1,1,2,2,2,2,3,3,3,4,4,5)
mtx.stats.pm$y = c(5,4,3,2,1,4,3,2,1,3,2,1,2,1,1)

mtx.PM.intact = ggplot(mtx.stats.pm, aes(x = x, y = y, fill = meansim)) +
  geom_tile() + scale_fill_viridis(na.value = "transparent",limits = (c(0,.3))) #+ facet_wrap(~ day)
mtx.PM.intact

#scrambled
rsa.VMPFC= rsa %>%  filter(roi == "lh-occ-pole", condition == "scrambled")
mtx.stats.pm = rsa.occ %>% group_by(condition,mtx) %>% summarise(meansim = mean(similarity))
mtx.stats.pm = data.frame(mtx.stats.pm)
mtx.stats.pm$x = c(1,1,1,1,1,2,2,2,2,3,3,3,4,4,5)
mtx.stats.pm$y = c(5,4,3,2,1,4,3,2,1,3,2,1,2,1,1)

mtx.PM.intact = ggplot(mtx.stats.pm, aes(x = x, y = y, fill = meansim)) +
  geom_tile() + scale_fill_viridis(na.value = "transparent",limits = (c(0,.3))) + theme(axis.text.y = element_blank(), axis.text.x = element_blank()) #+ facet_wrap(~ day)
mtx.PM.intact

