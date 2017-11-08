# RSA sum stats
library(tidyverse)
library(afex)

# back to pearsons but stricted beta scrub threshold. some subjects have closer to 15 excluded although some still only a couple
rsa =read.table("/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/singletrial_4_rsatoolbox/RSAmeans_SVSS_p1p2p3p4p5_pears_1.1betathresh_11_4_17.txt", header = TRUE, stringsAsFactors = TRUE)
# spearman RS matrices
# rsa =read.table("/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/singletrial_4_rsatoolbox/RSAmeans_SVSS_p1p2p3p4p5_spearman_11_4_17.txt", header = TRUE, stringsAsFactors = TRUE)
# SVSS all ps are separately coded
# rsa =read.table("/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/singletrial_4_rsatoolbox/RSAmeans_SVSS_p1p2p3p4p5_11_4_17.txt", header = TRUE, stringsAsFactors = TRUE)
rsa$position = as.factor(rsa$position)

# SVSS p4 and p5
# rsa =read.table("/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/singletrial_4_rsatoolbox/RSAmeans_SVSS_p4p5_11_4_17.txt", header = TRUE, stringsAsFactors = TRUE)
# SVSS p1 and p2
# rsa =read.table("/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/singletrial_4_rsatoolbox/RSAmeans_SVSS_p1p2_11_4_17.txt", header = TRUE, stringsAsFactors = TRUE)
# only p1 SVSS
# rsa =read.table("/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/singletrial_4_rsatoolbox/RSAmeans_SVSS_p1_11_3_17.txt", header = TRUE, stringsAsFactors = TRUE)
# all positions SVSS
# rsa =read.table("/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/singletrial_4_rsatoolbox/RSAmeans_SVSS_11_3_17.txt", header = TRUE, stringsAsFactors = TRUE)


# filter out bad subs
rsa = rsa %>% filter(sub != "s008") #sub != "s016")



rsa.sumstats = rsa %>% group_by(roi, condition) %>% summarise(mean = mean(similarity), sd = sd(similarity))

rsa.fixed = rsa %>% filter(condition != "random")
rsafixed.sumstats = rsa.fixed %>% group_by(roi, condition) %>% summarise(mean = mean(similarity), sd = sd(similarity))

all.rois = ggplot(data = rsa, aes(x = condition, y = similarity)) + geom_boxplot() + facet_wrap( ~roi, ncol=6)
all.rois


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
rsa[,6:9] = NULL
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


# looking for outliers
# rsa = rsa %>% filter(position == 5, condition == "scrambled", roi_group == "AT")
########################
# PMAT AT and HIPP only
rsa.PMAT = rsa %>% filter(roi_group != "Other")

# plot by cluster
clusters.rois = ggplot(data = rsa.PMAT, aes(x = factor(condition), y = similarity))  + geom_boxplot() + geom_point(aes(colour = factor(sub))) +  facet_wrap( ~roi_group) + scale_y_continuous(limits = c(-.2, .6))
clusters.rois
clusters.rois = ggplot(data = rsa.PMAT, aes(x = factor(condition), y = similarity))  + geom_boxplot() + facet_wrap( ~roi_group) #+ scale_y_continuous(limits = c(-.2, .6))
clusters.rois

##########################
# stats
# looking for main effect so only intact and scrambled
rsa.main = rsa.PMAT %>% filter(condition != "scrambled")
# with all 3 conditions
m1 = aov_ez("sub", "similarity", rsa.PMAT, within = c("condition", "roi_group"))
m1

# split df's by roi_group
# rsa.PM = rsa.main %>% filter(roi_group == "PM")
# rsa.AT = rsa.main %>% filter(roi_group == "AT")
# rsa.HIPP = rsa.main %>% filter(roi_group == "HIPP")
# split df's by roi_group
# rsa.PM = rsa %>% filter(roi_group == "PM")
# rsa.AT = rsa %>% filter(roi_group == "AT")
# rsa.HIPP = rsa %>% filter(roi_group == "HIPP")
#
# look at all conditions in HIPP
# RSA.HIPP.allconds = rsa.PMAT %>% filter(roi_group == "HIPP")
# pairwise ttests
# sample code
# pairwise.t.test(rsa.AT$similarity,rsa.PM$condition, p.adjust.method = "bonferroni")

#  anova
m1.PM = aov_ez("sub", "similarity", rsa.PM, within = c("condition", "roi","position")) # , "position"))
m1.PM
m1.AT = aov_ez("sub", "similarity", rsa.AT, within = c("condition", "roi", "position"))#, "position"))
m1.AT
m1.HIPP = aov_ez("sub", "similarity", rsa.HIPP, within = c("condition","roi", "position")) #, "position"))
m1.HIPP

#contrasts aren't valid with position in model
ls.HIPP = lsmeans(m1.HIPP,~ roi:condition,contr = "pairwise", adjust = NULL)
ls.HIPP
pairs(ls.HIPP,adjust = "bonferroni")
# c1 = list(c(1, 0, 0, 0, -1, 0, 0, 0))
# contrast(ls.HIPP,c1, adjust = "bonferroni")
# c2 = list(c(0, 1, 0, 0, 0, -1, 0, 0))
# contrast(ls.HIPP,c1)

###################
# AT by position
# look in left AT 
# left  only 
# rsa.AT.left = rsa.AT %>% filter(roi != "rh-prc", roi != "rh-TPole")
rsa.AT.left = rsa.AT %>% filter(roi == "lh-TPole")
m1.AT = aov_ez("sub", "similarity", rsa.AT.left, within = c("condition","position"))
m1.AT
ls.AT = lsmeans(m1.AT,~condition:position,contr = "pairwise", adjust = NULL)
ls.AT
c1 = list(con.p1 = c(1, -1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0),
          con.p2 = c(0, 0 ,1 ,-1 ,0 ,0 ,0 ,0 ,0 , 0),
          con.p3 = c(0,0,0,0,1,-1,0,0,0,0),
          con.p4 = c(0,0,0,0,0,0,1,-1,0,0),
          con.p5 = c(0,0,0,0,0,0,0,0,1,-1))
contrast(ls.AT, c1, adjust= "holm")



#############
# ls.AT = lsmeans(m1.AT,~ roi:condition,contr = "pairwise", adjust = NULL)
# ls.AT
# # # left t pole contrast
# c1 = list(c(1, 0, 1, 0, -1, 0, -1, 0))
# contrast(ls.AT,c1)
# # right t pole contrast
# c1 = list(c(0, 0, 0, 1, 0, 0, 0,-1))
# contrast(ls.AT,c1)
# 
# # lsmeans approach to difference means
# prec by position?
# if only 
rsa.PM = rsa.PM %>% filter(roi == "lh-Prec")
m1.PM = aov_ez("sub", "similarity", rsa.PM, within = c("condition","position" ))
m1.PM
ls.PM = lsmeans(m1.PM,~condition:position,contr = "pairwise", adjust = NULL)
ls.PM
pairs(ls.PM,adjust = "bonferroni")
# do all position between two conditions contrast
c1 = list(con.p1 = c(1, -1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0),
          con.p2 = c(0, 0 ,1 ,-1 ,0 ,0 ,0 ,0 ,0 , 0),
          con.p3 = c(0,0,0,0,1,-1,0,0,0,0),
          con.p4 = c(0,0,0,0,0,0,1,-1,0,0),
          con.p5 = c(0,0,0,0,0,0,0,0,1,-1))

# c1 = list(con.p1p2 = c(1/2, -1/2, 1/2 ,-1/2, 0 ,0, 0, 0, 0, 0 ),
#           con.p1 = c(1,-1,0,0,0,0,0,0,0,0),
#           con.p2 = c(0,-0,1,-1,0,0,0,0,0,0),
#           con.p4p5 = c(0,0,0,0,0,0,1/2,-1/2,1/2,-1/2))
# scrambled > intact at every region
# c1 = list(all = c(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1,  1,  1,  1,  1,  1,  1,  1,  1,  1 ),
#           prec_pcc = c(0, 0, 0, 1/4, 1/4, 0, 0, 1/4, 1/4, 0, 0, 0, 0, -1/4, -1/4, 0, 0, -1/4, -1/4, 0))
contrast(ls.PM,c1 ,adjust = "holm")
# # right t pole contrast
# c1 = list(c(0, 0, 0, 1, 0, 0, 0,-1))
# contrast(ls.AT,c1)

###########################################
# Plot PM by region
# plot by cluster + geom_line(group =1) # + geom_point(aes(colour = factor(sub)))
rsa.PM = rsa.main %>% filter(roi_group == "PM")

PM.plot = ggplot(data = rsa.PM, aes(x = factor(position), y = similarity)) + geom_boxplot()   + facet_wrap(condition~roi) + scale_y_continuous(limits = c(-.2, .6))
PM.plot

# Plot AT by region
# plot by cluster + geom_line(group =1)
AT.plot = ggplot(data = rsa.AT, aes(x = factor(position), y = similarity)) + geom_boxplot() + geom_point(aes(colour = factor(sub)))  + facet_wrap(roi~condition) + scale_y_continuous(limits = c(-.2, .6))
AT.plot

# Plot HIPP by region
# plot by cluster + geom_line(group =1)
HIPP.plot = ggplot(data = rsa.HIPP, aes(x = factor(position), y = similarity)) + geom_boxplot() + geom_point(aes(colour = factor(sub)))  + facet_wrap(~condition, ncol = 2) + scale_y_continuous(limits = c(-.2, .6))
HIPP.plot

##########################
# throwing RANDOM back in the mix
HIPP.plot.allconds = ggplot(data = RSA.HIPP.allconds, aes(x = factor(condition), y = similarity)) + geom_boxplot() + geom_point(aes(colour = factor(sub)))  + facet_wrap(~roi, ncol = 2)
HIPP.plot.allconds


# plot differences between intact and scrambled across positions
# PM regions
sumstats.pos.plot = rsa.PM %>% group_by(sub, condition, position) %>% summarise(simmean = mean(similarity))
sumstats.df = data.frame(sumstats.pos.plot)
sumstats.pos.plot = sumstats.df %>% group_by(condition,position) %>% summarise(mean = mean(simmean), sd = sd(simmean))
sumstats.pos.plot = data.frame(sumstats.pos.plot)
sumstats.pos.plot = mutate(sumstats.pos.plot, SE = sd/sqrt(12))

# create a bar graph
limits <- aes(ymax = sumstats.pos.plot$mean + sumstats.pos.plot$SE,
              ymin = sumstats.pos.plot$mean - sumstats.pos.plot$SE)

p.3way <- ggplot(data = sumstats.pos.plot, aes(x = factor(position), y = mean,
                                             fill = factor(condition)))
p.3way = p.3way + geom_bar(stat = "identity",
                                 position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.15) 
p.3way

# plot differences between intact and scrambled across positions
# AT regions
sumstats.pos.plot = rsa.AT %>% group_by(sub, condition, position) %>% summarise(simmean = mean(similarity))
sumstats.df = data.frame(sumstats.pos.plot)
# filter to look at just p4
# sumstats.df = sumstats.df %>% filter(position ==4)
sumstats.pos.plot = sumstats.df %>% group_by(condition,position) %>% summarise(mean = mean(simmean), sd = sd(simmean))
sumstats.pos.plot = data.frame(sumstats.pos.plot)
sumstats.pos.plot = mutate(sumstats.pos.plot, SE = sd/sqrt(12))

# create a bar graph
limits <- aes(ymax = sumstats.pos.plot$mean + sumstats.pos.plot$SE,
              ymin = sumstats.pos.plot$mean - sumstats.pos.plot$SE)

p.3way <- ggplot(data = sumstats.pos.plot, aes(x = factor(position), y = mean,
                                               fill = factor(condition)))
p.3way = p.3way + geom_bar(stat = "identity",
                           position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.15) 
p.3way

# plot differences between intact and scrambled across positions
# HIPP regions
# actually just in the body
rsa.HIPP.body = rsa.HIPP %>%  filter(roi != "lh-hipp-head", roi != "rh-hipp-head")
sumstats.pos.plot = rsa.HIPP.body %>% group_by(sub, condition, position) %>% summarise(simmean = mean(similarity))
sumstats.df = data.frame(sumstats.pos.plot)
# filter to look at just p4
# sumstats.df = sumstats.df %>% filter(position ==4)
sumstats.pos.plot = sumstats.df %>% group_by(condition,position) %>% summarise(mean = mean(simmean), sd = sd(simmean))
sumstats.pos.plot = data.frame(sumstats.pos.plot)
sumstats.pos.plot = mutate(sumstats.pos.plot, SE = sd/sqrt(12))

# create a bar graph
limits <- aes(ymax = sumstats.pos.plot$mean + sumstats.pos.plot$SE,
              ymin = sumstats.pos.plot$mean - sumstats.pos.plot$SE)

p.3way <- ggplot(data = sumstats.pos.plot, aes(x = factor(position), y = mean,
                                               fill = factor(condition)))
p.3way = p.3way + geom_bar(stat = "identity",
                           position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.15) 
p.3way

#### attempt 2 at difference plots
# plot differences between intact and scrambled across positions
# PM regions
# lsmeans approach to difference means
# m1.PM = aov_ez("sub", "similarity", rsa.PM, within = c("condition", "position"))
# m1.PM
# ls.PM = lsmeans(m1.PM,~ position:condition,contr = "pairwise", adjust = NULL)
# ls.PM = data.frame(summary(ls.PM))
# sumstats.pos.plot = ls.PM %>% spread(condition, lsmean ) 
# sumstats.pos.plot =  sumstats.pos.plot %>% group_by(position) %>% mutate(diff = scrambled - intact)
# diffs = sumstats.pos.plot[is.na(sumstats.pos.plot$intact),]
# diffs = data.frame(diffs)
# diffs$intact = sumstats.pos.plot[!is.na(sumstats.pos.plot$intact),6]
# diffs$diffs =  diffs$scrambled - diffs$intact
# ok i could plot that, but SE is all the same

##########################################
###### gaahhhh why was this so hard to figure out
# PM diffs
sumstats.pos.plot = rsa.PM %>% group_by(sub, condition, position) %>% summarise(simmean = mean(similarity))
sumstats.pos.plot = sumstats.pos.plot %>% spread(condition, simmean) 
sumstats.pos.plot = sumstats.pos.plot %>% mutate(diffs = scrambled - intact)
sumstats.pos.plot = sumstats.pos.plot %>% group_by(position) %>% summarise(diffmean = mean(diffs), sddiff = sd(diffs))
sumstats.pos.plot = sumstats.pos.plot %>% mutate(SE = sddiff/sqrt(12))
#
# create a bar graph
limits = aes(ymax = sumstats.pos.plot$diffmean + sumstats.pos.plot$SE,
              ymin = sumstats.pos.plot$diffmean - sumstats.pos.plot$SE)
diffs.PM = ggplot(data = sumstats.pos.plot, aes(x = factor(position), y = diffmean))
diffs.PM = diffs.PM + geom_bar(stat = "identity",
                           position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.15) 
diffs.PM

# AT difffs
sumstats.pos.plot = rsa.AT %>% group_by(sub, condition, position) %>% summarise(simmean = mean(similarity))
sumstats.pos.plot = sumstats.pos.plot %>% spread(condition, simmean) 
sumstats.pos.plot = sumstats.pos.plot %>% mutate(diffs = scrambled - intact)
sumstats.pos.plot = sumstats.pos.plot %>% group_by(position) %>% summarise(diffmean = mean(diffs), sddiff = sd(diffs))
sumstats.pos.plot = sumstats.pos.plot %>% mutate(SE = sddiff/sqrt(12))
#
# create a bar graph
limits = aes(ymax = sumstats.pos.plot$diffmean + sumstats.pos.plot$SE,
             ymin = sumstats.pos.plot$diffmean - sumstats.pos.plot$SE)
diffs.PM = ggplot(data = sumstats.pos.plot, aes(x = factor(position), y = diffmean))
diffs.PM = diffs.PM + geom_bar(stat = "identity",
                               position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.15) 
diffs.PM

# Hipp diffs
sumstats.pos.plot = rsa.HIPP.body %>% group_by(sub, condition, position) %>% summarise(simmean = mean(similarity))
sumstats.pos.plot = sumstats.pos.plot %>% spread(condition, simmean) 
sumstats.pos.plot = sumstats.pos.plot %>% mutate(diffs = scrambled - intact)
sumstats.pos.plot = sumstats.pos.plot %>% group_by(position) %>% summarise(diffmean = mean(diffs), sddiff = sd(diffs))
sumstats.pos.plot = sumstats.pos.plot %>% mutate(SE = sddiff/sqrt(12))
#
# create a bar graph
limits = aes(ymax = sumstats.pos.plot$diffmean + sumstats.pos.plot$SE,
             ymin = sumstats.pos.plot$diffmean - sumstats.pos.plot$SE)
diffs.PM = ggplot(data = sumstats.pos.plot, aes(x = factor(position), y = diffmean))
diffs.PM = diffs.PM + geom_bar(stat = "identity",
                               position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.15) 
diffs.PM

# All in one plot 
sumstats.pos.plot = rsa.main %>% group_by(sub, condition, position, roi_group) %>% summarise(simmean = mean(similarity))
sumstats.pos.plot = sumstats.pos.plot %>% spread(condition, simmean) 
sumstats.pos.plot = sumstats.pos.plot %>% mutate(diffs = scrambled - intact)
sumstats.pos.plot = sumstats.pos.plot %>% group_by(position,roi_group) %>% summarise(diffmean = mean(diffs), sddiff = sd(diffs))
sumstats.pos.plot = sumstats.pos.plot %>% mutate(SE = sddiff/sqrt(12))
#
# create a bar graph
limits = aes(ymax = sumstats.pos.plot$diffmean + sumstats.pos.plot$SE,
             ymin = sumstats.pos.plot$diffmean - sumstats.pos.plot$SE)
diffs.PM = ggplot(data = sumstats.pos.plot, aes(x = factor(position), y = diffmean, fill = factor(roi_group)))
diffs.PM = diffs.PM + geom_bar(stat = "identity",
                               position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.15) 
diffs.PM
