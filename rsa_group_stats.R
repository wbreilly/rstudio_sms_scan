# RSA sum stats
library(tidyverse)
library(afex)

# SVSS without positions and 1.1 threshold
rsa =read.table("/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/singletrial_4_rsatoolbox/RSAmeans_SVSS_1.1thresh_11_5_17.txt", header = TRUE, stringsAsFactors = TRUE)

# back to pearsons but stricted beta scrub threshold. some subjects have closer to 15 excluded although some still only a couple
# rsa =read.table("/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/singletrial_4_rsatoolbox/RSAmeans_SVSS_p1p2p3p4p5_pears_1.1betathresh_11_4_17.txt", header = TRUE, stringsAsFactors = TRUE)
# SVSS all ps are separately coded
# rsa =read.table("/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/singletrial_4_rsatoolbox/RSAmeans_SVSS_p1p2p3p4p5_11_4_17.txt", header = TRUE, stringsAsFactors = TRUE)
# rsa$position = as.factor(rsa$position)
# SVSS p4 and p5
# rsa =read.table("/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/singletrial_4_rsatoolbox/RSAmeans_SVSS_p4p5_11_4_17.txt", header = TRUE, stringsAsFactors = TRUE)
# SVSS p1 and p2
# rsa =read.table("/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/singletrial_4_rsatoolbox/RSAmeans_SVSS_p1p2_11_4_17.txt", header = TRUE, stringsAsFactors = TRUE)
# only p1 SVSS
# rsa =read.table("/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/singletrial_4_rsatoolbox/RSAmeans_SVSS_p1_11_3_17.txt", header = TRUE, stringsAsFactors = TRUE)
# all positions SVSS
# rsa =read.table("/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/singletrial_4_rsatoolbox/RSAmeans_SVSS_11_3_17.txt", header = TRUE, stringsAsFactors = TRUE)

rsa = rsa %>% filter(sub != "s008")
rsa.sumstats = rsa %>% group_by(roi, condition) %>% summarise(mean = mean(similarity), sd = sd(similarity))

rsa.fixed = rsa %>% filter(condition != "random")
rsafixed.sumstats = rsa.fixed %>% group_by(roi, condition) %>% summarise(mean = mean(similarity), sd = sd(similarity))

# all.rois = ggplot(data = rsa, aes(x = condition, y = similarity)) + geom_boxplot() + facet_wrap( ~roi, ncol=6)
# all.rois


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
# rsa[,6:9] = NULL
rsa[,5:8] = NULL
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

# plot by cluster
clusters.rois = ggplot(data = rsa.PMAT, aes(x = factor(condition), y = similarity))  + geom_boxplot() + geom_point(aes(colour = factor(sub))) +  facet_wrap( ~roi_group) + scale_y_continuous(limits = c(-.2, .6))
clusters.rois

##########################
# stats
m1 = aov_ez("sub", "similarity", rsa.PMAT, within = c("condition", "roi_group"))
m1
# not surprising

# looking for main effect so only intact and scrambled
rsa.main = rsa %>% filter(condition != "random")
# # # split df's by roi_group
# rsa.PM = rsa.main %>% filter(roi_group == "PM")
# rsa.AT = rsa.main %>% filter(roi_group == "AT")
# rsa.HIPP = rsa.main %>% filter(roi_group == "HIPP")
# look at all conds
rsa.PM = rsa %>% filter(roi_group == "PM")
rsa.AT = rsa %>% filter(roi_group == "AT")
rsa.HIPP = rsa %>% filter(roi_group == "HIPP")

# look at all conditions in HIPP
# RSA.HIPP.allconds = rsa %>% filter(roi_group == "HIPP")
# pairwise ttests
# sample code
# pairwise.t.test(rsa.AT$similarity,rsa.PM$condition, p.adjust.method = "bonferroni")

#  anova.. but this is doing same as paired ttest 
m1.PM = aov_ez("sub", "similarity", rsa.PM, within = c("condition", "roi"))
m1.PM
m1.AT = aov_ez("sub", "similarity", rsa.AT, within = c("condition", "roi"))
m1.AT
m1.HIPP = aov_ez("sub", "similarity", rsa.HIPP, within = c("condition","roi"))
m1.HIPP


###########################################
# Plot PM by region
# plot by cluster + geom_line(group =1)
PM.plot = ggplot(data = rsa.PM, aes(x = factor(condition), y = similarity)) + geom_boxplot() + geom_point(aes(colour = factor(sub)))  + facet_wrap(~roi,ncol = 5)
PM.plot

# Plot PM by region
# plot by cluster + geom_line(group =1)
AT.plot = ggplot(data = rsa.AT, aes(x = factor(condition), y = similarity)) + geom_boxplot() + geom_point(aes(colour = factor(sub)))  + facet_wrap(~roi, ncol = 2)
AT.plot

# Plot PM by region
# plot by cluster + geom_line(group =1)
HIPP.plot = ggplot(data = rsa.HIPP, aes(x = factor(condition), y = similarity)) + geom_boxplot() + geom_point(aes(colour = factor(sub)))  + facet_wrap(~roi, ncol = 2)
HIPP.plot

##########################
# throwing RANDOM back in the mix
# HIPP.plot.allconds = ggplot(data = RSA.HIPP.allconds, aes(x = factor(condition), y = similarity)) + geom_boxplot() + geom_point(aes(colour = factor(sub)))  + facet_wrap(~roi, ncol = 2)
# HIPP.plot.allconds

# look for somethin in hipp
#body only
rsa.HIPP.body = rsa.HIPP %>%  filter(roi != "lh-hipp-head", roi != "rh-hipp-head")
# rsa.HIPP.body.main = rsa.HIPP.body %>% filter(condition != "random")

m1.HIPP = aov_ez("sub", "similarity", rsa.HIPP.body, within = c("condition"))
m1.HIPP
ls.HIPP = lsmeans(m1.HIPP,~condition,contr = "pairwise", adjust = NULL)
ls.HIPP
contrast(ls.HIPP, alpha=0.05, method="pairwise", adjust= "holm")
# intact and random differ in hipp body, scrambled is near-trending

# look in left AT 
# left  only 
rsa.AT.left = rsa.AT %>% filter(roi != "rh-prc", roi != "rh-TPole")
m1.AT = aov_ez("sub", "similarity", rsa.AT, within = c("condition"))
m1.AT
ls.AT = lsmeans(m1.AT,~condition,contr = "pairwise", adjust = NULL)
ls.AT
contrast(ls.AT, alpha=0.05, method="pairwise", adjust= "holm")

# PM regions
rsa.PM.m = rsa.PM %>% filter(roi == "lh-Prec") #c("lh-phc-ant","rh-phc-ant")) # c("lh-Prec" ,"rh-Prec"))#!= "rh-ANG",roi != "lh-ANG",roi != "rh-phc-ant", roi != "lh-phc-ant")
m1.PM = aov_ez("sub", "similarity", rsa.PM.m, within = c("condition"))
m1.PM
ls.PM = lsmeans(m1.PM,~condition,contr = "pairwise", adjust = NULL)
contrast(ls.PM, alpha=0.05, method="pairwise", adjust= "bonferroni")

########
# using to check other script
# stat.check = rsa.PM %>%  group_by( condition,position) %>% summarise(meansim = mean(similarity))
# it passed

#########################################


rsa.diff.ps = rsa %>% group_by(sub,condition,roi_group) %>% filter(roi_group == "PM", roi == "rh-phc-ant", condition != "random") %>% summarise(simmean = mean(similarity))
# rsa.diff.ps = rsa %>% group_by(sub,condition,roi_group) %>% filter(roi_group != "Other", condition != "scrambled") %>%  summarise(simmean = mean(similarity))
rsa.diff.ps = rsa.diff.ps %>%  spread(condition, simmean) 
rsa.diff.ps = rsa.diff.ps %>% mutate(psdiffs = intact - scrambled)
rsa.diff.ps = data.frame(rsa.diff.ps)
# rsa.diff.VMPFC = rsa.diff.ps

# split by roi group
rsa.diff.PM = rsa.diff.ps %>% filter(roi_group == "PM")
rsa.diff.AT = rsa.diff.ps %>% filter(roi_group == "AT")
rsa.diff.HIPP = rsa.diff.ps %>% filter(roi_group == "HIPP")
rsa.diff.other = rsa.diff.ps %>% filter(roi_group == "Other")

# could also grab the data from the retrieval script in rtdiffs$rtdiff
# rsa.diff.PM$rtdiffs = c(251.423047 ,  1.467518 , 93.914570 , 85.780642 , 53.892886 ,201.064728 , 92.368439, 153.167665,  62.682078  ,98.558166 , 31.053571,  -9.381026)
# rsa.diff.AT$rtdiffs =  c(251.423047 ,  1.467518 , 93.914570 , 85.780642 , 53.892886 ,201.064728 , 92.368439, 153.167665,  62.682078  ,98.558166 , 31.053571,  -9.381026)
# rsa.diff.HIPP$rtdiffs = c(251.423047 ,  1.467518 , 93.914570 , 85.780642 , 53.892886 ,201.064728 , 92.368439, 153.167665,  62.682078  ,98.558166 , 31.053571,  -9.381026)
rsa.diff.PM$rtdiffs = rt.diffs$rtdiff
rsa.diff.HIPP$rtdiffs = rt.diffs$rtdiff
rsa.diff.AT$rtdiffs = rt.diffs$rtdiff
rsa.diff.other$rtdiffs = rt.diffs$rtdiff

cor.test(rsa.diff.PM$rtdiffs,rsa.diff.PM$psdiffs)
cor.test(rsa.diff.AT$rtdiffs,rsa.diff.AT$psdiffs)
cor.test(rsa.diff.HIPP$rtdiffs,rsa.diff.HIPP$psdiffs)
cor.test(rsa.diff.other$rtdiffs,rsa.diff.other$psdiffs)

# positive rt means faster on intact than scrambled 

corr.plot = ggplot(data = rsa.diff.PM,(aes(x = psdiffs, y= rtdiffs)))  + geom_point(aes(colour = factor(sub))) +   geom_smooth(method='lm')
corr.plot = corr.plot + ggtitle("r = 0.62, p-value = 0.03")
corr.plot
