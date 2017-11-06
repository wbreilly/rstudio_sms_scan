# RSA sum stats
library(tidyverse)
library(afex)

# back to pearsons but stricted beta scrub threshold. some subjects have closer to 15 excluded although some still only a couple
rsa =read.table("/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/singletrial_4_rsatoolbox/RSAmeans_SVSS_p1p2p3p4p5_pears_1.1betathresh_11_4_17.txt", header = TRUE, stringsAsFactors = TRUE)
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
# # split df's by roi_group
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
m1.PM = aov_ez("sub", "similarity", rsa.PM, within = c("condition", "roi"))
m1.PM
m1.AT = aov_ez("sub", "similarity", rsa.AT, within = c("condition", "roi"))
m1.AT
m1.HIPP = aov_ez("sub", "similarity", rsa.HIPP, within = c("condition","roi"))
m1.HIPP

ls.HIPP = lsmeans(m1.HIPP,~ roi:condition,contr = "pairwise", adjust = NULL)
ls1
pairs(ls1,adjust = "bonferroni")
c1 = list(c(1, 0, 0, 0, -1, 0, 0, 0))
contrast(ls1,c1, adjust = "bonferroni")
c2 = list(c(0, 1, 0, 0, 0, -1, 0, 0))
contrast(ls1,c1)

ls.AT = lsmeans(m1.AT,~ roi:condition,contr = "pairwise", adjust = NULL)
ls.AT
# left t pole contrast
c1 = list(c(1, 0, 1, 0, -1, 0, -1, 0))
contrast(ls.AT,c1)
# right t pole contrast
c1 = list(c(0, 0, 0, 1, 0, 0, 0,-1))
contrast(ls.AT,c1)

ls.PM = lsmeans(m1.PM,~ roi:condition,contr = "pairwise", adjust = NULL)
ls.PM
# scrambled > intact at every region
c1 = list(all = c(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1,  1,  1,  1,  1,  1,  1,  1,  1,  1 ),
          prec = c(0, 0, 0, 0, 1/2, 0, 0, 0, 1/2, 0, 0, 0, 0, 0, -1/2, 0, 0, 0, -1/2, 0))
contrast(ls.PM,c1) #,adjust = "bonferroni")
# right t pole contrast
c1 = list(c(0, 0, 0, 1, 0, 0, 0,-1))
contrast(ls.AT,c1)

###########################################
# Plot PM by region
# plot by cluster + geom_line(group =1)
PM.plot = ggplot(data = rsa.PM, aes(x = factor(position), y = similarity)) + geom_boxplot() + geom_point(aes(colour = factor(sub)))  + facet_wrap(roi ~condition)
PM.plot

# Plot PM by region
# plot by cluster + geom_line(group =1)
AT.plot = ggplot(data = rsa.AT, aes(x = factor(condition), y = similarity)) + geom_boxplot() + geom_point(aes(colour = factor(sub)))  + facet_wrap(~roi, ncol = 2)
AT.plot

# Plot PM by region
# plot by cluster + geom_line(group =1)
HIPP.plot = ggplot(data = rsa.HIPP, aes(x = factor(position), y = similarity)) + geom_boxplot() + geom_point(aes(colour = factor(sub)))  + facet_wrap(roi~condition)
HIPP.plot

##########################
# throwing RANDOM back in the mix
HIPP.plot.allconds = ggplot(data = RSA.HIPP.allconds, aes(x = factor(condition), y = similarity)) + geom_boxplot() + geom_point(aes(colour = factor(sub)))  + facet_wrap(~roi, ncol = 2)
HIPP.plot.allconds

# look for somethin in hipp
#body only
rsa.HIPP.body = RSA.HIPP.allconds %>%  filter(roi != "lh-hipp-head", roi != "rh-hipp-head")
rsa.HIPP.body.main = rsa.HIPP.body %>% filter(condition != "random")
m1.HIPP = aov_ez("sub", "similarity", rsa.HIPP.body, within = c("condition", "roi"))
m1.HIPP
ls.HIPP = lsmeans(m1.HIPP,~condition,contr = "pairwise", adjust = NULL)
ls.HIPP
contrast(ls.HIPP, alpha=0.05, method="pairwise", adjust= "bonferroni")
# intact and random differ in hipp body, scrambled is near-trending

# look in left AT 
# left  only 
rsa.AT.left = rsa.AT %>% filter(roi != "rh-prc", roi != "rh-TPole")
m1.AT = aov_ez("sub", "similarity", rsa.AT, within = c("condition"))
m1.AT
ls.AT = lsmeans(m1.AT,~condition,contr = "pairwise", adjust = NULL)
ls.AT
contrast(ls.AT, alpha=0.05, method="pairwise", adjust= "bonferroni")

# PM regions
rsa.PM.m = rsa.PM %>% filter(roi == "lh-Prec") #c("lh-phc-ant","rh-phc-ant")) # c("lh-Prec" ,"rh-Prec"))#!= "rh-ANG",roi != "lh-ANG",roi != "rh-phc-ant", roi != "lh-phc-ant")
m1.PM = aov_ez("sub", "similarity", rsa.PM.m, within = c("condition"))
m1.PM
ls.PM = lsmeans(m1.PM,~condition,contr = "pairwise", adjust = NULL)
contrast(ls.PM, alpha=0.05, method="pairwise", adjust= "bonferroni")

################################
# using to check other script
stat.check = rsa.PM %>%  group_by( condition,position) %>% summarise(meansim = mean(similarity))
