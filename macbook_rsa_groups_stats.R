# RSA sum stats
library(tidyverse)
library(afex)


# same verb RSA even for random
rsa.SVDP = read.table("/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/singletrial_4_rsatoolbox/november_17/RSAmeans_FS_ANTS_SVDPforrandom_allps_pears_1.15betathresh_3_3_18.txt", header = TRUE, stringsAsFactors = TRUE)
# change label for random verb based
rsa.SVDP$condition = ifelse(rsa.SVDP$condition == "random", "r_sameverb",ifelse(rsa.SVDP$condition == "intact","intact","scrambled")) 
rsa.SVDP = rsa.SVDP %>% filter(condition == "r_sameverb")
# SVSS with positions and beta threshold
rsa.SVSS = read.table("/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/singletrial_4_rsatoolbox/november_17/RSAmeans_FS_ANTS_SSsamepos_allps_pears_1.15betathresh_3_3_18.txt", header = TRUE, stringsAsFactors = TRUE)

# add random r_sameverb for one combined dataset
rsa = bind_rows(rsa.SVDP, rsa.SVSS)


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
rsa[,6:9] = NULL
# rsa[,5:8] = NULL
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

# plot by cluster # + geom_point(aes(colour = factor(sub)))
clusters.rois = ggplot(data = rsa.PMAT, aes(x = factor(condition), y = similarity))  + geom_boxplot(outlier.shape = NA)  +  facet_wrap( ~roi_group) + scale_y_continuous(limits = c(-.15, .5))
clusters.rois = clusters.rois + xlab("Condition") + theme(text = element_text(size = 30, face = "bold"))+ theme(strip.background = element_blank(),strip.text= element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())
clusters.rois
# conddition label order is the following: intact r_sameverb random scrambled
# network order is: AT HIPP PM
# ggsave("networkwide_condition_effects.pdf", plot = last_plot(), dpi = 600 )

##########################
# stats
m1 = aov_ez("sub", "similarity", rsa.PMAT, within = c("condition", "roi_group"))
m1
ls.m1 = lsmeans(m1,~condition|roi_group,contr = "pairwise", adjust = NULL)
ls.m1
contrast(ls.m1, alpha=0.05, method="pairwise", adjust= "holm")


# looking for main effect so only intact and scrambled
# rsa.main = rsa %>% filter(condition != "random")
# # # split df's by roi_group
# rsa.PM = rsa.main %>% filter(roi_group == "PM")
# rsa.AT = rsa.main %>% filter(roi_group == "AT")
# rsa.HIPP = rsa.main %>% filter(roi_group == "HIPP")


# look at all conds
rsa.PM = rsa %>% filter(roi_group == "PM")
rsa.AT = rsa %>% filter(roi_group == "AT")
rsa.HIPP = rsa %>% filter(roi_group == "HIPP")

# plot different conditions by network group
# PM 
PM.all = ggplot(data = rsa.PM, aes(x = factor(condition), y = similarity))  + geom_boxplot(outlier.shape = NA) + scale_y_continuous(limits = c(-.15, .5))
PM.all = PM.all + ylab("Similarity") + theme(text = element_text(size = 30, face = "bold"),axis.title.x=element_blank(),axis.text.x=element_text(angle = -45,hjust = 0)) + scale_x_discrete(labels = c("Intact", "Same Verb","Same Position","Scrambled"))
PM.all
# ggsave("pm.all.pdf", plot = last_plot(), dpi = 600 )

# AT
AT.all = ggplot(data = rsa.AT, aes(x = factor(condition), y = similarity))  + geom_boxplot(outlier.shape = NA) + scale_y_continuous(limits = c(-.15, .4))
AT.all = AT.all + ylab("Similarity") + theme(text = element_text(size = 30, face = "bold"),axis.title.x=element_blank(),axis.text.x=element_text(angle = -45,hjust = 0)) + scale_x_discrete(labels = c("Intact", "Same Verb","Same Position","Scrambled"))
AT.all
# ggsave("at.all.pdf", plot = last_plot(), dpi = 600 )

# HIPP
HIPP.all = ggplot(data = rsa.HIPP, aes(x = factor(condition), y = similarity))  + geom_boxplot(outlier.shape = NA) + scale_y_continuous(limits = c(-.1, .25))
HIPP.all = HIPP.all +ylab("Similarity") + theme(text = element_text(size = 30, face = "bold"),axis.title.x=element_blank(),axis.text.x=element_text(angle = -45,hjust = 0)) + scale_x_discrete(labels = c("Intact", "Same Verb","Same Position","Scrambled"))
HIPP.all
# ggsave("hipp.all.pdf", plot = last_plot(), dpi = 600 )

#  look which ROIs in network drive effects
m1.PM = aov_ez("sub", "similarity", rsa.PM, within = c("condition", "roi"))
m1.PM
ls.PM = lsmeans(m1.PM,~condition|roi,contr = "pairwise", adjust = NULL)
ls.PM
contrast(ls.PM, alpha=0.05, method="pairwise", adjust= "holm")

m1.AT = aov_ez("sub", "similarity", rsa.AT, within = c("condition", "roi"))
m1.AT
ls.AT = lsmeans(m1.AT,~condition|roi,contr = "pairwise", adjust = NULL)
ls.AT
contrast(ls.AT, alpha=0.05, method="pairwise", adjust= "holm")

m1.HIPP = aov_ez("sub", "similarity", rsa.HIPP, within = c("condition","roi"))
m1.HIPP
ls.HIPP = lsmeans(m1.HIPP,~condition|roi,contr = "pairwise", adjust = NULL)
ls.HIPP
contrast(ls.HIPP, alpha=0.05, method="pairwise", adjust= "holm")

rsa.VMPFC = rsa %>% filter(roi == "VMPFC")
m1.VMPFC = aov_ez("sub", "similarity", rsa.VMPFC, within = c("condition"))
m1.VMPFC
ls.VMPFC = lsmeans(m1.VMPFC,~condition,contr = "pairwise", adjust = NULL)
ls.VMPFC
contrast(ls.VMPFC, alpha=0.05, method="pairwise", adjust= "holm")

###########################################
# Plot PM by region
# plot by cluster + geom_line(group =1) #  + geom_point(aes(colour = factor(sub))) 
PM.plot = ggplot(data = rsa.PM, aes(x = factor(condition), y = similarity)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~roi,ncol = 5) + theme(text = element_text(size = 25, face = "bold"),axis.title.x=element_blank(),axis.text.x=element_text(angle = -90,hjust = 0))
PM.plot = PM.plot + ylab("Similarity") + scale_x_discrete(labels = c("Intact", "Same Verb","Same Position","Scrambled"))
PM.plot
# ggsave("pm.regiona.pdf", plot = last_plot(), dpi = 600 )


# Plot AT by region
# plot by cluster + geom_line(group =1) # + geom_point(aes(colour = factor(sub)))  
rsa.AT = rsa %>% filter(roi_group == "AT")
# rsa.AT = rsa.AT %>% filter(roi != "rh-TPole", roi != "rh-prc")
# val = c(1,2,3)
# rsa.AT$condition = factor(rsa.AT$condition, levels = rsa.AT$condition[order(val)])
# levels(rsa.AT$condition)
# rsa.AT$condition = as.factor(rsa.AT$condtion)
# AT.plot = ggplot(data = rsa.AT, aes(x = factor(condition), y = similarity, fill = condition)) + geom_boxplot() + facet_wrap(~roi, ncol = 3)
# AT.plot = AT.plot + coord_cartesian(ylim=c(-.25,.35)) + theme(text = element_text(size = 30, face = "bold")) + guides(fill = FALSE)
# AT.plot = AT.plot + labs(x = NULL) + scale_fill_brewer(palette= "Dark2")  + theme(strip.background = element_blank(),strip.text= element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())
# AT.plot
AT.plot = ggplot(data = rsa.AT, aes(x = factor(condition), y = similarity)) +coord_cartesian(ylim=c(-.3,.4))+ geom_boxplot(outlier.shape = NA) + facet_wrap(~roi) + theme(text = element_text(size = 25, face = "bold"),axis.title.x=element_blank(),axis.text.x=element_text(angle = -90,hjust = 0))
AT.plot =  AT.plot + ylab("Similarity") + scale_x_discrete(labels = c("Intact", "Same Verb","Same Position","Scrambled"))
AT.plot
# ggsave("at.regiona.pdf", plot = last_plot(), dpi = 600 )

# Plot HIPP by region
# plot by cluster + geom_line(group =1) # + geom_point(aes(colour = factor(sub))) 
rsa.HIPP = rsa %>% filter(roi_group == "HIPP")
# rsa.HIPP = rsa %>% filter(roi_group == "HIPP", roi != "lh-hipp-head", roi != "rh-hipp-head")
# val = c(1,2,3)
# rsa.HIPP$condition = factor(rsa.HIPP$condition, levels = rsa.HIPP$condition[order(val)])
# levels(rsa.HIPP$condition)
# rsa.HIPP$condition = as.factor(rsa.HIPP$condtion)
# HIPP.plot = ggplot(data = rsa.HIPP, aes(x = factor(condition), y = similarity, fill = condition)) + geom_boxplot()  + facet_wrap(~roi, ncol = 2)
# HIPP.plot = HIPP.plot + coord_cartesian(ylim=c(-.25,.25)) + theme(text = element_text(size = 30, face = "bold")) + guides(fill = FALSE)
# HIPP.plot = HIPP.plot + labs(x = NULL) + scale_fill_brewer(palette= "Dark2")  + theme(strip.background = element_blank(),strip.text= element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())
# HIPP.plot
# ggsave("hipp_body_SVDP.eps", plot = last_plot(), dpi = 600 )
HIPP.plot = ggplot(data = rsa.HIPP, aes(x = factor(condition), y = similarity)) +coord_cartesian(ylim=c(-.22,.27))+ geom_boxplot(outlier.shape = NA) + facet_wrap(~roi)  + theme(text = element_text(size = 25, face = "bold"),axis.title.x=element_blank(),axis.text.x=element_text(angle = -90,hjust = 0))
HIPP.plot = HIPP.plot +  ylab("Similarity") + scale_x_discrete(labels = c("Intact", "Same Verb","Same Position","Scrambled"))
HIPP.plot
# ggsave("hipp.region.pdf", plot = last_plot(), dpi = 600 )
##########################
#

########
# using to check other script
# stat.check = rsa.PM %>%  group_by( condition,position) %>% summarise(meansim = mean(similarity))
# it passed

#########################################
# Diff correlatins for intact and scrambled
# no behavior data for sub 23 
rsa.diff.ps = rsa %>% filter(sub != "s023")
# roi == "lh-phc-ant",
rsa.diff.ps = rsa.diff.ps %>% group_by(sub,condition,roi_group) %>% filter(roi_group == "PM", roi == "rh-Prec", condition != "scrambled", condition != "random") %>% summarise(simmean = mean(similarity))
# rsa.diff.ps = rsa.diff.ps %>% group_by(sub,condition,roi_group) %>% filter(roi_group != "Other", condition != "random") %>%  summarise(simmean = mean(similarity))
rsa.diff.ps = rsa.diff.ps %>%  spread(condition, simmean) 
rsa.diff.ps = rsa.diff.ps %>% mutate(psdiffs = intact - r_sameverb)
rsa.diff.ps = data.frame(rsa.diff.ps)
# rsa.diff.VMPFC = rsa.diff.ps

# split by roi group
rsa.diff.PM = rsa.diff.ps %>% filter(roi_group == "PM")
# rsa.diff.AT = rsa.diff.ps %>% filter(roi_group == "AT")
# rsa.diff.HIPP = rsa.diff.ps %>% filter(roi_group == "HIPP")
# rsa.diff.other = rsa.diff.ps %>% filter(roi_group == "Other")

# grab diffs
rsa.diff.PM$rtdiffs = rt.diffs.rand$rtdiff

# get rid of outlier
# rsa.diff.PM = rsa.diff.PM %>% filter(sub != "s007")

cor.test(rsa.diff.PM$rtdiffs,rsa.diff.PM$psdiffs)
# cor.test(rsa.diff.AT$rtdiffs,rsa.diff.AT$psdiffs)
# cor.test(rsa.diff.HIPP$rtdiffs,rsa.diff.HIPP$psdiffs)
# cor.test(rsa.diff.other$rtdiffs,rsa.diff.other$psdiffs)

# positive rt means faster on intact than scrambled 
# positive ps means higher ps for intact 

corr.plot = ggplot(data = rsa.diff.PM,(aes(x = psdiffs, y= rtdiffs)))  + geom_point(size = 5) +   geom_smooth(method='lm')
corr.plot = corr.plot + ggtitle("r = -.02, p = 0.94")  + theme(legend.position="none") + labs(x = "Pattern Similarity Difference", y = "Reaction Time Difference (ms)")
corr.plot = corr.plot + theme(text = element_text(size = 30, face = "bold"))
corr.plot.phc = corr.plot
corr.plot.phc 
# ggsave("rh_prec_corr_intact_sameverb.pdf", plot = last_plot(), dpi = 600 )

####################
# AT
# no behavior data for sub 23 
rsa.diff.ps = rsa %>% filter(sub != "s023")
#,  roi == "lh-TPole"
rsa.diff.ps = rsa.diff.ps %>% group_by(sub,condition,roi_group) %>% filter(roi_group == "AT",  roi == "rh-prc",condition != "random",condition != "scrambled") %>% summarise(simmean = mean(similarity))
# rsa.diff.ps = rsa %>% group_by(sub,condition,roi_group) %>% filter(roi_group != "Other", condition != "scrambled") %>%  summarise(simmean = mean(similarity))
rsa.diff.ps = rsa.diff.ps %>%  spread(condition, simmean) 
rsa.diff.ps = rsa.diff.ps %>% mutate(psdiffs = intact - r_sameverb)
rsa.diff.ps = data.frame(rsa.diff.ps)
# rsa.diff.VMPFC = rsa.diff.ps

# split by roi group
# rsa.diff.PM = rsa.diff.ps %>% filter(roi_group == "PM")
rsa.diff.AT = rsa.diff.ps %>% filter(roi_group == "AT")
# rsa.diff.HIPP = rsa.diff.ps %>% filter(roi_group == "HIPP")
# rsa.diff.other = rsa.diff.ps %>% filter(roi_group == "Other")

#get diffs
rsa.diff.AT$rtdiffs = rt.diffs.rand$rtdiff


# get rid of outlier 


# cor.test(rsa.diff.PM$rtdiffs,rsa.diff.PM$psdiffs)
cor.test(rsa.diff.AT$rtdiffs,rsa.diff.AT$psdiffs)
# cor.test(rsa.diff.HIPP$rtdiffs,rsa.diff.HIPP$psdiffs)
# cor.test(rsa.diff.other$rtdiffs,rsa.diff.other$psdiffs)

# plot
corr.plot = ggplot(data = rsa.diff.AT,(aes(x = psdiffs, y= rtdiffs)))  + geom_point(size = 5) +   geom_smooth(method='lm')
corr.plot = corr.plot + ggtitle("r = -.07, p = 0.79")  + theme(legend.position="none") + labs(x = "Pattern Similarity Difference", y = "Reaction Time Difference (ms)")
corr.plot = corr.plot + theme(text = element_text(size = 30, face = "bold"))
corr.plot
ggsave("rh_prc_diff_corr_intact_sameverb.pdf", plot = last_plot(), dpi = 600 )

####################
# HIPP
# no behavior data for sub 23 
# rsa.diff.ps = rsa %>% filter(sub != "s023")
# 
# rsa.diff.ps = rsa.diff.ps %>% group_by(sub,condition,roi_group) %>% filter(roi_group == "HIPP",  roi == "rh-hipp-body",condition != "random", condition != "r_sameverb") %>% summarise(simmean = mean(similarity))
# # rsa.diff.ps = rsa %>% group_by(sub,condition,roi_group) %>% filter(roi_group != "Other", condition != "scrambled") %>%  summarise(simmean = mean(similarity))
# rsa.diff.ps = rsa.diff.ps %>%  spread(condition, simmean)
# rsa.diff.ps = rsa.diff.ps %>% mutate(psdiffs = intact - scrambled)
# rsa.diff.ps = data.frame(rsa.diff.ps)
# # rsa.diff.VMPFC = rsa.diff.ps
# 
# # split by roi group
# # rsa.diff.PM = rsa.diff.ps %>% filter(roi_group == "PM")
# # rsa.diff.AT = rsa.diff.ps %>% filter(roi_group == "AT")
# # rsa.diff.HIPP = rsa.diff.ps %>% filter(roi_group == "HIPP")
# # rsa.diff.other = rsa.diff.ps %>% filter(roi_group == "Other")
# 
# # could also grab the data from the retrieval script in rtdiffs$rtdiff
# # rsa.diff.PM$rtdiffs = c(251.423047 ,  1.467518 , 93.914570 , 85.780642 , 53.892886 ,201.064728 , 92.368439, 153.167665,  62.682078  ,98.558166 , 31.053571,  -9.381026)
# # rsa.diff.AT$rtdiffs =  c(251.423047 ,  1.467518 , 93.914570 , 85.780642 , 53.892886 ,201.064728 , 92.368439, 153.167665,  62.682078  ,98.558166 , 31.053571,  -9.381026)
# # rsa.diff.HIPP$rtdiffs = c(251.423047 ,  1.467518 , 93.914570 , 85.780642 , 53.892886 ,201.064728 , 92.368439, 153.167665,  62.682078  ,98.558166 , 31.053571,  -9.381026)
# # rsa.diff.PM$rtdiffs = rt.diffs$rtdiff
# rsa.diff.HIPP$rtdiffs = rt.diffs$rtdiff
# # rsa.diff.AT$rtdiffs = rt.diffs$rtdiff
# # rsa.diff.other$rtdiffs = rt.diffs$rtdiff
# 
# #try without s018 leverage point
# rsa.diff.HIPP = rsa.diff.HIPP %>% filter(sub != "s018")
# 
# # cor.test(rsa.diff.PM$rtdiffs,rsa.diff.PM$psdiffs)
# # cor.test(rsa.diff.AT$rtdiffs,rsa.diff.AT$psdiffs)
# cor.test(rsa.diff.HIPP$rtdiffs,rsa.diff.HIPP$psdiffs)
# # cor.test(rsa.diff.other$rtdiffs,rsa.diff.other$psdiffs)
# 
# # plot
# corr.plot = ggplot(data = rsa.diff.HIPP,(aes(x = psdiffs, y= rtdiffs)))  + geom_point(size = 5) +   geom_smooth(method='lm')
# corr.plot = corr.plot + ggtitle("r = 0.46, p = 0.07")  + theme(legend.position="none") + labs(x = "Pattern Similarity Difference", y = "Reaction Time Difference (ms)")
# corr.plot = corr.plot + theme(text = element_text(size = 30, face = "bold")) + geom_point(aes(colour = factor(sub)))
# corr.plot.hipp = corr.plot
# corr.plot.hipp
# ggsave("_diff_corr.pdf", plot = last_plot(), dpi = 600 )

##################################################
##################################################
#################################################
# Diff correlatins for intact and random
# no behavior data for sub 23 
rsa.diff.ps = rsa %>% filter(sub != "s023")
# roi == "lh-phc-ant",
rsa.diff.ps = rsa.diff.ps %>% group_by(sub,condition,roi_group) %>% filter(roi_group == "PM", roi == "rh-Prec", condition != "random", condition != "r_sameverb") %>% summarise(simmean = mean(similarity))
# rsa.diff.ps = rsa.diff.ps %>% group_by(sub,condition,roi_group) %>% filter(roi_group != "Other", condition != "random") %>%  summarise(simmean = mean(similarity))
rsa.diff.ps = rsa.diff.ps %>%  spread(condition, simmean) 
rsa.diff.ps = rsa.diff.ps %>% mutate(psdiffs = intact - scrambled)
rsa.diff.ps = data.frame(rsa.diff.ps)
# rsa.diff.VMPFC = rsa.diff.ps

# split by roi group
rsa.diff.PM = rsa.diff.ps %>% filter(roi_group == "PM")
# rsa.diff.AT = rsa.diff.ps %>% filter(roi_group == "AT")
# rsa.diff.HIPP = rsa.diff.ps %>% filter(roi_group == "HIPP")
# rsa.diff.other = rsa.diff.ps %>% filter(roi_group == "Other")

# grab the diffs
rsa.diff.PM$rtdiffs = rt.diffs$rtdiff

# test the corr
cor.test(rsa.diff.PM$rtdiffs,rsa.diff.PM$psdiffs)


# positive rt means faster on intact than scrambled 
# positive ps means higher ps for intact 

corr.plot = ggplot(data = rsa.diff.PM,(aes(x = psdiffs, y= rtdiffs)))  + geom_point(size = 5) +   geom_smooth(method='lm')
corr.plot = corr.plot + ggtitle("r = 0.48, p = 0.06")  + theme(legend.position="none") + labs(x = "Pattern Similarity Difference", y = "Reaction Time Difference (ms)")
corr.plot = corr.plot + theme(text = element_text(size = 30, face = "bold"))
corr.plot.phc = corr.plot
corr.plot.phc 
# ggsave("phc_diff_corr_intact-rsameverb.pdf", plot = last_plot(), dpi = 600 )

####################
# AT ---- nothing happening with PRC but negative correlations with TPole for intact - random and intact - r_sameverb
# no behavior data for sub 23 
rsa.diff.ps = rsa %>% filter(sub != "s023")
#,  roi == "lh-TPole"
rsa.diff.ps = rsa.diff.ps %>% group_by(sub,condition,roi_group) %>% filter(roi_group == "AT",  roi == "rh-prc",condition != "r_sameverb",condition != "scrambled") %>% summarise(simmean = mean(similarity))
# rsa.diff.ps = rsa %>% group_by(sub,condition,roi_group) %>% filter(roi_group != "Other", condition != "scrambled") %>%  summarise(simmean = mean(similarity))
rsa.diff.ps = rsa.diff.ps %>%  spread(condition, simmean) 
rsa.diff.ps = rsa.diff.ps %>% mutate(psdiffs = intact - random)
rsa.diff.ps = data.frame(rsa.diff.ps)
# rsa.diff.VMPFC = rsa.diff.ps

# split by roi group
# rsa.diff.PM = rsa.diff.ps %>% filter(roi_group == "PM")
rsa.diff.AT = rsa.diff.ps %>% filter(roi_group == "AT")
# rsa.diff.HIPP = rsa.diff.ps %>% filter(roi_group == "HIPP")
# rsa.diff.other = rsa.diff.ps %>% filter(roi_group == "Other")

# add rt diffs
rsa.diff.AT$rtdiffs = rt.diffs.rand$rtdiff
# corr test
cor.test(rsa.diff.AT$rtdiffs,rsa.diff.AT$psdiffs)

# plot
corr.plot = ggplot(data = rsa.diff.AT,(aes(x = psdiffs, y= rtdiffs)))  + geom_point(size = 5) +   geom_smooth(method='lm')
corr.plot = corr.plot + ggtitle("r = -0.42, p = 0.10")  + theme(legend.position="none") + labs(x = "Pattern Similarity Difference", y = "Reaction Time Difference (ms)")
corr.plot = corr.plot + theme(text = element_text(size = 30, face = "bold"))
corr.plot
# ggsave("lh_Tpole_diff_corr_intact-random.pdf", plot = last_plot(), dpi = 600 )

####################
# HIPP ---------nothing going on here. Damn. 
# no behavior data for sub 23 
rsa.diff.ps = rsa %>% filter(sub != "s023")


rsa.diff.ps = rsa.diff.ps %>% group_by(sub,condition,roi_group) %>% filter(roi_group == "HIPP",  roi == "rh-hipp-body",condition != "r_sameverb", condition != "scrambled") %>% summarise(simmean = mean(similarity))
# rsa.diff.ps = rsa %>% group_by(sub,condition,roi_group) %>% filter(roi_group != "Other", condition != "scrambled") %>%  summarise(simmean = mean(similarity))
rsa.diff.ps = rsa.diff.ps %>%  spread(condition, simmean)
rsa.diff.ps = rsa.diff.ps %>% mutate(psdiffs = intact - random)
rsa.diff.ps = data.frame(rsa.diff.ps)
# rsa.diff.VMPFC = rsa.diff.ps

# split by roi group
# rsa.diff.PM = rsa.diff.ps %>% filter(roi_group == "PM")
# rsa.diff.AT = rsa.diff.ps %>% filter(roi_group == "AT")
# rsa.diff.HIPP = rsa.diff.ps %>% filter(roi_group == "HIPP")
# rsa.diff.other = rsa.diff.ps %>% filter(roi_group == "Other")

# could also grab the data from the retrieval script in rtdiffs$rtdiff
# rsa.diff.PM$rtdiffs = c(251.423047 ,  1.467518 , 93.914570 , 85.780642 , 53.892886 ,201.064728 , 92.368439, 153.167665,  62.682078  ,98.558166 , 31.053571,  -9.381026)
# rsa.diff.AT$rtdiffs =  c(251.423047 ,  1.467518 , 93.914570 , 85.780642 , 53.892886 ,201.064728 , 92.368439, 153.167665,  62.682078  ,98.558166 , 31.053571,  -9.381026)
# rsa.diff.HIPP$rtdiffs = c(251.423047 ,  1.467518 , 93.914570 , 85.780642 , 53.892886 ,201.064728 , 92.368439, 153.167665,  62.682078  ,98.558166 , 31.053571,  -9.381026)
# rsa.diff.PM$rtdiffs = rt.diffs$rtdiff

# s018 hipp ROIs are missing...? look into
rt.diffs.rand = rt.diffs.rand %>% filter(sub != "18")
rsa.diff.HIPP$rtdiffs = rt.diffs.rand$rtdiff
# rsa.diff.AT$rtdiffs = rt.diffs$rtdiff
# rsa.diff.other$rtdiffs = rt.diffs$rtdiff

#try without s018 leverage point
# rsa.diff.HIPP = rsa.diff.HIPP %>% filter(sub != "s018")

# cor.test(rsa.diff.PM$rtdiffs,rsa.diff.PM$psdiffs)
# cor.test(rsa.diff.AT$rtdiffs,rsa.diff.AT$psdiffs)
cor.test(rsa.diff.HIPP$rtdiffs,rsa.diff.HIPP$psdiffs)
# cor.test(rsa.diff.other$rtdiffs,rsa.diff.other$psdiffs)

# plot
# corr.plot = ggplot(data = rsa.diff.HIPP,(aes(x = psdiffs, y= rtdiffs)))  + geom_point(size = 5) +   geom_smooth(method='lm')
# corr.plot = corr.plot + ggtitle("r = 0.46, p = 0.07")  + theme(legend.position="none") + labs(x = "Pattern Similarity Difference", y = "Reaction Time Difference (ms)")
# corr.plot = corr.plot + theme(text = element_text(size = 30, face = "bold")) + geom_point(aes(colour = factor(sub)))
# corr.plot.hipp = corr.plot
# corr.plot.hipp
# ggsave("_diff_corr.pdf", plot = last_plot(), dpi = 600 )

#######################################################
# fake data plots
conditions = c("Intact", "Same Verb","Same Position","Scrambled")
intact_vals =   rnorm(n=17, mean=.3, sd=.1)
sameverb_vals =   rnorm(n=17, mean=.05, sd=.1)
samepos_vals =   rnorm(n=17, mean=.15, sd=.1)
scrambled_vals =   rnorm(n=17, mean=.15, sd=.1)
fakedf = data.frame(intact_vals,sameverb_vals,samepos_vals,scrambled_vals)
fakedf = fakedf %>% gather("condition","similarity")


fake.plot = ggplot(data = fakedf, aes(x = factor(condition), y = similarity)) +coord_cartesian(ylim=c(-.1,.4))+ geom_boxplot(outlier.shape = NA)  + theme(text = element_text(size = 25, face = "bold"),axis.title.x=element_blank(),axis.text.x=element_text(angle = -45,hjust = 0))
fake.plot =  fake.plot + ylab("Similarity") + scale_x_discrete(labels = c("Intact", "Same Verb","Same Position","Scrambled"))
fake.plot
ggsave("fake_AT_h2.pdf", plot = last_plot(), dpi = 600 )
# 