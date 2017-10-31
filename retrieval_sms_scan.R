# Walter Reilly
# sms2 analysis script
# started 8.8.17
# copied from sms 2 retrieval script

library(tidyverse)

d = read.csv("~/drive/grad_school/DML_WBR/Sequences_Exp3/sms_scan_drive/sms_scan_fmri_copy/all_sms_scan_retreival_dat/all_retrieval_10_30_17.csv", 
             na.strings="NaN", stringsAsFactors=FALSE)

d$rt = as.numeric(d$rt)
d$sub = as.numeric(d$sub)

#get rid of column info for each block
d = d[!(d$block == "" | d$block == "block"), ]


# add position info
d = mutate(d, position = rep(c(1:5),times = dim(d)[1]/5))

# exclude for same consecutive responses (3sds) or too many "no_resp"

# crunch the no_resp's
response.counts = count(group_by(d,sub,response))
noresp.exclude = response.counts %>% filter(response == "no_resp") 
x = mean(noresp.exclude$n)
y = sd(noresp.exclude$n)
noresp.exclude = noresp.exclude %>% mutate(zresp = (n - x)/y)
# looks like it isn't an issue so far. only 10 total from 4 subjects. 
# Might want to look at proportion of yes/no at some point

###################################################
# ID participants based on how frequently they give the same response 
# x =
# d$responsecopy = d$response
# x= 0 
# for (i in 2:dim(d)[1]){
#   if (d[i,1] == d[i-1,1]){
#     if (d[i,13] == d[i-1,13]){
#       d[i,14] = 1 + x; x = x+1
#     } else {
#       x = 0; d[i,14] = 0
#     }
#   } else {
#     x = 0; d[i,14] = 0
#   }
# }
# 
# 
# # count how many consecutive same responses
# nonono = d %>% group_by(sub) %>% summarise(mean = mean(V14), sd = sd(V14))
# nonono = nonono %>% mutate(zmean = scale(mean))


# add observation column to df
d = mutate(d, obs = 1:dim(d)[1])

# this didn't work for sms2 retr dat. NaN instead of NA?
exrows = d[!complete.cases(d$rt),14]
dsansnas = d[-c(exrows),]


# randoms need unified condition label
# create new column with numeric condition labels
# intact = 1, scrambled = 2, random = 3
dsansnas = dsansnas %>% mutate(con.num = condition)
dsansnas[,15] = ifelse(dsansnas[,15] == "intact", 1, ifelse(dsansnas[,15] == "scrambled", 2, 3))

# trim rts
sub.con.mean.rt = dsansnas %>%  group_by(sub,con.num) %>% 
  summarise(mean = mean(rt, na.rm = TRUE), SD = sd(rt, na.rm = TRUE))

dsansnas = dsansnas %>%  group_by(sub,con.num) %>% 
  mutate(mean.sub.con = mean(rt, na.rm = TRUE), SD.sub.con = sd(rt, na.rm = TRUE))

dsansnas = dsansnas %>% mutate(sd3 = 3 * SD.sub.con) 


# exlude too long
# 32 for the first 4 subjects. Not too bad
dsansnas = dsansnas %>% mutate(exclude.high = rt > mean.sub.con + sd3)
# exclude too short
# Zero, per usual
dsansnas = dsansnas %>% mutate(exclude.low = rt < mean.sub.con - sd3)
# exlude hella fast <  100ms
# only 3. Nice!
dsansnas = dsansnas %>% mutate(exclude.hellafast = rt < 100)

####put in new df
dclean = dsansnas %>% filter(exclude.low == FALSE, exclude.high == FALSE, exclude.hellafast == FALSE)

# sum stats
sumstats = dclean %>% group_by(con.num) %>% summarise(mean = mean(rt, na.rm = TRUE), SD = sd(rt, na.rm = TRUE))
sumstats.rep = dclean %>% group_by(repetition) %>% summarise(mean = mean(rt, na.rm = TRUE), SD = sd(rt, na.rm = TRUE))
sumstats.rep.con = dclean %>% group_by(repetition,con.num) %>% summarise(mean = mean(rt, na.rm = TRUE), SD = sd(rt, na.rm = TRUE))



# look at rt by position
#####################################################################################
sumstats.pos = dclean %>% group_by(con.num, position) %>% summarise(mean = mean(rt), SD = sd(rt))

p.pos =  ggplot(data=sumstats.pos, aes(x=position, y=mean, group=con.num, colour = con.num)) +
  geom_line(stat = "identity") + 
  ggtitle("SMS2: Condition * Position ")
p.pos

#############################################
# mean pos 4-5 by condition
#exclude all 1 positions

d25 = dclean %>% filter(!position == 1)

################################
# make plot of mean differences
d25sumstats.con = d25 %>% group_by(con.num) %>% summarise(mean = mean(rt), SD = sd(rt))

# trial counts
d25.counts = count(group_by(d25,con.num))

# add n to sumstats
d25sumstats.con[,4] = d25.counts[,2]

# add SE
d25sumstats.con = mutate(d25sumstats.con, SE = SD/sqrt(n))

# create a bar graph
limits <- aes(ymax = d25sumstats.con$mean + d25sumstats.con$SE,
              ymin = d25sumstats.con$mean - d25sumstats.con$SE)

p.con <- ggplot(data = d25sumstats.con, aes(x = factor(con.num), y = mean))
p.con = p.con + geom_bar(stat = "identity",position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.15) + 
  labs(x = "Conditions", y = "Pos 2-5 Mean RT") +
  ggtitle("Mean RT of Positions 2-5") +
  scale_x_discrete("Conditions", labels = 
                     c("1" = "Intact", "2" = "Scrambled-Fixed", 
                       "3" = "Scrambled-Random"))
p.con = p.con + coord_cartesian(ylim=c(500,1000))
p.con




# here it goes
# get some summary statistics
d25sumstats = d25 %>% group_by(con.num) %>% summarise(mean = mean(rt), SD = sd(rt))
d25sumstats.rep = d25 %>% group_by(repetition) %>% summarise(mean = mean(rt), SD = sd(rt))
d25sumstats.rep.con = d25 %>% group_by(repetition,con.num) %>% summarise(mean = mean(rt), SD = sd(rt))

# trial counts
d25.counts = count(d25, con.num)


# repetition by condition and only pos's 2-5
# get n
counts2 = count(group_by(d25,con.num, repetition))

# add n to sumstats
d25sumstats.rep.con[,5] = counts2[,3]

# add SE
d25sumstats.rep.con = mutate(d25sumstats.rep.con, SE = SD/sqrt(n))

# create a bar graph
limits <- aes(ymax = d25sumstats.rep.con$mean + d25sumstats.rep.con$SE,
              ymin = d25sumstats.rep.con$mean - d25sumstats.rep.con$SE)

p.rep.con <- ggplot(data = d25sumstats.rep.con, aes(x = factor(con.num), y = mean,
                                                    fill = factor(repetition)))
p.rep.con = p.rep.con + geom_bar(stat = "identity",position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.15) + 
  labs(x = "Conditions", y = "Pos 2-5 Mean RT") +
  ggtitle("Mean RT by Condition and Repetition") +
  scale_x_discrete("Conditions", labels = 
                     c("1" = "Intact", "2" = "Scrambled-Fixed", 
                       "3" = "Scrambled-Random"))
p.rep.con = p.rep.con + coord_cartesian(ylim=c(500,1200))
p.rep.con

############################
# condition and rrb
d25sumstats.rrb.con = d25 %>% group_by(rrb,con.num) %>% summarise(mean = mean(rt), SD = sd(rt))

# trial counts
d25.counts = count(d25, con.num)
# repetition by condition and only pos's 2-5
# get n
d25$rrb = as.factor(d25$rrb)
counts2 = count(group_by(d25,con.num, rrb))

# add n to sumstats
d25sumstats.rrb.con[,5] = counts2[,3]

# add SE
d25sumstats.rrb.con = mutate(d25sumstats.rrb.con, SE = SD/sqrt(n))

# create a bar graph
limits <- aes(ymax = d25sumstats.rrb.con$mean + d25sumstats.rrb.con$SE,
              ymin = d25sumstats.rrb.con$mean - d25sumstats.rrb.con$SE)

p.rrb.con <- ggplot(data = d25sumstats.rrb.con, aes(x = factor(con.num), y = mean,
                                                    fill = factor(rrb)))
p.rrb.con = p.rrb.con + geom_bar(stat = "identity",position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.15) + 
  labs(x = "Conditions", y = "Pos 2-5 Mean RT") +
  ggtitle("Mean RT by Condition and RRB") +
  scale_x_discrete("Conditions", labels = 
                     c("1" = "Intact", "2" = "Scrambled-Fixed", 
                       "3" = "Scrambled-Random"))
p.rrb.con = p.rrb.con + coord_cartesian(ylim=c(500,1100))
p.rrb.con

# # # # # # # # # # # # # # # #
# # make a bar graph of position by conditiong with SE
# # get n
# counts = count(group_by(dclean,con.num, position))
# 
# # add n to sumstats
# sumstats.pos[,5] = counts[,3]
# 
# # add SE
# sumstats.pos = mutate(sumstats.pos, SE = SD/sqrt(n))
# 
# # create a bar graph
# limits <- aes(ymax = sumstats.pos$mean + sumstats.pos$SE,
#               ymin = sumstats.pos$mean - sumstats.pos$SE)
# 
# p.pos.con <- ggplot(data = sumstats.pos, aes(x = factor(con.num), y = mean,
#                                              fill = factor(position)))
# p.pos.con = p.pos.con + geom_bar(stat = "identity",
#                                  position = position_dodge(0.9)) +
#   geom_errorbar(limits, position = position_dodge(0.9),
#                 width = 0.15) + 
#   labs(x = "Condition", y = "RT") +
#   ggtitle("RT by Position and Condition") +
#   scale_fill_discrete(name = "Position")  +
#   scale_x_discrete("Conditions", labels = 
#                      c("1" = "Intact", "2" = "Scrambled-Fixed", 
#                        "3" = "Scrambled-Random"))
# p.pos.con = p.pos.con + coord_cartesian(ylim=c(500,1200))
# p.pos.con
# ######################################
# # # # # # # # # # # # # # #
# condition means
# This is so dumb
sumstats.pos.se = dclean %>% group_by(con.num, sub,position) %>% summarise(mean = mean(rt))
# get sample means
con1.pos1 = sumstats.pos.se %>% filter(con.num ==1, position == 1) 
con1.pos2 = sumstats.pos.se %>% filter(con.num ==1, position == 2) 
con1.pos3 = sumstats.pos.se %>% filter(con.num ==1, position == 3) 
con1.pos4 = sumstats.pos.se %>% filter(con.num ==1, position == 4) 
con1.pos5 = sumstats.pos.se %>% filter(con.num ==1, position == 5) 

con2.pos1 = sumstats.pos.se %>% filter(con.num ==2, position == 1) 
con2.pos2 = sumstats.pos.se %>% filter(con.num ==2, position == 2) 
con2.pos3 = sumstats.pos.se %>% filter(con.num ==2, position == 3) 
con2.pos4 = sumstats.pos.se %>% filter(con.num ==2, position == 4) 
con2.pos5 = sumstats.pos.se %>% filter(con.num ==2, position == 5) 

con3.pos1 = sumstats.pos.se %>% filter(con.num ==3, position == 1) 
con3.pos2 = sumstats.pos.se %>% filter(con.num ==3, position == 2) 
con3.pos3 = sumstats.pos.se %>% filter(con.num ==3, position == 3) 
con3.pos4 = sumstats.pos.se %>% filter(con.num ==3, position == 4) 
con3.pos5 = sumstats.pos.se %>% filter(con.num ==3, position == 5) 



# get SD of sammple mean 
sumstats.pos[1,4] = sd(con1.pos1$mean)
sumstats.pos[2,4] = sd(con1.pos2$mean)
sumstats.pos[3,4] = sd(con1.pos3$mean)
sumstats.pos[4,4] = sd(con1.pos4$mean)
sumstats.pos[5,4] = sd(con1.pos5$mean)

sumstats.pos[6,4] = sd(con2.pos1$mean)
sumstats.pos[7,4] = sd(con2.pos2$mean)
sumstats.pos[8,4] = sd(con2.pos3$mean)
sumstats.pos[9,4] = sd(con2.pos4$mean)
sumstats.pos[10,4] = sd(con2.pos5$mean)

sumstats.pos[11,4] = sd(con3.pos1$mean)
sumstats.pos[12,4] = sd(con3.pos2$mean)
sumstats.pos[13,4] = sd(con3.pos3$mean)
sumstats.pos[14,4] = sd(con3.pos4$mean)
sumstats.pos[15,4] = sd(con3.pos5$mean)
# add SE
n = length(unique(dclean$sub))
sumstats.pos = mutate(sumstats.pos, SE = SD/sqrt(n))

# create a bar graph
limits <- aes(ymax = sumstats.pos$mean + sumstats.pos$SE,
              ymin = sumstats.pos$mean - sumstats.pos$SE)

p.pos.con <- ggplot(data = sumstats.pos, aes(x = factor(con.num), y = mean,
                                             fill = factor(position)))
p.pos.con = p.pos.con + geom_bar(stat = "identity",
                                 position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.15) + 
  labs(x = "Condition", y = "RT") +
  ggtitle("RT by Position and Condition") +
  scale_fill_discrete(name = "Position")  +
  scale_x_discrete("Conditions", labels = 
                     c("1" = "Intact", "2" = "Scrambled-Fixed", 
                       "3" = "Scrambled-Random"))
p.pos.con = p.pos.con + coord_cartesian(ylim=c(500,1250))
p.pos.con
# ######################################
# ANOVA
# type 3 SS w/in subs ANOVA
library(afex)
a1 = aov_ez("sub", "rt", d25, within = c("con.num","repetition", "rrb"))
a1
lsmeans(a1,"con.num",contr = "pairwise", adjust = "holm")

d25.con1and2 = d25 %>%  filter(con.num != 3)


# how about just intact and scrmabled
a2 = aov_ez("sub", "rt", d25.con1and2 , within = c("con.num","repetition", "rrb"))
a2
lsmeans(a2,"con.num",contr = "pairwise")


#####################
# look at diffs
#####################

d.diff = d 

d.diff = d.diff %>% mutate(con.num = condition)
d.diff[,15] = ifelse(d.diff[,15] == "intact", 1, ifelse(d.diff[,15] == "scrambled", 2, 3))

d.diff$seq = rep(1:3205, each = 5)

# mean and sd for ID'ing outliers 
d.diff = d.diff %>%  group_by(sub,condition) %>% 
  mutate(mean.sub.con = mean(rt, na.rm = TRUE), SD.sub.con = sd(rt, na.rm = TRUE))

d.diff = d.diff %>% mutate(sd3 = 3 * SD.sub.con) 

# exlude too long
d.diff = d.diff %>% mutate(exclude.high = rt > mean.sub.con + sd3)
# exclude too short
d.diff = d.diff %>% mutate(exclude.low = rt < mean.sub.con - sd3)


# identify too long position 1 RT and NA
test2 = d.diff %>% filter(exclude.high == TRUE | is.na(exclude.high)) %>% filter(position == 1)
# so ~125 position 1s that need to be exluded along with the other 4 positions in the sequence 
exrows.diff = which(d.diff$seq %in% c(test2$seq))
d.diff.clean = d.diff[-c(exrows.diff),]
# check that those sequences are gone
test3 = d.diff.clean %>% filter(exclude.high == TRUE | is.na(exclude.high)) %>% filter(position == 1)
# test 3 is empty so this was successful. 

# now remove the rest of long RTs
# only want complete sequences for this analysis 
# identify too long RT and NA
test4 = d.diff.clean %>% filter(exclude.high == TRUE | is.na(exclude.high))
# 140 trials that will be excluded along with their whole sequence actually makes 700 trials excluded... probably makes this not worthwhile
exrows.diff.2 = which(d.diff.clean$seq %in% c(test4$seq))
d.diff.clean = d.diff.clean[-c(exrows.diff.2),]
# check that those sequences are gone
test4 = d.diff.clean %>% filter(exclude.high == TRUE | is.na(exclude.high)) 



##### only remaining problem would be if there exists a pos 1 valid RT but not valid 2-5 RTs for a given sequence 



# now ready to look at sequence diffs without having mismatch of rows
d.diff.pos1 = d.diff.clean %>% filter(position == 1)
d.diff.25 = d.diff.clean %>% filter(!position == 1)
# for pos 2-5 remove NAs when calculating mean
mean25 = d.diff.25 %>% group_by(sub,block,repetition,condition,seq) %>%  summarise(mean25 = mean(rt)) 
pos1 = d.diff.pos1 %>% group_by(sub,block,repetition,condition,seq) %>% summarise(pos1 = mean(rt))
# these don't have equal numebrs of rows, find the discrepancy
mismatch = mean25[,"seq"] == pos1[,"seq"]
######### mismatch only consists of TRUE, indicating that the seq label for pos 1 matches for all pos 2-5 means


diffmeans = mean25
diffmeans[,7] = pos1[,6]
diffmeans = diffmeans %>% mutate(diffs = pos1 - mean25)
diffmeans = diffmeans %>% mutate(sancheck = diffs < 0)
sum(diffmeans$sancheck, na.rm =TRUE)
####  1936 sequences with faster 2-5 than pos 1 (positive diff), and 950 with the opposite 
# most of those should be from scrambled-random


#also make repetition a factor
diffmeans$repetition = as.factor(diffmeans$repetition)
diffmeans = diffmeans %>% mutate(con.num = condition)
diffmeans[,10] = ifelse(diffmeans[,10] == "intact", 1, ifelse(diffmeans[,10] == "scrambled", 2, 3))


################################## now for the businesss
library(afex)
diffmeans$con.num = as.factor(diffmeans$con.num)
dm1 = aov_ez("sub", "diffs", diffmeans, within = c("con.num","repetition"))
dm1
lsmeans(dm1,"con.num",contr = "pairwise", adjust = "holm")  


# t test between conditions
pairwise.t.test(diffmeans$diffs, diffmeans$con.num, p.adjust="holm", pool.sd = T)


# # # # # # # # # # # # # # graph it # # # # # #
diffmeans.sumstats.rep.con = diffmeans %>% group_by(repetition,con.num) %>% summarise(mean = mean(diffs), SD = sd(diffs))

# trial counts
diffmeans.counts = count(diffmeans, con.num)


# repetition by condition and only pos's 2-5
# get n
counts2 = count(group_by(diffmeans,con.num, repetition))

# add n to sumstats
diffmeans.sumstats.rep.con[,5] = counts2[,3]

# add SE
diffmeans.sumstats.rep.con = mutate(diffmeans.sumstats.rep.con, SE = SD/sqrt(n))

# create a bar graph
limits <- aes(ymax = diffmeans.sumstats.rep.con$mean + diffmeans.sumstats.rep.con$SE,
              ymin = diffmeans.sumstats.rep.con$mean - diffmeans.sumstats.rep.con$SE)

diffsp.rep.con <- ggplot(data = diffmeans.sumstats.rep.con, aes(x = factor(con.num), y = mean,
                                                                fill = factor(repetition)))
diffsp.rep.con = diffsp.rep.con + geom_bar(stat = "identity",position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.15) + 
  labs(x = "Conditions", y = "RT facilitation") +
  ggtitle("RT Facilitation by Condition and Repetition") +
  scale_x_discrete("Conditions", labels = 
                     c("1" = "Intact", "2" = "Scrambled-Fixed", 
                       "3" = "Scrambled-Random"))
diffsp.rep.con = diffsp.rep.con + coord_cartesian(ylim=c(0,350))
diffsp.rep.con

# too noisy to look at between repetition effects 
diffmeans.sumstats.con = diffmeans %>% group_by(con.num) %>% summarise(mean = mean(diffs), SD = sd(diffs))

# repetition by condition and only pos's 2-5
# get n
diff.con.counts = count(group_by(diffmeans,con.num))

# add n to sumstats
diffmeans.sumstats.con[,4] = diff.con.counts[,2]

# add SE
diffmeans.sumstats.con = mutate(diffmeans.sumstats.con, SE = SD/sqrt(n))

# create a bar graph
limits <- aes(ymax = diffmeans.sumstats.con$mean + diffmeans.sumstats.con$SE,
              ymin = diffmeans.sumstats.con$mean - diffmeans.sumstats.con$SE)

diffsp.con <- ggplot(data = diffmeans.sumstats.con, aes(x = factor(con.num), y = mean))
diffsp.con = diffsp.con + geom_bar(stat = "identity",position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.15) + 
  labs(x = "Conditions", y = "RT facilitation") +
  ggtitle("RT Facilitation by Condition") +
  scale_x_discrete("Conditions", labels = 
                     c("1" = "Intact", "2" = "Scrambled-Fixed", 
                       "3" = "Scrambled-Random"))
diffsp.con = diffsp.con + coord_cartesian(ylim=c(0,350))
diffsp.con


##########################################
# random stuff from looking back through data after initial analyses

# why are rrb3 rts so much faster? 

# paired t test for difference between intact and scrambled in one rrb
d25.con1and2.rrb1 = d25.con1and2[which(d25.con1and2$rrb == "rrb1"),]
d25.con1and2.rrb2 = d25.con1and2[which(d25.con1and2$rrb == "rrb2"),]
d25.con1and2.rrb3 = d25.con1and2[which(d25.con1and2$rrb == "rrb3"),]
mean.stats = d25.con1and2.rrb1 %>% group_by(sub, con.num) %>%  summarise(mean = mean(rt), SD = sd(rt))
intact.sub.mean = mean.stats[which(mean.stats$con.num == "1"),3]
scrambled.sub.mean = mean.stats[which(mean.stats$con.num == "2"),3]

t.test(as.matrix(intact.sub.mean), as.matrix(scrambled.sub.mean), paired = T )

maineffect =  as.matrix(scrambled.sub.mean) - as.matrix(intact.sub.mean)


#######
lsmip(a1, con.num ~ repetition)

a6 = aov_ez("sub", "rt", d25.con1and2, within = c("con.num","repetition", "rrb"))
a6
lsmeans(a6,"con.num",contr = "pairwise")
lsmip(a6, con.num ~ rrb)

a3 = aov_ez("sub", "rt", d25.con1and2.rrb1, within = c("con.num","repetition"))
a3
lsmeans(a3,"con.num",contr = "pairwise")
lsmip(a3, con.num ~ rrb)

a4 = aov_ez("sub", "rt", d25.con1and2.rrb2, within = c("con.num","repetition"))
a4
lsmeans(a4,"con.num",contr = "pairwise")
lsmip(a4, con.num ~ rrb)

a5 = aov_ez("sub", "rt", d25.con1and2.rrb3, within = c("con.num","repetition"))
a5
lsmeans(a5,"con.num",contr = "pairwise")
lsmip(a5, con.num ~ rrb)

# 16.71 * 5
# 48.84/83.55
# from a5 pairwise contrast mean difference and SE gives effect size of .58
# why is this so different from ges? 

r1 <- lsmeans(a3, ~con.num )
r1

# compute cohen's d using halle's function
# get summary stats for conditions and subjects 
library(halle)
d.stats = d25.con1and2 %>% group_by(sub,con.num) %>%  summarize(mean = mean(rt))
intact.sub.mean = d.stats[which(d.stats$con.num == "1"),3]
scrambled.sub.mean = d.stats[which(d.stats$con.num == "2"),3]
cohensd.frame = as.data.frame(c(intact.sub.mean, scrambled.sub.mean)) 
colnames(cohensd.frame) =  c("var1","var2")
compute_cohens_d(cohensd.frame)

# what about effect size calculated for each rrb
# get summary stats for conditions and subjects 
d.stats = d25.con1and2.rrb2 %>% group_by(sub,con.num) %>%  summarize(mean = mean(rt))
intact.sub.mean = d.stats[which(d.stats$con.num == "1"),3]
scrambled.sub.mean = d.stats[which(d.stats$con.num == "2"),3]
cohensd.frame = as.data.frame(c(intact.sub.mean, scrambled.sub.mean)) 
colnames(cohensd.frame) =  c("var1","var2")
compute_cohens_d(cohensd.frame)

#same for diffs
diffmeans.con1and2 =  diffmeans %>%  filter(con.num != 3)
d.stats = diffmeans.con1and2 %>% group_by(sub,con.num) %>%  summarize(mean = mean(diffs))
intact.sub.mean = d.stats[which(d.stats$con.num == "1"),3]
scrambled.sub.mean = d.stats[which(d.stats$con.num == "2"),3]
cohensd.frame = as.data.frame(c(intact.sub.mean, scrambled.sub.mean)) 
colnames(cohensd.frame) =  c("var1","var2")
compute_cohens_d(cohensd.frame)

#####
# how about effect size for conditions 1 and 3
d25.con1and3 = d25 %>% filter(con.num != 2)
d.stats = d25.con1and3 %>% group_by(sub,con.num) %>%  summarize(mean = mean(rt))
intact.sub.mean = d.stats[which(d.stats$con.num == "1"),3]
scrambled.sub.mean = d.stats[which(d.stats$con.num == "3"),3]
cohensd.frame = as.data.frame(c(intact.sub.mean, scrambled.sub.mean)) 
colnames(cohensd.frame) =  c("var1","var2")
compute_cohens_d(cohensd.frame)

################################################
# look at test data 
sms2_test_acc <- read_csv("~/walter/dml/sms2/sms2_test_acc",col_names = FALSE)
sms2_test_acc = sms2_test_acc %>% mutate(z.test2 = scale(X2))
sms2_test_acc[,4] = x
sms2_test_acc = sms2_test_acc %>% mutate(sub = 8:33)
cor.test(sms2_test_acc$z.test2,sms2_test_acc$V4)
cor.test(sms2_test_acc$X1,sms2_test_acc$V4)
# learning test is not diagnostic whatsoever 

#######################################################################
# look at responses
dsansnas$response.num = ifelse(dsansnas$response == "yes", 1,0)
dsansnas.con1 = dsansnas %>% filter(con.num == 1)
dsansnas.con2 = dsansnas %>% filter(con.num == 2)
dsansnas.con3 = dsansnas %>% filter(con.num == 3)

sum.response.stats1 = dsansnas.con1 %>% group_by(seq_num) %>% summarise(mean = mean(response.num), sd = sd(response.num))
sum.response.stats2 = dsansnas.con2 %>% group_by(seq_num) %>% summarise(mean = mean(response.num),sd = sd(response.num))
sum.response.stats3 = dsansnas.con3 %>% group_by(seq_num) %>% summarise(mean = mean(response.num),sd = sd(response.num))

responsedf = as.data.frame(cbind(sum.response.stats1$seq_num,sum.response.stats1$mean,sum.response.stats2$mean,sum.response.stats3$mean))
colnames(responsedf) =  c("seq_num","intact","scrambled","scrambled-random")
responsedf = responsedf %>% gather("con","resp", 2:4)

# these didn't work. Not sure why.
# rdf1 = aov_ez("seq_num","resp",responsedf, between = "con")
# rdf1 = aov_car(resp ~ con + (seq_num/con),responsedf)

rdf1 = aov(resp ~ con + (seq_num),responsedf)
summary(rdf1)
# perfect. no difference in yes/nos between conditions

###################################################
# ID participants based on how frequently they give the same response 
# x = 
x= 0 
for (i in 2:dim(dsansnas)[1]){
  if (dsansnas[i,1] == dsansnas[i-1,1]){
    if (dsansnas[i,22] == dsansnas[i-1,22]){
      dsansnas[i,23] = 1 + x; x = x+1
    } else {
      x = 0; dsansnas[i,23] = 0
    }
  } else {
    x = 0; dsansnas[i,23] = 0
  }
}


# here goes
nonono = dsansnas %>% group_by(sub) %>% summarise(mean = mean(V23), sd = sd(V23))
nonono = nonono %>% mutate(zmean = scale(mean))
# just get rid of sub 10


# look at subs who didn't have main effect
maineffect = d25 %>% filter(con.num != 3) %>% group_by(sub,con.num) %>% summarise(mean = mean(rt))
maineffect = maineffect %>% spread(con.num,mean)
colnames(maineffect) = c("sub","intact","scrambled")
maineffect = maineffect %>% mutate(intact.minus.scrambled = intact - scrambled)
subs.no.main.effect = c(31, 10, 33, 28, 9, 25)


########################################################################
# mixed effect model? 
# ICC to see if seq accounts for large amount of variance 
d25$seq_num = as.factor(d25$seq_num)
d25$item = as.factor(d25$item)
d25$con.num = as.factor(d25$con.num)
d25$repetition = as.factor(d25$repetition)
d25$sub = as.factor(d25$sub)

# mm2 = mixed(rt ~ con.num * rrb + ( 1 | sub) + (seq_num || item),d25, expand_re = TRUE, progress = TRUE)

mm1 = lmer(rt ~ con.num * rrb + ( con.num * rrb || sub) + (1  | seq_num/item),d25)
#failed to converge
mm2 = lmer(rt ~ con.num * rrb + ( con.num + rrb || sub) + (1  | seq_num/item),d25)
# failed to converge
mm3 = lmer(rt ~ con.num * rrb + ( con.num * rrb || sub) + (1  | item),d25)
#failed to converge
mm4 = lmer(rt ~ con.num * rrb + ( con.num + rrb || sub) + (1  | item),d25)
# model is nearly unidentifiable 
# mx4 = mm4 = mixed(rt ~ con.num * rrb + ( con.num + rrb || sub) + (1  | item),d25)
# RAM'd out

mm5 = lmer(rt ~ con.num * repetition * rrb + ( con.num + repetition + rrb || sub) + (1  | item),d25)
# failed to converge
mm6 = lmer(rt ~ con.num * repetition * rrb + ( con.num + repetition + rrb || sub)  ,d25)
# failed to converge
mm7 = lmer(rt ~ con.num * repetition * rrb + ( con.num + rrb || sub) + (1  | item),d25)
# not helpful to include repetition
mm8 = lmer(rt ~ con.num + ( con.num + rrb + repetition || sub) + (1  | item),d25)
# ok 
mm9 = lmer(rt ~ con.num + ( con.num + rrb + repetition || sub) + (item || seq_num) ,d25)
# failed

mm10 = lmer(rt ~ 1 +  (item  || seq_num),d25)

# hold that thought :  con.num + ( con.num + rrb || sub) +



###############
# misc for talk
d.pos1 = dsansnas %>% filter(position == 1)
library(afex)
d.pos1$con.num = as.factor(d.pos1$con.num)
d.pos1$rrb = as.factor(d.pos1$rrb)
d.pos1$repetition = as.factor(d.pos1$repetition)

apos1 = aov_ez("sub", "rt", d.pos1, within = c("con.num","repetition", "rrb"))
apos1
lsmeans(apos1,"con.num",contr = "pairwise", adjust = "holm")  
