# recognition data script
# started 1.8.18

library(tidyverse)
library(afex)

# import recognition data
d = read.csv("~/drive/grad_school/DML_WBR/Sequences_Exp3/sms_scan_drive/recognition_for_r/all_recog_1_8_18.csv",na.strings="NaN", stringsAsFactors=FALSE)

# import condition keys
df.key = read.csv("~/drive/grad_school/DML_WBR/Sequences_Exp3/sms_scan_drive/recognition_for_r/all_recog_keys_1_8_18.csv",na.strings="NaN", stringsAsFactors=FALSE)

#get rid of header info
d = d[!(d$item == ""), ]
d = d[!(d$sub == "sub"), ]
# do these things
d$sub = as.factor(as.numeric(as.character(d$sub)))
df.key$sub = as.factor(as.numeric(as.character(df.key$sub)))
# rt matches the .dat
d$rt = as.numeric(as.character(d$rt))

# should use join... but would need to go back and add item to the df key in matlab... bind works just not as elegant 
# anything that changes order would disrupt this scheme..  but I don't think missing trials would
# d = inner_join(d,df.key, by = "sub","item") should do it
d$condition = df.key$condition

d$response = as.factor(d$response)
d$condition = as.factor(d$condition)


# need to categorize lures dangit...got ahead of myself

# sum stats
 sum.stats = d %>% group_by(sub,condition) %>% count(response)

 # get ride of lure for now
 sum.stats = sum.stats %>% filter(condition != "lure")
 
 # anova on count data? maybe should be doing x^2...?
 m1 = aov_ez("sub", "n", sum.stats, within = c("condition","response"))
 m1
 # that's not right... plot first
 plot.dat = d %>% group_by(condition) %>% count(response)
 plot.dat
 # cool scrambled has more remember and fewer misses (new) than intact!
 # wanna scatter data so not using plot.dat :(
 
a.plot = ggplot(data = sum.stats, aes(x = factor(response), y = n)) + geom_point(aes(colour = factor(sub)))  + facet_wrap(~condition)
a.plot
 
