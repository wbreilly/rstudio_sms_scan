# recognition data script
# started 1.8.18

library(tidyverse)
library(afex)

# import recognition data
d = read.csv("~/drive/grad_school/DML_WBR/Sequences_Exp3/sms_scan_drive/recognition_for_r/all_recog_1_8_18.csv",na.strings="NaN", stringsAsFactors=FALSE)
# import condition keys
df.key = read.csv("~/drive/grad_school/DML_WBR/Sequences_Exp3/sms_scan_drive/recognition_for_r/all_recog_keys_1_11_18.csv",na.strings="NaN", stringsAsFactors=FALSE)


#get rid of header info
d = d[!(d$item == ""), ]
d = d[!(d$sub == "sub"), ]
# do these things
d$sub = as.factor(as.numeric(as.character(d$sub)))
df.key$sub = as.factor(as.numeric(as.character(df.key$sub)))
# rt matches the .dat
d$rt = as.numeric(as.character(d$rt))

# bind works just not as elegant 
# d = inner_join(d,df.key, by = "sub","item") should do it
d = inner_join(d,df.key)


d$response = as.factor(d$response)
d$condition = as.factor(d$condition)


# need to categorize lures dangit...got ahead of myself

# sum stats
 sum.stats = d %>% group_by(sub) %>% count(condition)
 # appropriate trial numbers for all conditions and participants are present

#  # plot first
#  plot.dat = d %>% group_by(condition) %>% count(response)
#  plot.dat
#  # cool scrambled has more remember and fewer misses (new) than intact!
#  # wanna scatter data so not using plot.dat :(
#  
# a.plot = ggplot(data = sum.stats, aes(x = factor(response), y = n)) + geom_point(aes(colour = factor(sub)))  + facet_wrap(~condition)
# a.plot
 
 d$hits_intact = ifelse(d$condition == "intact" & d$response == "rem",1,0)
 d$hits_scrambled = ifelse(d$condition == "scrambled" & d$response == "rem", 2,0)
 d$hits_random = ifelse(d$condition == "random" & d$response == "rem", 3,0)
 
 d$fa_intact = ifelse(d$condition == "lure_intact" & d$response == "rem", 4,0)
 d$fa_scrambled = ifelse(d$condition == "lure_scrambled" & d$response == "rem", 5,0)
 d$fa_random = ifelse(d$condition == "lure_random" & d$response == "rem", 6,0)
 
  # Unite
 d = d %>% mutate(hits_fas = hits_intact + hits_scrambled + hits_random + fa_intact + fa_scrambled + fa_random)
 d$hits_fas = if_else(d$hits_fas == "1","hit_intact", if_else(d$hits_fas == "2", "hit_scrambled", if_else(d$hits_fas == "3", "hit_random",
                      if_else(d$hits_fas == "4", "fa_intact", if_else(d$hits_fas == "5", "fa_scrambled", if_else(d$hits_fas == "6", "fa_random","0"))))))
 
 # sum stats
 sum.stats = d %>% group_by(sub) %>% count(hits_fas)
 
 #convert number of hits and false alarms to appropriate to proportions
 # 30 for intact/scram hits, 15 for random hits. For FAs, 18, 18 and 9 
 sum.stats = data.frame(sum.stats)
 sum.stats1 = sum.stats %>% filter(hits_fas == "hit_intact") %>% mutate(props = n/30)
 sum.stats2 = sum.stats %>% filter(hits_fas == "hit_scrambled") %>% mutate(props = n/30)
 sum.stats3 = sum.stats %>% filter(hits_fas == "hit_random") %>% mutate(props = n/15)
 sum.stats4 = sum.stats %>% filter(hits_fas == "fa_intact") %>% mutate(props = n/18)
 sum.stats5 = sum.stats %>% filter(hits_fas == "fa_scrambled") %>% mutate(props = n/18)
 sum.stats6 = sum.stats %>% filter(hits_fas == "fa_random") %>% mutate(props = n/9)
 
 # trying to put all those dfs together now. Join wants to make new columns, but it's close
 # ddprime = full_join(sum.stats1, sum.stats2, by = c("sub","hits_fas"))
 ddprime = bind_rows(sum.stats1,sum.stats2,sum.stats3,sum.stats4,sum.stats5,sum.stats6)
 # boom
 
# z score 
 ddprime = ddprime %>% mutate(zhitsfas = qnorm(props))
 
 # dprime = function(hit, fa) {
 #   qnorm(hit) - qnorm(fa)
 # }
 
 
 

 