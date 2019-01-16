vox = read.table("/Users/wbr/walter/fmri/rsfc_glasser_ROIs/12_11_18_halle_glasser_roi_dat.txt", header = TRUE, stringsAsFactors = TRUE)

library(tidyverse)

sum_vox = vox %>% group_by(roi) %>% summarize(meanvox = mean(nvox))

write.csv(sum_vox,"/Users/wbr/walter/fmri/rsfc_glasser_ROIs/vox_means.txt")

connames  = read.table("conn_names.csv")
glass_lut = read.table("glasser_lut.txt",header = TRUE, stringsAsFactors = TRUE)

glass_lut$keep = glass_lut$roi %in% connames$V1

write.csv(glass_lut,"rois_to_keep.csv")
