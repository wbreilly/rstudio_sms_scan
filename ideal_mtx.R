#  idealized rsa mtx


# values for staircase matlab style matrix indexing

fill_vals = c(1, .5, 0, 0, 0, 1, .5,0,0,1,.5,0,1,.5,1)
mtx.ideal = as.data.frame(fill_vals)
mtx.ideal$x = c(1,1,1,1,1,2,2,2,2,3,3,3,4,4,5)
mtx.ideal$y = c(5,4,3,2,1,4,3,2,1,3,2,1,2,1,1)

mtx.PM.intact = ggplot(mtx.ideal, aes(x = x, y = y, fill = fill_vals)) +
  geom_tile() + scale_fill_viridis(na.value = "transparent",option = "inferno",limits = (c(0,1)))  + theme(axis.text.y = element_blank(), axis.text.x = element_blank())
mtx.PM.intact = mtx.PM.intact + theme(axis.line=element_blank(),
                                      axis.text.x=element_blank(),
                                      axis.text.y=element_blank(),
                                      axis.ticks=element_blank(),
                                      axis.title.x=element_blank(),
                                      axis.title.y=element_blank(),
                                      legend.position="none",
                                      legend.title = element_blank(),
                                      panel.background=element_blank(),
                                      panel.border=element_blank(),
                                      panel.grid.major=element_blank(),
                                      panel.grid.minor=element_blank(),
                                      plot.background=element_blank(),
                                      plot.title = element_blank(),
                                      text = element_text(size = 36, face = "bold"),
                                      legend.position = c(0.75, 0.8))
mtx.PM.intact
ggsave("design_mtx.eps", plot = last_plot(), dpi = 600 )
