theme_pub <- function(base_size = 12, base_family = "Helvetica")
{
    half_line <- base_size/4
    theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(
        axis.title=element_text(size=6),
        axis.text=element_text(size=6, colour="black"),
        axis.line = element_line(colour = "black"),
        legend.position="none",
        legend.text=element_text(size=6),
        legend.title=element_text(size=7),
        legend.key = element_blank(),
        legend.key.size=unit(0.5, units="lines"),
        legend.background = element_blank(),
        line=element_line(colour = "black", size = 0.12, linetype = 1, lineend = "butt"),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black", size = 0.35, linetype="solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing.y=unit(0.55, units="lines"),
        panel.spacing.x=unit(0.2, units="lines"),
        plot.margin=unit(c(2,2,0.1,0.1), units="mm"), # modified 7/11/17 to correct errors
        plot.title=element_text(size=10),
        rect=element_rect(fill="white", colour="black", size=0.35, linetype=1),
        strip.background = element_blank(),
        strip.text.y=element_text(size=6, margin = margin(t = half_line, b = half_line)),
        strip.text.x=element_text(size=6, margin = margin(t = half_line, b = half_line)),
        complete = TRUE      
    )
}
