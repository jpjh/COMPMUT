theme_pub <- function(base_size = 6, base_family = "Helvetica")
{
    half_line <- 0.5*base_size
    theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(
        axis.text=element_text(colour="black"),
        axis.line = element_line(colour = "black"),
        legend.position="none",
        legend.title=element_text(size=rel(1.2)),
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
        plot.margin=unit(c(2,2,0.1,0.1), units="mm"),
        plot.title=element_text(size=rel(1.6)),
        rect=element_rect(fill="white", colour="black", size=0.35, linetype=1),
        strip.background = element_blank(),
        strip.text.y=element_text(size=rel(1.2), 
                                  margin = margin(t = half_line, b = half_line)),
        strip.text.x=element_text(size=rel(1.2), 
                                  margin = margin(t = half_line, b = half_line)),
        complete = TRUE      
    )
}
