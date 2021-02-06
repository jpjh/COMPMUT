## Function for calculating fitness
# Requires tidyverse

calculate_fitness <- function(df, id_cols, start_dilfac=0.01, end_dilfac=1){
  df$gm_ml <- with(df, (1000/spread) * (10^dilution) * count_white)
  df$sm_ml <- with(df, (1000/spread) * (10^dilution) * count_blue)
  
  # Note: .x indicates starts, .y indicates ends.
  df_w <- left_join(
    df %>% filter(timepoint=="start"),
    df %>% filter(timepoint=="end"),
    by=all_of(id_cols),
    suffix=c(".start",".end")
  ) %>% select(all_of(id_cols), gm_ml.start, sm_ml.start, gm_ml.end, sm_ml.end)
  
  # calculate m per ml
  # m is ln(end_counts/start_counts)
  # For KB experiments (default), start_dilfac is 0.01 (1% transferred) and end_dilfac is 1
  # For soil, counts are calculated per gram of soil:
  #    Overnight mixes are diluted 1:10 before adding 0.1 ml per 10 g.
  #    So it's the equivalent of 0.01 ml starting mix per 10 g.
  #      so start_dilfac is 0.01/10 = 0.001.
  #    For the end counts, it's a volume of 11 ml from 10 g.
  #      so end_dilfac is 11/10 = 1.1.
  
  df_w %>% 
    mutate(
      start_gm  = start_dilfac * gm_ml.start,
      end_gm    = end_dilfac * gm_ml.end,
      start_sm  = start_dilfac * sm_ml.start,
      end_sm    = end_dilfac * sm_ml.end,
      
      m_gm      = log(end_gm/start_gm),
      m_sm      = log(end_sm/start_sm),
      W_gm      = m_gm/m_sm
    ) %>%
    select(all_of(id_cols), W_gm)
}
