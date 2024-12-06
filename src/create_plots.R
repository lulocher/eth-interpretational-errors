library(arrow)
library(ggplot2)
library(dplyr)
library(wesanderson)

###
# Read and prepare data
###

path_data <- './data/data_interpretational_errors.parquet'
df <- read_parquet(path_data)

df_mutated <- mutate(
  df,
  total='Total',
  application_type = factor(application_type, levels=c('MR', 'RCT', 'Other')),
  discipline = factor(discipline, levels=c('Epidemiology', 'Medicine', 'Economics', 'Political Science')),
  framework_manual = case_when(
    framework == 'Pragmatic' ~ causal_framework,
    estimand_of_interest == 'LATE' ~ 'Ideal - Interest in LATE',
    T ~ 'Ideal - Targeted ATE'      
  ),
  estimand_specification = ifelse(estimand_explicitly_specified, 'Explicit specification', 'No explicit specification'),
  slippage_at_least_once = (!is.na(identity_slippage_abstract) & identity_slippage_abstract) |  (!is.na(identity_slippage_results) & identity_slippage_results) | (!is.na(identity_slippage_discussion) & identity_slippage_discussion),
  slippage_w_claim_type_abstract = ifelse(
    !is.na(identity_slippage_abstract),
    paste(ifelse(identity_slippage_abstract, 'Yes', 'No'), ' (', class_claims_abstract, ')', sep = ''),
    NA
  ),
  slippage_w_claim_type_results = ifelse(
    !is.na(identity_slippage_results),
    paste(ifelse(identity_slippage_results, 'Yes', 'No'), ' (', class_claims_results, ')', sep = ''),
    NA
  ),
  slippage_w_claim_type_discussion = ifelse(
    !is.na(identity_slippage_discussion),
    paste(ifelse(identity_slippage_discussion, 'Yes', 'No'), ' (', class_claims_discussion, ')', sep = ''),
    NA
  )
)

df_filtered <- filter(df_mutated, targeted_estimand == 'LATE')

###
# Specify variables / functions
###

colors <- wes_palette('GrandBudapest2')

create_similar_colors <- function(original_color){
  rgb_color <- hex2RGB(original_color)
  hcl_color <- as(rgb_color, "polarLUV")
  color1 <- polarLUV(hcl_color@coords[1] - 8, hcl_color@coords[2], hcl_color@coords[3])
  color2 <- polarLUV(hcl_color@coords[1] - 16, hcl_color@coords[2], hcl_color@coords[3])
  color1_hex <- hex(color1)
  color2_hex <- hex(color2)
  return(
    c(color1_hex, color2_hex)
  )
}

construct_data_framework <- function(data, grouping_var){
  data['group'] <- data[grouping_var]
  
  counts <- summarise(
    group_by_at(data, c('causal_framework', 'group')),
    n=n()
  )
  
  counts_mutated <- mutate(
    mutate(group_by(counts, group), sum = sum(n)),
    perc=(n / sum) * 100
  )
  
  counts_mutated['framework_manual'] <- counts_mutated['causal_framework']
  
  return(counts_mutated)
}

create_bar_plot <- function(data, ylab=element_blank()){ # Drop axis description
  scheme <- c('#BE8AAC', colors[1], colors[2])
  ggplot(data=data, aes(x=group, y=n, fill=framework_manual, label=sprintf("%.0f%% (%d)", perc, n))) +
    geom_bar(stat='identity') +
    scale_fill_manual(
      values=scheme
    ) +
    labs(x=element_blank(), y=ylab, fill='Causal framework') +
    coord_flip() +
    theme_minimal() +
    theme(
      legend.text = element_text(hjust = 0),
      legend.position = 'top',
      legend.direction = 'horizontal'
    )
}

generate_framework_plot <- function(alt_setting=F){
  df_all <- construct_plot_data(df, 'total', alt_setting)
  df_discipline <- construct_plot_data(df, 'discipline', alt_setting)
  df_application_types <- construct_plot_data(df, 'application_type', alt_setting)
  
  plot_all <- create_bar_plot(df_all)
  plot_disciplines <- create_bar_plot(df_discipline)
  plot_application_types <- create_bar_plot(df_application_types, 'Number of articles')
  
  plot_frameworks <- (guide_area() / plot_all / plot_disciplines / plot_application_types) + plot_layout(guides = 'collect', heights = c(1, 1, 4, 3))
  
  return(plot_frameworks)
}

construct_data_slippage <- function(data, grouping_var){
  data_grouped <- group_by_at(data, grouping_var)
  
  summarised <- summarise(
    data_grouped,
    n_total = n(),
    n_slip_abstract = sum(identity_slippage_abstract, na.rm=T),
    n_no_slip_abstract = sum(!identity_slippage_abstract, na.rm=T),
    n_slip_results = sum(identity_slippage_results, na.rm=T),
    n_no_slip_results = sum(!identity_slippage_results, na.rm=T),
    n_slip_discussion = sum(identity_slippage_discussion, na.rm=T),
    n_no_slip_discussion = sum(!identity_slippage_discussion, na.rm=T)
  )
  
  summarised_extended <- mutate(
    summarised, 
    n_no_iv_abstract = n_total - n_slip_abstract - n_no_slip_abstract,
    n_no_iv_results = n_total - n_slip_results - n_no_slip_results,
    n_no_iv_discussion = n_total - n_slip_discussion - n_no_slip_discussion
  )
  
  counts <- mutate(
    rbind(
      mutate(select(summarised_extended, c(n=n_slip_abstract, group=grouping_var)), section='Abstract', slippage = 'Yes'),
      mutate(select(summarised_extended, c(n=n_no_slip_abstract, group=grouping_var)), section='Abstract', slippage = 'No'),
      mutate(select(summarised_extended, c(n=n_slip_results, group=grouping_var)), section='Results', slippage = 'Yes'),
      mutate(select(summarised_extended, c(n=n_no_slip_results, group=grouping_var)), section='Results', slippage = 'No'),
      mutate(select(summarised_extended, c(n=n_slip_discussion, group=grouping_var)), section='Discussion', slippage = 'Yes'),
      mutate(select(summarised_extended, c(n=n_no_slip_discussion, group=grouping_var)), section='Discussion', slippage = 'No')
    ), 
    section = factor(section, c('Results', 'Discussion', 'Abstract')),
    slippage = factor(slippage, c('Yes', 'No'))
  )
  
  counts_mutated <- mutate(
    group_by(counts, section, group),
    total = sum(n),
    perc = (n / total)
  )
  
  return (counts_mutated)
}

construct_data_claim_class <- function(data, grouping_var){
  data_mutated <- mutate(group_by_at(data, grouping_var), n_total_group = n())
  
  summarised_abstract <- summarise(
    group_by_at(data_mutated, c(grouping_var, 'slippage_w_claim_type_abstract')),
    n=n(),
  )
  
  summarised_results <- summarise(
    group_by_at(data, c(grouping_var, 'slippage_w_claim_type_results')),
    n=n()
  )
  summarised_discussion <- summarise(
    group_by_at(data, c(grouping_var, 'slippage_w_claim_type_discussion')),
    n=n()
  )
  
  counts = filter(
    rbind(
      mutate(select(summarised_abstract, c(n, group=grouping_var, claim_type = slippage_w_claim_type_abstract)), section = 'Abstract'),
      mutate(select(summarised_results, c(n, group=grouping_var, claim_type = slippage_w_claim_type_results)), section = 'Results'),
      mutate(select(summarised_discussion, c(n, group=grouping_var, claim_type = slippage_w_claim_type_discussion)), section = 'Discussion')
    ), !is.na(claim_type)
  )
  
  all_combinations = expand.grid(
    group = unique(data[[grouping_var]]),
    claim_type=c('No (a)', 'No (b)', 'No (c)', 'Yes (d)', 'Yes (e)', 'Yes (f)'),
    section=c('Abstract', 'Results', 'Discussion')
  )
  
  counts_filled <- mutate(
    left_join(
      all_combinations,
      counts,
      by=c('group', 'claim_type', 'section')
    ),
    n = case_when(
      is.na(n) ~ 0,
      TRUE ~ n
    ),
    section = factor(section, c('Results', 'Discussion', 'Abstract')),
    claim_type = factor(claim_type, c('No (a)', 'No (b)', 'No (c)', 'Yes (d)', 'Yes (e)', 'Yes (f)')) # Specify ordering in plot
  )
  
  counts_mutated <- mutate(
    group_by(counts_filled, section, group),
    total = sum(n),
    perc = (n / total)
  )
  
  return (counts_mutated)
}

create_slippage_plot <- function(data, fill_var, guide_title, y_lab=element_blank()){
  if (fill_var=='slippage'){
    scheme <- c(colors[1:2])
    breaks <- c('Yes', 'No')
  } else {
    similar_slip <- create_similar_colors(colors[1])
    similar_no_slip <- create_similar_colors(colors[2])
    scheme <- c(colors[2], similar_no_slip, colors[1], similar_slip)
    breaks <- c(c('No (a)', 'No (b)', 'No (c)', 'Yes (d)', 'Yes (e)', 'Yes (f)'))
  }
  data['fill'] = data[fill_var]
  
  plot_out <- ggplot(data=data, aes(x=section, y=perc, fill=fill, group=fill, label=sprintf("%.2f (%d)", perc, n))) + # TODO: Drop labels
    geom_bar(stat = 'identity', width=0.6, alpha=1) +
    geom_area(aes(x = c('Results' = 1.3, 'Discussion' = 1.7)[section]), alpha=0.6) +
    geom_area(aes(x = c('Discussion' = 2.3, 'Abstract' = 2.7)[c(
      section[(length(section)*(2/3) + 1):length(section)], section[1:(length(section)*(2/3))]
    )]), alpha=0.6) +
    scale_fill_manual(values=scheme, breaks = breaks, guide = guide_legend(nrow = 1)) +
    scale_y_continuous(labels = label_percent()) +
    labs(x=element_blank(), y=y_lab, fill=guide_title) +
    theme_minimal() +
    theme(
      legend.text = element_text(hjust = 0),   # Left-align the legend text
      legend.position = 'top',
      legend.direction = 'horizontal',
      axis.text.x = element_text(size=7)
    )
  
  if (length(unique(data$group)) > 1){plot_out <- plot_out + facet_wrap(~group, nrow=1)}
  return (plot_out)
}

###
# Create figures
###

# Framework figures
plot_frameworks_main <- generate_framework_plot()

# Slippage figures
guide_title_slippage <- guide_title_class <- 'Identity slippage'

data_total <- construct_data_slippage(df_filtered, 'total')
data_explicitly_specified <- construct_data_slippage(df_filtered, 'estimand_specification')
data_framework <- construct_data_slippage(df_filtered, 'framework')
data_disciplines <- construct_data_slippage(df_filtered, 'discipline')
data_application_types <- construct_data_slippage(df_filtered, 'application_type')

plot_total_slippage <- create_slippage_plot(data_total, 'slippage',guide_title_slippage, y_lab='% of studies')
plot_specification_slippage <- create_slippage_plot(data_explicitly_specified, 'slippage', guide_title_slippage) 
plot_framework_slippage <- create_slippage_plot(data_framework, 'slippage', guide_title_slippage)

plot_slippage_combined <- (plot_total_slippage | (plot_specification_slippage / plot_framework_slippage))  + plot_layout(widths = c(1, 2))
plot_main_slippage <- guide_area() / plot_slippage_combined  + plot_layout(guides = 'collect', heights = c(1, 9)) + plot_annotation(tag_levels = 'a', tag_suffix = ')')

plot_disciplines_slippage <- create_slippage_plot(data_disciplines, 'slippage', guide_title_slippage, y_lab='% of studies')
plot_application_types_slippage <- create_slippage_plot(data_application_types, 'slippage', guide_title_slippage, y_lab='% of studies')
bottom_row <- (plot_application_types_slippage | plot_spacer()) + plot_layout(widths = c(3.25, 1))
plot_disciplines_applications_slippage <- (guide_area() / plot_disciplines_slippage / bottom_row) + plot_layout(guides='collect', heights = c(1, 9, 9)) + plot_annotation(tag_levels = 'a', tag_suffix = ')')

# Granular slippage figures
data_total_claim_class <- construct_data_claim_class(df_filtered, 'total')
data_explicitly_specified_claim_class <- construct_data_claim_class(df_filtered, 'estimand_specification')
data_framework_claim_class <- construct_data_claim_class(df_filtered, 'framework')
data_disciplines_claim_class <- construct_data_claim_class(df_filtered, 'discipline')
data_application_types_claim_class <- construct_data_claim_class(df_filtered, 'application_type')

plot_total_claim_class <- create_slippage_plot(data_total_claim_class, 'claim_type', guide_title_class, y_lab='% of studies')
plot_specification_claim_class <- create_slippage_plot(data_explicitly_specified_claim_class, 'claim_type', guide_title_class)
plot_framework_claim_class <- create_slippage_plot(data_framework_claim_class, 'claim_type', guide_title_class)

plot_claim_class_combined <- (plot_total_claim_class | (plot_specification_claim_class / plot_framework_claim_class)) + plot_layout(widths = c(1, 2))
plot_main_claim_class <- guide_area() / plot_claim_class_combined  + plot_layout(guides = 'collect', heights = c(1, 9)) + plot_annotation(tag_levels = 'a', tag_suffix = ')')

plot_disciplines_claim_class <- create_slippage_plot(data_disciplines_claim_class, 'claim_type', guide_title_class, y_lab='% of studies')
plot_application_types_claim_class <- create_slippage_plot(data_application_types_claim_class, 'claim_type', guide_title_class,y_lab='% of studies')
bottom_row <- (plot_application_types_claim_class | plot_spacer()) + plot_layout(widths = c(3.25, 1))
plot_disciplines_applications_claim_class <- (guide_area() / plot_disciplines_claim_class / bottom_row) + plot_layout(guides='collect', heights = c(1, 9, 9)) + plot_annotation(tag_levels = 'a', tag_suffix = ')')

###
# Save figures
###

default_width = 2200
default_height = 1200

ggsave('./figures/plot_framework_main.png', plot_frameworks_main, width=default_width, height=default_height, units='px')
ggsave('./figures/plot_slippage_main.png', plot_main_slippage, width=default_width, height=default_height, units='px')
ggsave('./figures/plot_slippage_discipline_and_application_type.png', plot_disciplines_applications_slippage, width=default_width, height=default_height, units='px')
ggsave('./figures/plot_claim_class_main.png', plot_main_claim_class, width=default_width, height=default_height, units='px')
ggsave('./figures/plot_claim_class_discipline_and_application_type.png', plot_disciplines_applications_claim_class, width=default_width, height=default_height, units='px')