# Figure 2

library(tidyverse)
library(multcomp)
library(survival)
library(coxme)
library(broom)

library(ggbeeswarm)
library(scico)
library(scales)
library(ggpubr)
library(cowplot)


source("clean_data.R")
source("cumulative_branches.R")

# load data
total <- clean_data()
total_filt <- total %>% filter(Location != "Unclear")
#total_df <- as.data.frame(total_filt) # survminer can't handle tibbles


## Statistics ----

# Branching events by type WT vs KO
count_exp <- total %>% 
  group_by(Branchtype, Culture, Genotype, CellSum) %>% 
  summarise(events = n(),
            Neurons = sum(unique(Neurons))) %>% 
  mutate(norm_events = events/Neurons ) %>% 
  ungroup()



# Testing model assumptions
m0_comb <- aov(norm_events ~ Genotype, data = count_exp)
autoplot(m0_comb, which = 1:6, ncol = 2, label.size = 3, color = Genotype)


# Welchs's t-test 
m0_combined <- compare_means(norm_events ~ Genotype, group.by = "Branchtype", 
                             data = count_exp, method = "t.test", p.adjust.method = "holm")



## Survival analysis
m_combined <- coxph(Surv(Lifetime, CompleteData) ~ 
                      Branchtype * Location + Genotype,
                    weights = Prob,
                    data = total)

# Diagnostics
diag_surv_combined<- cox.zph(m_combined)
survminer::ggcoxzph(diag_surv_combined) # Schoenfeld test
survminer::ggcoxdiagnostics(m_combined, type = "deviance", linear.predictions = FALSE) # 

# Results table
tidy(m_combined, exponentiate = T)



m_combined3 <- coxph(Surv(Lifetime, CompleteData) ~ 
                      Branchtype * Location + Genotype + cluster(Culture),
                    weights = Prob,
                    data = total)

tidy(m_combined3, exponentiate = T)


# Mixed coxph does not converge.. 
me_combined <- coxme(Surv(Lifetime, CompleteData) ~ 
                      Branchtype * Location + Genotype + (1|Culture),
                    weights = Prob,
                    data = total)

summary(me_combined)


## Plotting setup ----
Branchtheme <- theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none", #legend.justification = "top", 
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()
        #panel.grid.major.x = element_blank(),
        #text = element_text(family = "Source Sans Pro")
  )



## Colors
# for Genotype & Location
sci_pal = "batlow"

# For Branchtypes
sci_pal2 <- "roma"



## Plotting functions
# formations over time
formations_plot <- function(data) {
  ggplot({{data}}, 
         aes(x = Bins, y = Formed/CellSum, fill = Branchtype)) +
    Branchtheme + 
    
    geom_col() + 
    facet_grid(Location ~ Genotype) +
    
    theme(panel.grid.major.x = element_blank(),
          legend.position = c(0.5, 0.49), legend.direction = "horizontal",
          legend.key.size = unit(10, "pt"),
          #legend.margin = margin(0,1,0,1, unit = "cm")
    ) +
    scale_x_discrete(labels = c("", "", "6h", "", "", "12 h", "", "", "18 h", "", "", "24 h")) +
    scale_y_continuous(limits = c(0, 0.3)) +
    labs(
      x = "Time",
      y = "Branch initiations per cell",
      title = "",
      fill = ""
    ) 
}

# Cumulative branches

cumulative_plot <- function(data) {
ggplot({{data}}, 
       aes(x = Bins, y = Accumulated_normalized, fill = Branchtype)) +
    Branchtheme + 
  
  geom_col() + 
  facet_grid(Location ~ Genotype) +
  
  theme(panel.grid.major.x = element_blank(),
        legend.position = c(0.5, 0.49), legend.direction = "horizontal",
        legend.key.size = unit(10, "pt"),
        #legend.margin = margin(0,1,0,1, unit = "cm")
        ) +
    scale_x_discrete(labels = c("", "", "6h", "", "", "12 h", "", "", "18 h", "", "", "24 h")) +
    scale_y_continuous(limits = c(0, 1.4)) +
    labs(
    x = "Time",
    y = "Branches per cell",
    title = "",
    fill = ""
  ) 
}

# Events per group Scatter
events_scatter <- function(data) {
ggplot({{data}}, aes(x = Genotype, y = norm_events, color = Genotype)) + 
  geom_quasirandom(position = "dodge", alpha = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", color = "black", fill = NA, width = 0.7, size = 0.6) + 
  stat_summary(fun.data = "mean_se", geom = "errorbar", color = "black", width = 0.4, size = 0.3) + 
  
  facet_wrap(~Branchtype, nrow = 1, strip.position = "bottom") +
  Branchtheme + 
  theme(panel.grid.major.x = element_blank(), axis.text.x = element_blank(), # remove direct Genotype labels
        legend.position = "bottom", legend.box.spacing = unit(-0.9, 'cm') # remove extra space between legend & plot
        ) +
  labs(
    x = "",
    y = "Branch formations\n per cell in 24h",
    color = ""
  )
}  


## Survival plotting with CI
survival_CI <- function(data) {
  ggplot(data = {{data}}, aes(x = time, color = set, fill = set)) + 
    geom_step(aes(y = estimate), size = 1) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, color = NA) +   
    Branchtheme +
    scale_y_continuous(labels = scales::percent, 
                       limits = c(0,1),
                       #expand = expansion(mult = c(0.01, 0.1))
    ) +
    scale_x_continuous(breaks = c(0, 6, 12, 18, 24), 
                       labels = c("0 h", "6 h", "12 h", "18 h", "24 h"),
                       #limits = c(0,24),
                       #expand = expansion(mult = c(0.01, 0.05))
    ) +
    
    labs(x = "Follow up time",
         y = " Survival",
         color = "",
         fill = "")  
}



#Selfmade_facet_grid

CI_data <- function(loc, btype) {

df_AxFil  <- data.frame(Genotype = factor(c("WT", "KO"), levels=levels(total$Genotype)),
                           Location = factor(c(loc, loc), levels = levels(total$Location)),
                           Branchtype = factor(c(btype, btype), levels = levels(total$Branchtype))) 

CI_AxFil <- tidy(survfit(m_combined, newdata = df_AxFil)) %>% 
  pivot_longer(!time:n.censor,
               names_to = c(".value", "set"),
               names_pattern = "(.+).(.+)$") %>% 
  mutate(set = fct_recode(set, WT = "1", KO = "2")) # rename to WT & KO

survival_CI(CI_AxFil) +  
    scale_color_scico_d(palette = sci_pal, begin = 0, end = 0.7) + 
    scale_fill_scico_d(palette = sci_pal, begin = 0, end = 0.7) + 
 
  # only plot Facet groups on top-row 
  if_else(loc == "Axon", list(labs(subtitle = btype)), list(labs())) + 
  # only plot axes on outsides of plot
    labs(x = ifelse(loc == "Neurite", "Follow-up time", "" ),
         y = ifelse(btype == "Filopodium", "Survival", "" )) +
  
  if_else(btype %in% c("Mixed", "Lamellipodium", "Splitting"), list(theme(axis.text.y = element_blank())), list(theme())) +
  if_else(loc == "Axon", list(theme(axis.text.x = element_blank())), list(theme())) +
  
  theme(plot.subtitle = element_text(hjust = 0.5), plot.margin = margin(0,0,0,0,"pt"))
}





### Plotting ----

# A) Cumulative events by type & Location
Cumulative <- cumulative_branches(total_filt)

(A0 <- formations_plot(Cumulative)+
    scale_color_scico_d(palette = sci_pal2, begin = 0, end = 1) + 
    scale_fill_scico_d(palette = sci_pal2, begin = 0, end = 1))

(A <- cumulative_plot(Cumulative) +
    scale_color_scico_d(palette = sci_pal2, begin = 0, end = 1) + 
    scale_fill_scico_d(palette = sci_pal2, begin = 0, end = 1))


# B) Branching events by type

B <- events_scatter(count_exp) +
  stat_pvalue_manual(m0_combined, label = "p.signif", 
                     y = 2.2, tip.length = 0, hide.ns = T) +
  scale_color_scico_d(palette = sci_pal, begin = 0, end = 0.7) + 
  scale_fill_scico_d(palette = sci_pal, begin = 0, end = 0.7)



# C) DAG 
C <- ggdraw() + draw_image(file.path("figures", "figure_2","combined_DAG.png"), scale = 0.85)

C_up <- ggdraw() + draw_image(file.path("figures", "figure_2","combined_DAG_vert.png"), scale = 0.7)



## Survival analysis
# All plots
AxFil <- CI_data("Axon", "Filopodium")
AxMix <- CI_data("Axon", "Mixed")
AxLam <- CI_data("Axon", "Lamellipodium")
AxSpl <- CI_data("Axon", "Splitting")
NeuFil <- CI_data("Neurite", "Filopodium")
NeuMix <- CI_data("Neurite", "Mixed")
NeuLam <- CI_data("Neurite", "Lamellipodium")
NeuSpl <- CI_data("Neurite", "Splitting") +
            theme(legend.position = c(0.85, 0.9), legend.direction = "vertical",
            legend.key.size = unit(10, "pt"))

# combine to one facet using cowplot
D1 <- plot_grid(AxFil, AxMix, AxLam, AxSpl, scale = 1, 
               nrow = 1, align = "h", axis = "bt", rel_widths = c(1.15, 1,1,1))

D2 <- plot_grid(NeuFil, NeuMix, NeuLam, NeuSpl, scale = 1, 
               nrow = 1, align = "h", axis = "bt", rel_widths = c(1.15, 1,1,1))


D <- plot_grid(D1, D2, nrow = 2, align = "v", axis = "bt", rel_heights = c(1.1,1)) 
# legend, Axon/Neurite
D



## Merge into final figure

right <- plot_grid(B, C, nrow = 2, rel_heights = c(2,1), labels = c("B", "C")) 

top <- plot_grid(A, right, ncol = 2,  rel_widths = c(1.05, 1), scale = 0.98) 

spacer <- ggplot() + Branchtheme

Fig_2 <- plot_grid(top, spacer, D, nrow = 3, labels = c("A", "", "D"), rel_heights = c(1.1,0.05, 1)) 



ggsave(file.path("figures", "figure_2", "Fig_2.png"), Fig_2, device = "png", scale = 1, width = 210, height = 200, units = "mm" )

ggsave(file.path("figures", "figure_2", "Fig_2.pdf"), Fig_2, device = "pdf", scale = 1, width = 210, height = 210, units = "mm" )


## second option

top0 <- plot_grid(A0, A, ncol = 2, labels = c("A", "B"))
mid <- plot_grid(B, C_up, ncol = 2, labels = c("C", "D"), rel_widths = c(1.5,1))


Fig_2b <- plot_grid(top0, mid, spacer, D, nrow = 4, labels = c("", "", "E"), rel_heights = c(1.1,1, 0.05, 1)) 



ggsave(file.path("figures", "figure_2", "Fig_2b.png"), Fig_2b, device = "png", scale = 1, width = 210, height = 240, units = "mm" )

