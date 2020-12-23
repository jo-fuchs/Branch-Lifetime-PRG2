## Figure 1 ----

library(tidyverse)
library(survival)
library(broom)
library(ggfortify)

library(ggbeeswarm)
library(scico)
library(scales)
library(ggpubr)
library(cowplot)

source("clean_data.R")


# load data
total <- clean_data()



## Statistics ----

# Summarize Lifetime data by experiment and specified additional column
mean_lifetime_by_experiment <- function(total_data, by_column) {
  total_data %>% 
    group_by(Culture, {{by_column}}) %>% 
    mutate(Lifetime_trans = 1/Lifetime) %>%  # to adjust for right-skewed data
    summarise(Lifetime = mean(Lifetime_trans)) %>% 
    ungroup()
}


## raw lifetimes not accounting for right-censoring of data --

## A) Genotype --
# Mean lifetime per experiment
genotype_exp <- mean_lifetime_by_experiment(total, Genotype)


# Testing model assumptions
m0_geno <- aov(Lifetime ~ Genotype, data = genotype_exp)
autoplot(m0_geno, which = 1:4, ncol = 2, label.size = 3, colour = "Genotype")


# Testing with Welch's test (non-equal variances)
m0_genotype <- compare_means(Lifetime ~ Genotype, 
                             data = genotype_exp, method = "t.test")



## E) Location --
location_exp <- mean_lifetime_by_experiment(total, Location)


# Testing model assumptions
m0_loc <- aov(Lifetime ~ Location, data = location_exp)
autoplot(m0_loc, which = 1:4, ncol = 2, label.size = 3, colour = "Location")


# Testing with Wilcox test (non-normality & non-equal variances)
m0_location <- compare_means(Lifetime ~ Location, 
                               data = filter(location_exp, Location != "Unclear"), method = "t.test")



## I) Branchtype --
branchtype_exp <- mean_lifetime_by_experiment(total, Branchtype)


# Testing model assumptions
m0_btype <- aov(Lifetime ~ Branchtype, data = branchtype_exp)
autoplot(m0_btype, which = 1:4, ncol = 2, label.size = 3, colour = "Branchtype")


# Post Hoc comparisons: pairwise Welch's test (non-equal variances)
m0_branchtype <- compare_means(Lifetime ~ Branchtype, 
                               data = branchtype_exp, method = "t.test", 
                               p.adjust.method = "holm")




### Survival analysis

## c) Genotype --
m_genotype <-  coxph(Surv(Lifetime, CompleteData) ~ Genotype,
                   weights = Prob,
                    #cluster = Culture,
                     data = total)
# Diagnostics
diag_surv_genotype<- cox.zph(m_genotype)
survminer::ggcoxzph(diag_surv_genotype) # Schoenfeld test
survminer::ggcoxdiagnostics(m_genotype, type = "deviance", linear.predictions = FALSE) # 

# Results table
tidy(m_genotype, exponentiate = T)



## G) Location --
m_location <-  coxph(Surv(Lifetime, CompleteData) ~ Location ,
                     weights = Prob,
                     data = total)

# Diagnostics
diag_surv_location<- cox.zph(m_location)
survminer::ggcoxzph(diag_surv_location) # Schoenfeld test
survminer::ggcoxdiagnostics(m_location, type = "deviance", linear.predictions = FALSE) # 


# Results table
tidy(m_location, exponentiate = T)



## K) Branchtype --
m_branchtype <-  coxph(Surv(Lifetime, CompleteData) ~ Branchtype,
                      weights = Prob,
                      data = total)

# Diagnostics
diag_surv_branchtype<- cox.zph(m_branchtype)
survminer::ggcoxzph(diag_surv_branchtype) # Schoenfeld test
survminer::ggcoxdiagnostics(m_branchtype, type = "deviance", linear.predictions = FALSE) # 

# Results table
tidy(m_branchtype, exponentiate = T)



## Plotting setup ----
Branchtheme <- theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        plot.title.position = "plot", plot.caption.position = "plot",
        legend.position = "none", legend.justification = "top", 
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()
        #panel.grid.major.x = element_blank(),
        #text = element_text(family = "Source Sans Pro")
  )



## Colors from scico-package --
# Crameri, F. (2018). Scientific colour maps. Zenodo. 
# http://doi.org/10.5281/zenodo.1243862

# for Genotype & Location
sci_pal = "batlow"

# For Branchtypes
sci_pal2 <- "roma"




## custom plotting functions --
# Lifetime beeswarm superplot 
# Lord, SL et al. (2020) SuperPlots: Communicating reproducibility and variability in cell biology. 
# J Cell Biol 1 June 2020; doi: https://doi.org/10.1083/jcb.202001064)

lifetime_swarm <- function(per_exp, total, X_axis) {
  ggplot(per_exp, aes(x = {{X_axis}}, y = 1/Lifetime, color = {{X_axis}})) +  # 1/Lifetime to reverse normalization for stats
  geom_point() + 

  stat_summary(fun = mean, geom = "crossbar", color = "black", width = 0.6, size = 0.3) + 
  stat_summary(fun.data = "mean_se", geom = "errorbar", color = "black", width = 0.4, size = 0.3) + 
    
  geom_quasirandom(aes(y = Lifetime), data = {{total}}, alpha = 0.01, shape = 16) + # non-collapsing
  Branchtheme + theme(panel.grid.major.x = element_blank(), 
  ) + 
  
  scale_y_continuous(breaks = c(0, 6, 12, 18, 24), 
                         labels = c("0 h", "6 h", "12 h", "18 h", "24 h")) +
    
  labs(
    x = "",
    y = "Lifetime per branch"
  )
  }


# Scatterplot highlighting each indiviual branch event, colored by specified column
lifetime_scatter <- function(total, Group) {
  ggplot({{total}}, aes(x = Collapse, y = Formation, color = {{Group}})) + 
  geom_point(alpha = 0.1, shape = 16) + 
  
  Branchtheme +
  labs(
    x = "Timepoint of Collapse",
    y = "Timepoint of Formation"
  ) +
    
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24), 
                       labels = c("0 h", "6 h", "12 h", "18 h", "24 h")) +
  scale_y_continuous(breaks = c(0, 6, 12, 18, 24), 
                     labels = c("0 h", "6 h", "12 h", "18 h", "24 h"))
}



## Survival curve with CI
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
         y = " Survival") 
}




## A: WT vs. KO -----

# lifetime beeswarm

A1 <- lifetime_swarm(genotype_exp, total, Genotype) +
  
  stat_pvalue_manual(m0_genotype, label = "p.signif", 
                     y = 23, tip.length = 0) +
  
  scale_color_scico_d(palette = sci_pal, begin = 0, end = 0.7)



# total overview axon vs. dendrite by Branchtypes and Genotype
A2 <- lifetime_scatter(total, Genotype) +
  #scale_color_viridis_d()
  scale_color_scico_d(palette = sci_pal, begin = 0, end = 0.7)





#A3c: CoxPH + Confidence intervals
df_Genotype  <- data.frame(Genotype = factor(c("WT", "KO"), levels=levels(total$Genotype)))

CI_genotype <- tidy(survfit(m_genotype, newdata=df_Genotype)) %>% 
  pivot_longer(!time:n.censor,
               names_to = c(".value", "set"),
               names_pattern = "(.+).(.+)$")


(A3c <- survival_CI(CI_genotype) +  
  scale_color_scico_d(palette = sci_pal, begin = 0, end = 0.7) + 
  scale_fill_scico_d(palette = sci_pal, begin = 0, end = 0.7)
)


## Location ----

(B1 <- lifetime_swarm(filter(location_exp, Location != "Unclear"), 
               filter(total, Location != "Unclear"), 
               Location) +
   stat_pvalue_manual(m0_location, label = "p.signif", 
                      y = 23, tip.length = 0) +
  scale_color_scico_d(palette = sci_pal, begin = 0.25, end = 0.5))

# total overview axon vs. dendrite by Branchtypes and Genotype
B2 <- lifetime_scatter(filter(total, Location != "Unclear"), Location) +
  #scale_color_viridis_d()
  scale_color_scico_d(palette = sci_pal, begin = 0.25, end = 0.5)




#B3c: CoxPH + CIst
df_Location  <- data.frame(Location = factor(c("Axon", "Neurite"), levels=levels(total$Location)))

CI_location <- tidy(survfit(m_location, newdata=df_Location)) %>% 
  pivot_longer(!time:n.censor,
               names_to = c(".value", "set"),
               names_pattern = "(.+).(.+)$")


(B3c <- survival_CI(CI_location) +  
    scale_color_scico_d(palette = sci_pal, begin = 0.25, end = 0.5) + 
    scale_fill_scico_d(palette = sci_pal, begin = 0.25, end = 0.5)
)



## Branchtype ----

(C1 <- lifetime_swarm(branchtype_exp, 
               total, 
               Branchtype) +
  stat_pvalue_manual(m0_branchtype, label = "p.signif", 
                     y = c(21, 22.5, 24, 18, 19.5 ), tip.length = 0, hide.ns = T) +

    scale_color_scico_d(palette = sci_pal2, begin = 0, end = 1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
 )

# total overview axon vs. dendrite by Branchtypes and Genotype
(C2 <- lifetime_scatter(total, Branchtype) +
  #scale_color_viridis_d()
  scale_color_scico_d(palette = sci_pal2, begin = 0, end = 1))


#B3c: CoxPH + CIst
df_type  <- data.frame(Branchtype = factor(c("Filopodium", "Mixed", "Lamellipodium", "Splitting"), 
                                           levels=levels(total$Branchtype)),
                       set = as.character(1: length(levels(total$Branchtype))))

CI_type <- tidy(survfit(m_branchtype, newdata=df_type)) %>% 
  pivot_longer(!time:n.censor,
               names_to = c(".value", "set"),
               names_pattern = "(.+).(.+)$") %>% 
  left_join(df_type)


(C3c <- survival_CI(CI_type) +  
    scale_color_scico_d(palette = sci_pal2, begin = 0, end = 1) + 
    scale_fill_scico_d(palette = sci_pal2, begin = 0, end = 1))



# Load DAG .png files
A4 <- ggdraw() + draw_image(file.path("figures", "figure_1","PRG2.png"), scale = 0.9)
B4 <- ggdraw() + draw_image(file.path("figures", "figure_1","Location.png"), scale = 0.9)
C4 <- ggdraw() + draw_image(file.path("figures", "figure_1","Type.png"), scale = 0.9)


## Merge it in one figure ----
A <- plot_grid(A1, A2, A3c, A4, scale = 0.9, labels = c("A", "B", "C", "D"), 
               rel_widths = c(1.5, 2, 2, 1), nrow = 1, align = "h", axis = "bt")

B <- plot_grid(B1, B2, B3c, B4, scale = 0.9, labels = c("E", "F", "G", "H"), 
               rel_widths = c(1.5, 2, 2, 1), nrow = 1, align = "h", axis = "bt")

C <- plot_grid(C1, C2, C3c, C4, scale = 0.9, labels = c("I", "J", "K", "L"), 
               rel_widths = c(1.5, 2, 2, 1), nrow = 1, align = "h", axis = "bt")

Fig_1 <- plot_grid(A, B, C, nrow = 3, rel_heights = c(1, 1, 1.2))



ggsave(file.path("figures", "figure_1", "Fig_1.png"), Fig_1, device = "png", scale = 1, width = 210, height = 210, units = "mm" )

ggsave(file.path("figures", "figure_1", "Fig_1.pdf"), Fig_1, device = "pdf", scale = 1, width = 210, height = 210, units = "mm" )



