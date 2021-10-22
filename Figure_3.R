## Figure 3 ----
## now figure 4

library(tidyverse)
library(survival)
library(broom)
library(ggfortify)
library(coxme)

library(ggbeeswarm)
library(scico)
library(scales)
library(ggpubr)
library(cowplot)




source(file.path("functions_dataset_2" , "clean_data.R"))
source(file.path("functions_dataset_2" , "cumulative_branches.R"))


# load data
total <- clean_data()


## Statistics ----


### Survival analysis (total treatment effect)

## Treatment --
m_treatment <-  coxph(Surv(Lifetime, CompleteData) ~ Treatment,
                     weights = Prob,
                     data = total)

# Diagnostics
diag_surv_treatment<- cox.zph(m_treatment)
survminer::ggcoxzph(diag_surv_treatment) # Schoenfeld test
survminer::ggcoxdiagnostics(m_treatment, type = "deviance", linear.predictions = FALSE) # 

# Results table
tidy(m_treatment, exponentiate = T)




## Branching events by type and Treatment
count_exp <- total %>% 
  group_by(Branchtype, Culture, Treatment, CellSum) %>% 
  summarise(events = n(),
            Neurons = sum(unique(Neurons))) %>% 
  mutate(norm_events = events/Neurons,
         tf_norm_events = sqrt(norm_events)) %>% # to reach normality
  ungroup() %>% group_by(Culture) %>% 
  mutate(norm_norm_events = norm_events/mean(norm_events)) %>% ungroup

# Testing model assumptions
m0_comb <- aov(tf_norm_events ~ Treatment, data = count_exp)
autoplot(m0_comb, which = 1:4, ncol = 2, label.size = 3, colour = "Treatment")
shapiro.test(m0_comb$residuals)
car::leveneTest(m0_comb)

# Welchs's t-test 
m0_combined <- compare_means(tf_norm_events ~ Treatment, group.by = "Branchtype", ref.group = "Control", 
                             data = count_exp, method = "t.test", p.adjust.method = "holm")
# no comparison is significant, but effects are in accordance to literature



## Survival analysis (treatment effect not mediated by precursor types)
m_combined <- coxph(Surv(Lifetime, CompleteData) ~ 
                      Branchtype + Location + Treatment,
                    weights = Prob,
                    data = total)

# mixed effect version
me_combined <- coxme(Surv(Lifetime, CompleteData) ~ 
                       Branchtype + Location + Treatment + (1|Culture),
                     weights = Prob,
                     data = total)

# Diagnostics
diag_surv_combined<- cox.zph(m_combined)
survminer::ggcoxzph(diag_surv_combined) # Schoenfeld test
survminer::ggcoxdiagnostics(m_combined, type = "deviance", linear.predictions = FALSE) # 

# Results table
tidy(m_combined, exponentiate = T)

# hazard estimates do not change with mixed model, the p-values are lower but still significant at precursor & neurite comparisons
summary(me_combined)





## Plotting setup ----
Branchtheme <- theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none", 
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()
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
    facet_grid(Location ~ Treatment) +
    
    theme(panel.grid.major.x = element_blank(),
          legend.position = c(0.5, 0.49), legend.direction = "horizontal",
          legend.key.size = unit(10, "pt"),
          plot.title = element_text(hjust = 0.5)
    ) +
    scale_x_discrete(labels = c("", "", "6h", "", "", "12 h", "", "", "18 h", "", "", "24 h", "", "", "30 h")) +
    scale_y_continuous(limits = c(0, 1)) +
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
    facet_grid(Location ~ Treatment) +
    
    theme(panel.grid.major.x = element_blank(),
          legend.position = c(0.5, 0.49), legend.direction = "horizontal",
          legend.key.size = unit(10, "pt"),
          plot.title = element_text(hjust = 0.5)
    ) +
    scale_x_discrete(labels = c("", "", "6h", "", "", "12 h", "", "", "18 h", "", "", "24 h", "", "", "30 h")) +
    scale_y_continuous(limits = c(0, 6)) +
    labs(
      x = "Time",
      y = "Branches per cell",
      title = "",
      fill = ""
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
   # facet_wrap(~Treatment, ncol = 3) +
    
    scale_x_continuous(breaks = c(0, 6, 12, 18, 24, 30), 
                       labels = c("0 h", "6 h", "12 h", "18 h", "24 h", "30 h")) +
    scale_y_continuous(breaks = c(0, 6, 12, 18, 24, 30), 
                       labels = c("0 h", "6 h", "12 h", "18 h", "24 h", "30 h"))
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
    scale_x_continuous(breaks = c(0, 6, 12, 18, 24, 30), 
                       labels = c("0 h", "6 h", "12 h", "18 h", "24 h", "30 h"),
                       #limits = c(0,24),
                       #expand = expansion(mult = c(0.01, 0.05))
    ) +
    
    labs(x = "Follow up time",
         y = " Survival") 
}




# Events per group Scatter
events_scatter <- function(data) {
  ggplot({{data}}, aes(x = Treatment, y = norm_norm_events, color = Treatment)) + 
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
    scale_x_continuous(breaks = c(0, 6, 12, 18, 24, 30), 
                       labels = c("0 h", "6 h", "12 h", "18 h", "24 h", "30 h"),
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
  
  df_AxFil  <- data.frame(Treatment = factor(c("Control", "Netrin", "FGF"), levels=levels(total$Treatment)),
                          Location = factor(c(loc, loc, loc), levels = levels(total$Location)),
                          Branchtype = factor(c(btype, btype, btype), levels = levels(total$Branchtype))) 
  
  CI_AxFil <- tidy(survfit(m_combined, newdata = df_AxFil)) %>% 
    pivot_longer(!time:n.censor,
                 names_to = c(".value", "set"),
                 names_pattern = "(.+).(.+)$") %>% 
    mutate(set = fct_recode(set, "Control" = "1", "Netrin" = "2", "FGF" = "3")) # rename to WT & KO
  
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

# A) Cumulative events by type & Location (excluding "unclear")
Cumulative <- cumulative_branches(total)

(A0 <- formations_plot(filter(Cumulative, Location != "Unclear"))+
    labs(title = "Initiations of branches over time") +
    scale_color_scico_d(palette = sci_pal2, begin = 0, end = 1) + 
    scale_fill_scico_d(palette = sci_pal2, begin = 0, end = 1))

(A <- cumulative_plot(filter(Cumulative, Location != "Unclear")) +
    labs(title = "Accumulation of branches over time") +
    scale_color_scico_d(palette = sci_pal2, begin = 0, end = 1) + 
    scale_fill_scico_d(palette = sci_pal2, begin = 0, end = 1))



# total overview axon vs. dendrite by Branchtypes and Genotype
B <- lifetime_scatter(total, Treatment) +
  scale_color_scico_d(palette = sci_pal, begin = 0, end = 0.7)





#C: CoxPH + Confidence intervals
df_Treatment  <- data.frame(Treatment = factor(c("Control", "Netrin", "FGF"), levels=levels(total$Treatment)))

CI_treatment <- tidy(survfit(m_treatment, newdata=df_Treatment)) %>% 
  pivot_longer(!time:n.censor,
               names_to = c(".value", "set"),
               names_pattern = "(.+).(.+)$")


(C <- survival_CI(CI_treatment) +  
    scale_color_scico_d(palette = sci_pal, begin = 0, end = 0.7) + 
    scale_fill_scico_d(palette = sci_pal, begin = 0, end = 0.7)
)



#D) Forest plot of individual model
D <- survminer::ggforest(m_treatment, data = total)





# E) DAG placeholder
spacer <- ggplot() + Branchtheme


# F) Forest plot of final model
EF <- survminer::ggforest(m_combined, data = total)


## G) multifactorial survival analysis
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
G1 <- plot_grid(AxFil, AxMix, AxLam, AxSpl, scale = 1, 
                nrow = 1, align = "h", axis = "bt", rel_widths = c(1.15, 1,1,1))

G2 <- plot_grid(NeuFil, NeuMix, NeuLam, NeuSpl, scale = 1, 
                nrow = 1, align = "h", axis = "bt", rel_widths = c(1.15, 1,1,1))


G <- plot_grid(G1, G2, nrow = 2, align = "v", axis = "bt", rel_heights = c(1.1,1)) 



## Merge into final figure
top <- plot_grid(A0, A, ncol = 2, labels = c("A", "B"))
mid <- plot_grid(B, C, D, spacer, ncol = 4, labels = c("C", "D", "E", "F"), rel_widths = c(1,1, 1.5, 0.5))
bottom <- plot_grid(spacer, EF, ncol = 2, labels = c("G", "H", rel_widths = c(0.5, 2)))

Fig_3 <- plot_grid(top, mid, bottom, G, nrow = 4, labels = c("", "","", "H"), rel_heights = c(1.5,1, 1,1.5)) 



ggsave(file.path("figures", "figure_3", "Fig_3_raw.png"), Fig_3, device = "png", scale = 1, width = 210, height = 240, units = "mm" )
ggsave(file.path("figures", "figure_3", "Fig_3.pdf"), Fig_3, device = "pdf", scale = 1, width = 210, height = 240, units = "mm" )



# Supplementary
# Branching events by type

(Sup1 <- events_scatter(count_exp) +
    stat_pvalue_manual(m0_combined, label = "p.signif", 
                       y = 2.2, tip.length = 0, hide.ns = T) +
    labs(title = "Formation of branches by precursor") +
    scale_color_scico_d(palette = sci_pal, begin = 0, end = 0.7) + 
    scale_fill_scico_d(palette = sci_pal, begin = 0, end = 0.7)) +
  scale_y_continuous(limits = c(0,2.2))

ggsave(file.path("figures", "supplement", "Sup_2B_raw.png"), Sup1, device = "png", scale = 1, width = 160, height = 80, units = "mm" )
ggsave(file.path("figures", "supplement",  "Sup_2B.pdf"), Sup1, device = "pdf", scale = 1, width = 160, height = 80, units = "mm" )




