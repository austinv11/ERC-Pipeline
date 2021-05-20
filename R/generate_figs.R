# Code to generate supplemental figures

library(ggplot2)
library(readxl)
library(writexl)
library(nVennR)
library(ggvenn)
library(dplyr)
library(gridExtra)
library(ape)
library(picante)
library(Cairo)
library(EnvStats)
library(cowplot)
library(patchwork)

options(digits = 3)
CairoWin()

generic_corr_plot <- function() {
  x <- 25:125
  x <- jitter(sample(x, 25, replace=TRUE), amount = 25)
  
  ycorr <- jitter(sapply(x, function(v) {v}), amount = 25)
  yrand <- jitter(sample(25:125, 25, replace = TRUE), amount = 25)
  
  blank_theme <- theme(axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank(),
                       panel.border = element_blank(),
                       plot.title = element_text(hjust = 0.5),
                       axis.line = element_line())
  
  plot1 <- ggplot(data.frame(x=x, y=ycorr), aes(x=x, y=y)) + geom_point() + xlab("Protein 1 Rate") + theme_bw() + blank_theme + ylab("Protein 2 Rate") + ggtitle("Strong Positive Correlation") + geom_smooth(color = "blue", size = 1, method = "lm", se = FALSE)
  plot2 <- ggplot(data.frame(x=x, y=yrand), aes(x=x, y=y)) + geom_point() + xlab("Protein 1 Rate") + theme_bw() + blank_theme + ylab("Protein 2 Rate") + ggtitle("No Correlation") + geom_hline(color = "blue", size = 1, aes(yintercept = 50)) + ylim(0, 100) 
  
  return(grid.arrange(plot1, plot2, ncol=2))
}

plot <- generic_corr_plot()
ggsave(file="generic_corr_plot.png", plot=plot, type="cairo-png")


plot_ace2_xy <- function(protein) {
  full_df <- read_excel(paste("fig_data", "ace2_xys",
                            paste(protein, "_ace2_data.xlsx", sep = ""),
                            sep = "/"), sheet = "Full Data")
  mya_df <- read_excel(paste("fig_data", "ace2_xys",
                            paste(protein, "_ace2_data.xlsx", sep = ""),
                            sep = "/"), sheet = "30MYA")
  
  full_df <- full_df %>% mutate(
    order = case_when(
      Order == "Primates" ~ "Primates",
      Order == "Rodentia" ~ "Rodentia",
      Order == "Chiroptera" ~ "Chiroptera",
      Order == "Artiodactyla" ~ "Artiodactyla",
      Order == "Carnivora" ~ "Carnivora",
      T ~ "Other"
    )
  )
  
  mya_df <- mya_df %>% mutate(
    order = case_when(
      Order == "Primates" ~ "Primates",
      Order == "Rodentia" ~ "Rodentia",
      Order == "Chiroptera" ~ "Chiroptera",
      Order == "Artiodactyla" ~ "Artiodactyla",
      Order == "Carnivora" ~ "Carnivora",
      T ~ "Other"
    )
  )
  
  plot1 <- ggplot(full_df, aes(x = `ACE2 Rate Rank`, y = `Prot Rate Rank`, color=order)) + 
    geom_point()+ scale_color_manual("Order", values = c("purple", "blue", "#00BA38", "black", "red", "#E58700")) +
    ggtitle(paste("Original Branch Data (ACE2 vs ", protein, ")", sep = "")) + ylab(paste(protein, "Rate Rank")) + geom_smooth(color = "black", size = 1, method = "lm", se = FALSE) +
    theme_half_open() +
    background_grid()
  
  plot2 <- ggplot(mya_df, aes(x = `ACE2 Rate Rank`, y = `Prot Rate Rank`, color=order))  + 
    geom_point()+ scale_color_manual("Order", values = c("purple", "blue", "#00BA38", "black", "red", "#E58700")) +
    ggtitle(paste("30MY Branch Data (ACE2 vs ", protein, ")", sep = "")) + ylab(paste(protein, "Rate Rank")) + geom_smooth(color = "black", size = 1, method = "lm", se = FALSE) +
    theme_half_open() +
    background_grid()

  return(plot1 | plot2)
}

plot_time_xy <- function(protein, return_indiv=FALSE) {
  full_df <- read_excel(paste("fig_data", "btime_xys",
                              paste(protein, "_btime_data.xlsx", sep = ""),
                              sep = "/"), sheet = "Full Data")
  mya_df <- read_excel(paste("fig_data", "btime_xys",
                             paste(protein, "_btime_data.xlsx", sep = ""),
                             sep = "/"), sheet = "30MYA")
  stats_df <- read_excel(paste("fig_data", "btime_xys",
                               paste(protein, "_btime_data.xlsx", sep = ""),
                               sep = "/"), sheet = "Summary")
  
  full_df <- full_df %>% mutate(
    order = case_when(
      Order == "Primates" ~ "Primates",
      Order == "Rodentia" ~ "Rodentia",
      Order == "Chiroptera" ~ "Chiroptera",
      Order == "Artiodactyla" ~ "Artiodactyla",
      Order == "Carnivora" ~ "Carnivora",
      T ~ "Other"
    )
  )
  
  mya_df <- mya_df %>% mutate(
    order = case_when(
      Order == "Primates" ~ "Primates",
      Order == "Rodentia" ~ "Rodentia",
      Order == "Chiroptera" ~ "Chiroptera",
      Order == "Artiodactyla" ~ "Artiodactyla",
      Order == "Carnivora" ~ "Carnivora",
      T ~ "Other"
    )
  )
  
  full_label <- substitute(
    paste(rho, " = ", estimate, ", P = ", pvalue),
    list(estimate = signif(stats_df[stats_df$Treatment == 'Full Taxa',c(2, 3)][[1]], 2), pvalue = signif(stats_df[stats_df$Treatment == 'Full Taxa',c(2, 3)][[2]], 2))
  )
  my_label <- substitute(
    paste(rho, " = ", estimate, ", P = ", pvalue),
    list(estimate = signif(stats_df[stats_df$Treatment == '30MY Taxa',c(2, 3)][[1]], 2), pvalue = signif(stats_df[stats_df$Treatment == '30MY Taxa',c(2, 3)][[2]], 2))
  )
  
  plot1 <- ggplot(full_df, aes(x = `Time`, y = `Rate`, color=order)) + 
    geom_point()+ scale_color_manual("Order", values = c("purple", "blue", "#00BA38", "black", "red", "#E58700")) +
    ggtitle(paste("Original Branch Data (", protein, " Rate)", sep = "")) + ylab(paste(protein, "Rate")) +
    xlim(0, 200) + ylim(0, max(full_df$Rate)*1.25) + xlab("Terminal Branch Time (MY)") +
    geom_smooth(color = "black", size = 0.75, method = "lm", se = FALSE, linetype="dashed") +
    draw_label(full_label, x = 100, y = max(full_df$Rate)*1.25, size = 19, fontface="bold") +
    theme_half_open(font_size=19, rel_large=0.8) 
  plot2 <- ggplot(mya_df, aes(x = `Time`, y = `Rate`, color=order))  + 
    geom_point()+ scale_color_manual("Order", values = c("purple", "blue", "#00BA38", "black", "red", "#E58700")) +
    ggtitle(paste("30MY Branch Data (", protein, " Rate)", sep = "")) + ylab(paste(protein, "Rate")) +
    xlim(0, 200) + ylim(0, max(mya_df$Rate)*1.25) + xlab("Terminal Branch Time (MY)") +
    geom_smooth(color = "black", size = 0.75, method = "lm", se = FALSE, linetype="dashed") +
    draw_label(my_label, x = 100, y = max(mya_df$Rate)*1.25, size = 19, fontface="bold") +
    theme_half_open(font_size=19, rel_large=0.8) 
  
  if (return_indiv) {
    return(list(full.plot=plot1, my.plot=plot2))
  } else {
    return(plot1 | plot2)
  }
}


time_regression_for_prots <- function(protein) {
  full_df <- read_excel(paste("fig_data", "btime_xys",
                              paste(protein, "_btime_data.xlsx", sep = ""),
                              sep = "/"), sheet = "Full Data")
  mya_df <- read_excel(paste("fig_data", "btime_xys",
                             paste(protein, "_btime_data.xlsx", sep = ""),
                             sep = "/"), sheet = "30MYA")
  
  mod <- lm(`Rate Rank` ~ `Time Rank`, data = full_df)
  fullaov <- aov(`Rate Rank` ~ `Time Rank`, data = full_df)
  
  
  mod2 <- lm(`Rate Rank` ~ `Time Rank`, data = mya_df)
  myaov <- aov(`Rate Rank` ~ `Time Rank`, data = mya_df)
  
  return(list(full.lm = summary(mod), my.lm = summary(mod2), 
              full.aov = summary(fullaov), my.aov = summary(myaov)))
}


pcortests <- function(protein) {
  full_df <- read_excel(paste("fig_data", "ace2_xys",
                              paste(protein, "_ace2_data.xlsx", sep = ""),
                              sep = "/"), sheet = "Full Data")
  mya_df <- read_excel(paste("fig_data", "ace2_xys",
                             paste(protein, "_ace2_data.xlsx", sep = ""),
                             sep = "/"), sheet = "30MYA")
  
  mod <- lm(`Prot Rate Rank` ~ `Time Rank`, 
            data = full_df)
  
  mod2 <- lm(`ACE2 Rate Rank` ~ `Time Rank`, 
             data = full_df)
  
  combined.mod <- lm(prot_resid ~ ace2_resid, 
                     data=data.frame(prot_resid=resid(mod), ace2_resid=resid(mod2)))
  
  mymod <- lm(`Prot Rate Rank` ~ `Time Rank`, 
              data = mya_df)
  
  mymod2 <- lm(`ACE2 Rate Rank` ~ `Time Rank`, 
               data = mya_df)
  
  mycombined.mod <- lm(prot_resid ~ ace2_resid, 
                       data=data.frame(prot_resid=resid(mymod), ace2_resid=resid(mymod2)))
  
  
  return(list(full.prot.resid.bt.cor = cor.test(resid(mod), full_df$`Time`, method = "spearman"),
              full.ace2.resid.bt.cor = cor.test(resid(mod2), full_df$`Time`, method = "spearman"),
              my.prot.resid.bt.cor = cor.test(resid(mymod), mya_df$`Time`, method = "spearman"),
              my.ace2.resid.bt.cor = cor.test(resid(mymod2), mya_df$`Time`, method = "spearman"),
              full.resid.lm = summary(combined.mod), my.resid.lm=summary(mycombined.mod)))
}

regression_for_prots_no_order <- function(protein) {
  full_df <- read_excel(paste("fig_data", "ace2_xys",
                              paste(protein, "_ace2_data.xlsx", sep = ""),
                              sep = "/"), sheet = "Full Data")
  mya_df <- read_excel(paste("fig_data", "ace2_xys",
                             paste(protein, "_ace2_data.xlsx", sep = ""),
                             sep = "/"), sheet = "30MYA")
  
  mod <- lm(`Prot Rate Rank` ~ `ACE2 Rate Rank` + `Time Rank`, 
            data = full_df)
  
  fullaov <- aov(`Prot Rate Rank` ~ `ACE2 Rate Rank` + `Time Rank`, 
                 data = full_df)
  
  mod2 <- lm(`Prot Rate Rank` ~ `ACE2 Rate Rank` + `Time Rank`, 
             data = mya_df)
  
  myaov <- aov(`Prot Rate Rank` ~ `ACE2 Rate Rank` + `Time Rank`, 
               data = mya_df)
  
  return(list(full.lm=summary(mod), my.lm=summary(mod2), 
              full.aov = summary(fullaov), my.aov = summary(myaov)))
}


regression_for_prots <- function(protein) {
  full_df <- read_excel(paste("fig_data", "ace2_xys",
                              paste(protein, "_ace2_data.xlsx", sep = ""),
                              sep = "/"), sheet = "Full Data")
  mya_df <- read_excel(paste("fig_data", "ace2_xys",
                             paste(protein, "_ace2_data.xlsx", sep = ""),
                             sep = "/"), sheet = "30MYA")
  
  mod <- lm(`Prot Rate Rank` ~ `ACE2 Rate Rank` + Order + `Time Rank`, 
            data = full_df)
  
  fullaov <- aov(`Prot Rate Rank` ~ `ACE2 Rate Rank` + Order + `Time Rank`, 
                 data = full_df)
  
  mod2 <- lm(`Prot Rate Rank` ~ `ACE2 Rate Rank` + Order + `Time Rank`, 
            data = mya_df)
  
  myaov <- aov(`Prot Rate Rank` ~ `ACE2 Rate Rank` + Order + `Time Rank`, 
                 data = mya_df)
  
  return(list(full.lm=summary(mod), my.lm=summary(mod2), 
              full.aov = summary(fullaov), my.aov = summary(myaov)))
}

contrasts_test <- function(protein) {
  full_df <- read_excel(paste("fig_data", "ace2_xys",
                              paste(protein, "_ace2_data.xlsx", sep = ""),
                              sep = "/"), sheet = "Full Data")
  my_df <- read_excel(paste("fig_data", "ace2_xys",
                             paste(protein, "_ace2_data.xlsx", sep = ""),
                             sep = "/"), sheet = "30MYA")
  
  full_tree <- read.tree(paste("fig_data", "ace2_trees",
                                paste(protein, "_full_time.nwk", sep = ""),
                                sep = "/"))
  
  mya_tree <- read.tree(paste("fig_data", "ace2_trees",
                               paste(protein, "_30my_time.nwk", sep = ""),
                               sep = "/"))
  
  full_rates <- full_df$`Prot Rate`
  names(full_rates) <- full_df$Taxa
  full_ace2_rates <- full_df$`ACE2 Rate`
  names(full_ace2_rates) <- full_df$Taxa
  
  my_rates <- my_df$`Prot Rate`
  names(my_rates) <- my_df$Taxa
  my_ace2_rates <- my_df$`ACE2 Rate`
  names(my_ace2_rates) <- my_df$Taxa
  
  Contrast.full.prot <- pic(full_rates, full_tree)
  Contrast.full.ace2 <- pic(full_ace2_rates, full_tree)
  
  Contrast.my.prot <- pic(my_rates, mya_tree)
  Contrast.my.ace2 <- pic(my_ace2_rates, mya_tree)
  
  full_mod <- lm(Contrast.full.prot ~ Contrast.full.ace2 - 1)
  my_mod <- lm(Contrast.my.prot ~ Contrast.my.ace2 - 1)
  
  full_cor <- cor.table(data.frame(prot=Contrast.full.prot, ace2=Contrast.full.ace2),
                        cor.method = "spearman", cor.type = "contrast")
  my_cor <- cor.table(data.frame(prot=Contrast.my.prot, ace2=Contrast.my.ace2),
                      cor.method = "spearman", cor.type = "contrast")
  
  plot(Contrast.full.prot, Contrast.full.ace2)
  
  return(list(full=summary(full_mod), my=summary(my_mod), 
              full_cor=full_cor, my_cor=my_cor))
}

resid_to_btime_corr <- function(protein) {
  full_df <- read_excel(paste("fig_data", "btime_xys",
                              paste(protein, "_btime_data.xlsx", sep = ""),
                              sep = "/"), sheet = "Full Data")
  mya_df <- read_excel(paste("fig_data", "btime_xys",
                             paste(protein, "_btime_data.xlsx", sep = ""),
                             sep = "/"), sheet = "30MYA")
  
  mod <- lm(`Rate` ~ `Time`, data = full_df)
  mod2 <- lm(`Rate` ~ `Time`, data = mya_df)
  
  full_df <- data.frame(Residual=resid(mod), Time=full_df$`Time`)
  mya_df <- data.frame(Residual=resid(mod2), Time=mya_df$`Time`)
  
  full.test <- cor.test(full_df$Residual, full_df$Time, method="spearman")
  mya.test <- cor.test(mya_df$Residual, mya_df$Time, method="spearman")
  
  full.plot <- ggplot(full_df, aes(x = Time, y = Residual)) + geom_point() +
    ggtitle(paste("Time vs", prot, "Residuals")) + ylab("Rate vs Time Residuals") +
    xlab("Time") + ylim(-0.008, 0.008) +
    geom_smooth(color = "black", size = 1, method = "lm", se = FALSE) +
    theme_half_open() +
    background_grid()
  #mya.plot <- ggplot(mya_df, aes(x = Time, y = Residual)) + geom_point() +
   # ggtitle(paste("Time vs", prot, "Residuals")) + ylab("Rate vs Time Residuals") +
    #xlab("Time Rank") + ylim(-0.008, 0.008) +
    #geom_smooth(color = "black", size = 1, method = "lm", se = FALSE)
  
  return(list(full.rho=full.test$estimate, full.p=full.test$p.value,
              mya.rho=mya.test$estimate, mya.p=mya.test$p.value,
              plot=full.plot))
}

resid_to_resid_XY <- function(protein) {
  full_df <- read_excel(paste("fig_data", "ace2_xys",
                              paste(protein, "_ace2_data.xlsx", sep = ""),
                              sep = "/"), sheet = "Full Data")
  
  mod <- lm(`Prot Rate` ~ `Time`, data = full_df)
  mod2 <- lm(`ACE2 Rate` ~ `Time`, data = full_df)

  full_df <- data.frame(Residual.1=resid(mod), Residual.2=resid(mod2))
  
  plot <- ggplot(full_df, aes(x = Residual.1, y = Residual.2)) + geom_point() +
    ggtitle(paste("ACE2 vs", protein, "Residuals")) + ylab(paste(protein, "vs Time Residuals")) +
    xlab("ACE2 vs Time Residuals") + 
    geom_smooth(color = "black", size = 1, method = "lm", se = FALSE) +
    theme_half_open() +
    background_grid()
  
  return(list(plot=plot))
}

make_p_label <- function(pval) {
    return(paste("P=", signif(pval, 2), sep = ""))
}

compare_extended_branches <- function(protein) {
  full_df <- read_excel(paste("fig_data", "btime_xys",
                              paste(protein, "_btime_data.xlsx", sep = ""),
                              sep = "/"), sheet = "Full Data")
  mya20_df <- read_excel(paste("fig_data", "btime_xys",
                               paste(protein, "_btime_data.xlsx", sep = ""),
                               sep = "/"), sheet = "20MYA")
  mya30_df <- read_excel(paste("fig_data", "btime_xys",
                             paste(protein, "_btime_data.xlsx", sep = ""),
                             sep = "/"), sheet = "30MYA")
  
  full_vs_20 <- inner_join(full_df, mya20_df, by = "Taxa")
  full_vs_20["TimeDiff"] <- full_vs_20["Time.x"] - full_vs_20["Time.y"]
  full_vs_20["RateDiff"] <- full_vs_20["Rate.y"] - full_vs_20["Rate.x"]
  full_vs_20 <- dplyr::filter(full_vs_20, abs(TimeDiff) > 0)
  full_vs_20["Treatment"] <- "0MY vs 20MY"
  
  full_vs_30 <- inner_join(full_df, mya30_df, by = "Taxa")
  full_vs_30["TimeDiff"] <- full_vs_30["Time.x"] - full_vs_30["Time.y"]
  full_vs_30["RateDiff"] <- full_vs_30["Rate.y"] - full_vs_30["Rate.x"]
  full_vs_30 <- dplyr::filter(full_vs_30, abs(TimeDiff) > 0)
  full_vs_30["Treatment"] <- "0MY vs 30MY"
  
  `20_vs_30` <- inner_join(mya20_df, mya30_df, by = "Taxa")
  `20_vs_30`["TimeDiff"] <- `20_vs_30`["Time.x"] - `20_vs_30`["Time.y"]
  `20_vs_30`["RateDiff"] <- `20_vs_30`["Rate.y"] - `20_vs_30`["Rate.x"]
  `20_vs_30` <- dplyr::filter(`20_vs_30`, abs(TimeDiff) > 0)
  `20_vs_30`["Treatment"] <- "20MY vs 30MY"
  
  combined <- rbind(full_vs_20, full_vs_30, `20_vs_30`)
  combined$Treatment <- factor(combined$Treatment, 
                               levels = c("0MY vs 20MY", "20MY vs 30MY", "0MY vs 30MY"))
  
  full_vs_20.test <- wilcox.test(full_vs_20$Rate.x, full_vs_20$Rate.y, paired = TRUE)
  full_vs_30.test <- wilcox.test(full_vs_30$Rate.x, full_vs_30$Rate.y, paired = TRUE)
  `20_vs_30.test` <- wilcox.test(`20_vs_30`$Rate.x, `20_vs_30`$Rate.y, paired = TRUE)
  
  plot <- ggplot(combined, aes(x=Treatment, y=RateDiff)) + geom_boxplot() + geom_hline(yintercept = 0, linetype = "dashed", color="blue") +
    ggtitle(paste("Difference of Extended Branch Rates of", protein)) + 
    ylab("Difference in Rate") + xlab("Branches Compared") + 
    theme(axis.text.x=element_text(size=10)) + stat_n_text(size = 6) +
    stat_summary(geom = 'text', size = 6,
                 label = c(make_p_label(full_vs_20.test$p.value),
                           make_p_label(`20_vs_30.test`$p.value),
                           make_p_label(full_vs_30.test$p.value)), 
                 vjust = -1, fun.y = function(y) {
                   range.y <- range(combined$RateDiff, na.rm = TRUE)
                   pos <- range.y[1] - diff(range.y) * 0.1 # Based on stat_n_text source
                   return(pos)
                 }) +
    theme_half_open(font_size=19, rel_large=0.9) +
    background_grid(major = "y")

  
  return(list(full.20.p=full_vs_20.test$p.value,
              full.30.p=full_vs_30.test$p.value,
              `20.30.p`=`20_vs_30.test`$p.value,
              plot=plot))
}

dist_of_rate_diffs <- function() {
  full_df <- read_excel("fig_Data/rate_bt_corr/btime_rate_data.xlsx", sheet = "Orig Rates")
  mya20_df <- read_excel("fig_Data/rate_bt_corr/btime_rate_data.xlsx", sheet = "20MY Rates")
  mya30_df <-read_excel("fig_Data/rate_bt_corr/btime_rate_data.xlsx", sheet = "30MY Rates")
  
  test_results <- make_empty_df(c("Protein", "MY0.MY20.p", "MY20.MY30.p", "MY0.MY30.p"))
  for (prot in unique(mya30_df$Protein)) {
    full <- full_df[which(full_df$Protein == prot),]
    mya20 <- mya20_df[which(mya20_df$Protein == prot),]
    mya30 <- mya30_df[which(mya30_df$Protein == prot),]
    
    full_vs_20 <- inner_join(full, mya20, by = "Taxon")
    full_vs_20["TimeDiff"] <- full_vs_20["Time.x"] - full_vs_20["Time.y"]
    full_vs_20["RateDiff"] <- full_vs_20["Rate.y"] - full_vs_20["Rate.x"]
    full_vs_20 <- dplyr::filter(full_vs_20, abs(TimeDiff) > 0)
    full_vs_20["Treatment"] <- "0MY vs 20MY"
    
    full_vs_30 <- inner_join(full, mya30, by = "Taxon")
    full_vs_30["TimeDiff"] <- full_vs_30["Time.x"] - full_vs_30["Time.y"]
    full_vs_30["RateDiff"] <- full_vs_30["Rate.y"] - full_vs_30["Rate.x"]
    full_vs_30 <- dplyr::filter(full_vs_30, abs(TimeDiff) > 0)
    full_vs_30["Treatment"] <- "0MY vs 30MY"
    
    `20_vs_30` <- inner_join(mya20, mya30, by = "Taxon")
    `20_vs_30`["TimeDiff"] <- `20_vs_30`["Time.x"] - `20_vs_30`["Time.y"]
    `20_vs_30`["RateDiff"] <- `20_vs_30`["Rate.y"] - `20_vs_30`["Rate.x"]
    `20_vs_30` <- dplyr::filter(`20_vs_30`, abs(TimeDiff) > 0)
    `20_vs_30`["Treatment"] <- "20MY vs 30MY"
    
    full_vs_20.test <- wilcox.test(full_vs_20$Rate.x, full_vs_20$Rate.y, paired = TRUE)$p.value
    full_vs_30.test <- wilcox.test(full_vs_30$Rate.x, full_vs_30$Rate.y, paired = TRUE)$p.value
    `20_vs_30.test` <- wilcox.test(`20_vs_30`$Rate.x, `20_vs_30`$Rate.y, paired = TRUE)$p.value
    
    test_results <- rbind(test_results, data.frame(
      Protein=prot,
      MY0.MY20.p=full_vs_20.test,
      MY20.MY30.p=`20_vs_30.test`,
      MY0.MY30.p=full_vs_30.test,
      MY0.MY20.median.change=median(full_vs_20$RateDiff),
      MY20.MY30.median.change=median(`20_vs_30`$RateDiff),
      MY0.MY30.median.change=median(full_vs_30$RateDiff)
    ))
  }
  
  # FDR Correction
  test_results$MY0.MY20.p.fdr <- p.adjust(test_results$MY0.MY20.p, method = "BH")
  test_results$MY0.MY30.p.fdr <- p.adjust(test_results$MY0.MY30.p, method = "BH")
  test_results$MY20.MY30.p.fdr <- p.adjust(test_results$MY20.MY30.p, method = "BH")
  
  
  breaks <- sapply(seq(0, 1, 0.05), function(x){return(ifelse(x == 1, x, x))})
  plot1.p <- ggplot(test_results, aes(x = MY0.MY20.p)) + geom_histogram(breaks = breaks, fill="#2596be", alpha=0.85) +
    stat_bin(aes(y=..count.., label=..count..), geom="text", vjust=-.5, breaks = breaks) + 
    ggtitle("Paired Wilcoxon Tests P-values of Change in Protein Rate, 0MY vs 20MY") + 
    xlab("Wilcoxon P-value") + geom_vline(xintercept=0.05, col='red') + ylab("Number of Proteins") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    theme_half_open() +
    background_grid(major = "y")
  plot2.p <- ggplot(test_results, aes(x = MY20.MY30.p)) + geom_histogram(breaks = breaks, fill="#2596be", alpha=0.85) +
    stat_bin(aes(y=..count.., label=..count..), geom="text", vjust=-.5, breaks = breaks) + 
    ggtitle("Paired Wilcoxon Tests P-values of Change in Protein Rate, 20MY vs 30MY") + 
    xlab("Wilcoxon P-value") + geom_vline(xintercept=0.05, col='red') + ylab("Number of Proteins") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    theme_half_open() +
    background_grid(major = "y")
  plot3.p <- ggplot(test_results, aes(x = MY0.MY30.p)) + geom_histogram(breaks = breaks, fill="#2596be", alpha=0.85) +
    stat_bin(aes(y=..count.., label=..count..), geom="text", vjust=-.5, breaks = breaks) + 
    ggtitle("Paired Wilcoxon Tests P-values of Change in Protein Rate, 0MY vs 30MY") + 
    xlab("Wilcoxon P-value") + geom_vline(xintercept=0.05, col='red') + ylab("Number of Proteins") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    theme_half_open() +
    background_grid(major = "y")
  
  plot1 <- ggplot(test_results, aes(x = MY0.MY20.p.fdr)) + geom_histogram(breaks = breaks, fill="#2596be", alpha=0.85) +
    stat_bin(aes(y=..count.., label=..count..), geom="text", vjust=-.5, breaks = breaks) + 
    ggtitle("Paired Wilcoxon Tests P-values of Change in Protein Rate, 0MY vs 20MY") + 
    xlab("Wilcoxon P-value (FDR Adjusted)") + geom_vline(xintercept=0.05, col='red') + ylab("Number of Proteins") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    theme_half_open() +
    background_grid(major = "y")
  plot2 <- ggplot(test_results, aes(x = MY20.MY30.p.fdr)) + geom_histogram(breaks = breaks, fill="#2596be", alpha=0.85) +
    stat_bin(aes(y=..count.., label=..count..), geom="text", vjust=-.5, breaks = breaks) + 
    ggtitle("Paired Wilcoxon Tests P-values of Change in Protein Rate, 20MY vs 30MY") + 
    xlab("Wilcoxon P-value (FDR Adjusted)") + geom_vline(xintercept=0.05, col='red') + ylab("Number of Proteins") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    theme_half_open() +
    background_grid(major = "y")
  plot3 <- ggplot(test_results, aes(x = MY0.MY30.p.fdr)) + geom_histogram(breaks = breaks, fill="#2596be", alpha=0.85) +
    stat_bin(aes(y=..count.., label=..count..), geom="text", vjust=-.5, breaks = breaks) + 
    ggtitle("Paired Wilcoxon Tests P-values of Change in Protein Rate, 0MY vs 30MY") + 
    xlab("Wilcoxon P-value (FDR Adjusted)") + geom_vline(xintercept=0.05, col='red') + ylab("Number of Proteins") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    theme_half_open() +
    background_grid(major = "y")
  
  return(list(ORIG=list(MY0.MY20=plot1.p, MY20.MY30=plot2.p, MY0.MY30=plot3.p, test.data=test_results),
    FDR=list(MY0.MY20=plot1, MY20.MY30=plot2, MY0.MY30=plot3, test.data=test_results)))
}

ace2_erc_comparison <- function() {
  df <- read_excel("fig_data/ace2_top40_ercs.xlsx")
  #vennobj <- plotVenn(
  #  list(`Original ERCs`=df$Protein...2[1:20],
  #       `BT-Corrected ERCs`=df$Protein...7[1:20],
  #       `30MY-Adjusted ERCs`=df$Protein...12[1:20])
  #)
  #showSVG(plotVenn(nVennObj = plotVenn(nVennObj = vennobj)), outFile="ace2_top20erc_overlaps.svg", fontScale = 2, borderWidth = 3, labelRegions = FALSE)
  #getVennRegion(vennobj, c("(1) ACE2", "(2) CAT", "(3) CERS3"))
  
  p <- ggvenn(list(`Original ERCs`=df$Protein...2[1:20],
              `BT-Corrected ERCs`=df$Protein...7[1:20],
              `30MY-Adjusted ERCs`=df$Protein...12[1:20]))
  ggsave(paste("figs/", "ACE2_top20_overlaps.png", sep = ""), plot=p, type="cairo-png")
}

make_empty_df <- function(columns) {
  df <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(df) <- columns
  return(df)
}

correlate_resids_to_time <- function() {
  df <- read_excel("fig_Data/rate_bt_corr/btime_rate_data.xlsx", sheet = "Orig Rates")
  prots <- unique(df$Protein)
  results <- make_empty_df(c("Protein", "Rho", "P"))
  for (i in 1:length(prots)) {
    data <- df[which(as.character(df$Protein) == as.character(prots[[i]])),]
    model <- lm(`Rate` ~ `Time`, data = data)
    
    model_df <- data.frame(Residual=resid(model), Time=data$`Time`)
    
    full.test <- cor.test(model_df$Residual, model_df$Time, method="spearman")
    
    results <- rbind(results, data.frame(
      Protein=prots[[i]],
      Rho=full.test$estimate,
      P=full.test$p.value
    ))
  }
  return(results)
}

PROTS = c("ACE2", "GEN1", "XCR1", "CLU", "TMEM63C", "IFNAR2", "APOB", "F5", 
          "PLA2R1", "CAT", "CERS3")


time_regression_df <- make_empty_df(c("Protein", "full.intercept", "full.intercept.p", "full.time", "full.time.p", "full.adj_rsq", "full.p", "full.aov.p",
                                      "mya.intercept", "mya.intercept.p", "mya.time", "mya.time.p", "mya.adj_rsq", "mya.p", "mya.aov.p"))


ace2_regression_df <- make_empty_df(c("Protein", "full.adj_rsq", "full.p", "full.aov.ace2.p", "full.aov.order.p", "full.aov.time.p",
                                      "mya.adj_rsq", "mya.p", "mya.aov.ace2.p", "mya.aov.order.p", "mya.aov.time.p"))

ace2_regression_coefs_df <- make_empty_df(c("Protein", "data", "coefficient.name", "coefficient", "coefficient.p"))


ace2_regression_no_order_df <- make_empty_df(c("Protein", "full.intercept", "full.intercept.p", "full.ace2", "full.ace2.p", "full.time", "full.time.p", "full.adj_rsq", "full.p", "full.aov.ace2.p", "full.aov.time.p",
                                               "mya.intercept", "mya.intercept.p", "mya.ace2", "mya.ace2.p", "mya.time", "mya.time.p", "mya.adj_rsq", "mya.p", "full.aov.ace2.p", "full.aov.time.p"))

ace2_contrasts_df <- make_empty_df(c("Protein", "full.contrast.r", "full.contrast.p",
                                     "mya.contrast.r", "full.contrast.p"))

prot_resid_time_corr_df <- make_empty_df(c("Protein", "full.rho", "full.p", 
                                           "mya.rho", "mya.p"))

prot_branch_extended_df <- make_empty_df(c("Protein", "full.20.p", "full.30.p", "20.30.p"))

for (prot in PROTS[2:length(PROTS)]) {
  p <- plot_ace2_xy(prot)
  ggsave(paste("figs/ACE2_vs_", prot, "_XY_plot.png", sep = ""), plot=p, type="cairo-png")
  
  res <- regression_for_prots_no_order(prot)
  ace2_regression_no_order_df <- rbind(ace2_regression_no_order_df, data.frame(
    Protein = prot,
    full.intercept = res$full.lm$coefficients[[1]],
    full.intercept.p = res$full.lm$coefficients[[10]],
    full.ace2 = res$full.lm$coefficients[[2]],
    full.ace2.p = res$full.lm$coefficients[[11]],
    full.time = res$full.lm$coefficients[[3]],
    full.time.p = res$full.lm$coefficients[[12]],
    full.adj_rsq = res$full.lm$adj.r.squared,
    full.p = pf(res$full.lm$fstatistic[1],
                res$full.lm$fstatistic[2],
                res$full.lm$fstatistic[3],
                lower.tail=FALSE),
    full.aov.ace2.p=res$full.aov[[1]][["Pr(>F)"]][[1]],
    full.aov.time.p=res$full.aov[[1]][["Pr(>F)"]][[2]],
    
    mya.intercept = res$my.lm$coefficients[[1]],
    mya.intercept.p = res$my.lm$coefficients[[10]],
    mya.ace2 = res$my.lm$coefficients[[2]],
    mya.ace2.p = res$my.lm$coefficients[[11]],
    mya.time = res$my.lm$coefficients[[3]],
    mya.time.p = res$my.lm$coefficients[[12]],
    mya.adj_rsq = res$my.lm$adj.r.squared,
    mya.p = pf(res$my.lm$fstatistic[1],
               res$my.lm$fstatistic[2],
               res$my.lm$fstatistic[3],
               lower.tail=FALSE),
    mya.aov.ace2.p=res$my.aov[[1]][["Pr(>F)"]][[1]],
    mya.aov.time.p=res$my.aov[[1]][["Pr(>F)"]][[2]]
  ))
  
  res <- regression_for_prots(prot)
  ace2_regression_df <- rbind(ace2_regression_df, data.frame(
    Protein = prot,
    full.adj_rsq = res$full.lm$adj.r.squared,
    full.p = pf(res$full.lm$fstatistic[1],
                res$full.lm$fstatistic[2],
                res$full.lm$fstatistic[3],
                lower.tail=FALSE),
    full.aov.ace2.p=res$full.aov[[1]][["Pr(>F)"]][[1]],
    full.aov.order.p=res$full.aov[[1]][["Pr(>F)"]][[2]],
    full.aov.time.p=res$full.aov[[1]][["Pr(>F)"]][[3]],
    
    
    mya.adj_rsq = res$my.lm$adj.r.squared,
    mya.p = pf(res$my.lm$fstatistic[1],
               res$my.lm$fstatistic[2],
               res$my.lm$fstatistic[3],
               lower.tail=FALSE),
    mya.aov.ace2.p=res$my.aov[[1]][["Pr(>F)"]][[1]],
    mya.aov.order.p=res$my.aov[[1]][["Pr(>F)"]][[2]],
    mya.aov.time.p=res$my.aov[[1]][["Pr(>F)"]][[3]]
  ))
  
  for (i in 1:nrow(res$full.lm$coefficients)) {
    ace2_regression_coefs_df <- rbind(ace2_regression_coefs_df, data.frame(
      Protein = prot,
      data = "Full",
      coefficient.name = rownames(res$full.lm$coefficients)[[i]],
      coefficient = res$full.lm$coefficients[i, 1],
      coefficient.p = res$full.lm$coefficients[i, 4]
    ))
  }
  
  for (i in 1:nrow(res$my.lm$coefficients)) {
    ace2_regression_coefs_df <- rbind(ace2_regression_coefs_df, data.frame(
      Protein = prot,
      data = "30MY",
      coefficient.name = rownames(res$my.lm$coefficients)[[i]],
      coefficient = res$my.lm$coefficients[i, 1],
      coefficient.p = res$my.lm$coefficients[i, 4]
    ))
  }
  
  res <- contrasts_test(prot)
  ace2_contrasts_df <- rbind(ace2_contrasts_df, data.frame(
    Protein = prot,
    full.contrast.r = res$full_cor$r[[2]],
    full.contrast.r = res$full_cor$P[[2]],
    mya.contrast.r = res$my_cor$r[[2]],
    mya.contrast.r = res$my_cor$P[[2]]
  ))
  
  res <- resid_to_resid_XY(prot)
  ggsave(paste("figs/", prot, "_ACE2_resids_XY_plot.png", sep = ""), plot=res$plot, type="cairo-png")
} 

for (prot in PROTS) {
  p <- plot_time_xy(prot)
  ggsave(paste("figs/", prot, "_vs_btime", "_XY_plot.png", sep = ""), plot=p, type="cairo-png")
  
  res <- time_regression_for_prots(prot)
  time_regression_df <- rbind(time_regression_df, data.frame(
    Protein = prot,
    full.intercept = res$full.lm$coefficients[[1]],
    full.intercept.p = res$full.lm$coefficients[[7]],
    full.time = res$full.lm$coefficients[[2]],
    full.time.p = res$full.lm$coefficients[[8]],
    full.adj_rsq = res$full.lm$adj.r.squared,
    full.p = pf(res$full.lm$fstatistic[1],
                res$full.lm$fstatistic[2],
                res$full.lm$fstatistic[3],
                lower.tail=FALSE),
    full.aov.p=res$full.aov[[1]][["Pr(>F)"]][[1]],
    
    mya.intercept = res$my.lm$coefficients[[1]],
    mya.intercept.p = res$my.lm$coefficients[[7]],
    mya.time = res$my.lm$coefficients[[2]],
    mya.time.p = res$my.lm$coefficients[[8]],
    mya.adj_rsq = res$my.lm$adj.r.squared,
    mya.p = pf(res$my.lm$fstatistic[1],
                res$my.lm$fstatistic[2],
                res$my.lm$fstatistic[3],
                lower.tail=FALSE),
    mya.aov.p=res$my.aov[[1]][["Pr(>F)"]][[1]]
  ))
  
  res <- resid_to_btime_corr(prot)
  ggsave(paste("figs/", prot, "_resid_vs_time_XY_plot.png", sep = ""), plot=res$plot, type="cairo-png")
  prot_resid_time_corr_df <- rbind(prot_resid_time_corr_df, data.frame(
    Protein = prot,
    full.rho = res$full.rho,
    full.p = res$full.p,
    mya.rho = res$mya.rho,
    mya.p = res$mya.p
  ))
  
  res <- compare_extended_branches(prot)
  ggsave(paste("figs/", prot, "_btime_extended.png", sep = ""), plot=res$plot, type="cairo-png")
  prot_branch_extended_df <- rbind(prot_branch_extended_df, data.frame(
    Protein = prot,
    full.20.p = res$full.20.p,
    full.30.p = res$full.30.p,
    `20.30.p` = res$`20.30.p`
  ))
} 

res <- dist_of_rate_diffs()
ggsave("figs/pval_dist_rate_diff_MY0_MY20.png", plot=res$MY0.MY20, type="cairo-png")
ggsave("figs/pval_dist_rate_diff_MY20_MY30.png", plot=res$MY20.MY30, type="cairo-png")
ggsave("figs/pval_dist_rate_diff_MY0_MY30.png", plot=res$MY0.MY30, type="cairo-png")
write_xlsx(res$test.data, "figs/branch_extension_wilcoxon_results.xlsx")


write_xlsx(correlate_resids_to_time(), "figs/time_resids_to_time_corr.xlsx")

write_xlsx(time_regression_df, "figs/time_regression.xlsx")
write_xlsx(ace2_regression_df, "figs/ace2_regression.xlsx")
write_xlsx(ace2_regression_coefs_df, "figs/ace2_regression_coef.xlsx")
write_xlsx(ace2_regression_no_order_df, "figs/ace2_regression_no_order.xlsx")
write_xlsx(ace2_contrasts_df, "figs/ace2_contrasts.xlsx")
write_xlsx(prot_resid_time_corr_df, "figs/prot_resid_time_corr.xlsx")
write_xlsx(prot_branch_extended_df, "figs/prot_branch_extension.xlsx")

# Final pub figs
## BT-Rate plots
ace2_res <- plot_time_xy("ACE2", TRUE)
gen1_res <- plot_time_xy("GEN1", TRUE)
xcr1_res <- plot_time_xy("XCR1", TRUE)
clu_res <- plot_time_xy("CLU", TRUE)
apob_res <- plot_time_xy("APOB", TRUE)
ifnar2_res <- plot_time_xy("IFNAR2", TRUE)
original_plot <- (ace2_res$full.plot | gen1_res$full.plot) / (xcr1_res$full.plot | clu_res$full.plot) / (apob_res$full.plot | ifnar2_res$full.plot)
compare_plot <- (ace2_res$full.plot | ace2_res$my.plot) / (gen1_res$full.plot | gen1_res$my.plot) / (xcr1_res$full.plot | xcr1_res$my.plot) / (clu_res$full.plot | clu_res$my.plot)
ggsave("final files/S3_Orig_Rate_BT.png", plot=original_plot, type="cairo-png", dpi = 300, width = 20, height = 15, units = "in")
ggsave("final files/S4_Orig_vs_30MY_Rate_BT.png", plot=compare_plot, type="cairo-png", dpi = 300, width = 20, height = 20, units = "in")

## Rate extension boxplots
ace2_res <- compare_extended_branches("ACE2")
gen1_res <- compare_extended_branches("GEN1")
xcr1_res <- compare_extended_branches("XCR1")
clu_res <- compare_extended_branches("CLU")
apob_res <- compare_extended_branches("APOB")
ifnar2_res <- compare_extended_branches("IFNAR2")

plot <- (ace2_res$plot | gen1_res$plot) / (xcr1_res$plot | clu_res$plot) / (apob_res$plot | ifnar2_res$plot)
ggsave("final files/S6_rate_diff_boxplots.png", plot=plot, type="cairo-png", dpi=300, width = 20, height = 20, units = "in")

## P-value distributions
res <- dist_of_rate_diffs()
plot <- res$ORIG$MY0.MY20 / res$ORIG$MY20.MY30 / res$ORIG$MY0.MY30
ggsave("final files/S7_wilcox_distribution.png", plot=plot, type="cairo-png", dpi=300, width = 10, height = 15, units = "in")

