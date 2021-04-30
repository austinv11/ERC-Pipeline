library(ggplot2)
library(readxl)
library(writexl)
# library(nVennR)
library(ggvenn)
library(dplyr)
library(gridExtra)
library(ape)
library(picante)
library(Cairo)
library(EnvStats)

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
    ggtitle(paste("Original Branch Data (ACE2 vs ", protein, ")", sep = "")) + ylab(paste(protein, "Rate Rank")) + geom_smooth(color = "black", size = 1, method = "lm", se = FALSE)
  plot2 <- ggplot(mya_df, aes(x = `ACE2 Rate Rank`, y = `Prot Rate Rank`, color=order))  + 
    geom_point()+ scale_color_manual("Order", values = c("purple", "blue", "#00BA38", "black", "red", "#E58700")) +
    ggtitle(paste("30MY Branch Data (ACE2 vs ", protein, ")", sep = "")) + ylab(paste(protein, "Rate Rank")) + geom_smooth(color = "black", size = 1, method = "lm", se = FALSE)

  return(grid.arrange(plot1, plot2, ncol=2))
}

plot_time_xy <- function(protein) {
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
  
  plot1 <- ggplot(full_df, aes(x = `Time`, y = `Rate`, color=order)) + 
    geom_point()+ scale_color_manual("Order", values = c("purple", "blue", "#00BA38", "black", "red", "#E58700")) +
    ggtitle(paste("Original Branch Data (Time vs ", protein, ")", sep = "")) + ylab(paste(protein, "Rate")) +
    xlim(0, 200) +
    labs(caption=paste("rho=", stats_df[stats_df$Treatment == 'Full Taxa',c(2, 3)][[1]],
                       " p=", stats_df[stats_df$Treatment == 'Full Taxa',c(2, 3)][[2]], sep = "")) + geom_smooth(color = "black", size = 0.75, method = "lm", se = FALSE, linetype="dashed")
  plot2 <- ggplot(mya_df, aes(x = `Time`, y = `Rate`, color=order))  + 
    geom_point()+ scale_color_manual("Order", values = c("purple", "blue", "#00BA38", "black", "red", "#E58700")) +
    ggtitle(paste("30MY Branch Data (Time vs ", protein, ")", sep = "")) + ylab(paste(protein, "Rate")) +
    xlim(0, 200) + 
    labs(caption=paste("rho=", stats_df[stats_df$Treatment == '30MY Taxa',c(2, 3)][[1]],
                       " p=", stats_df[stats_df$Treatment == '30MY Taxa',c(2, 3)][[2]], sep = "")) + geom_smooth(color = "black", size = 0.75, method = "lm", se = FALSE, linetype="dashed")
  
  ggsave(paste("figs/", protein, "_vs_time_0MY.png", sep = ""), plot=plot1, type="cairo-png")
  ggsave(paste("figs/", protein, "_vs_time_30MY.png", sep = ""), plot=plot2, type="cairo-png")
  
  return(grid.arrange(plot1, plot2, ncol=2))
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
    geom_smooth(color = "black", size = 1, method = "lm", se = FALSE)
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
    geom_smooth(color = "black", size = 1, method = "lm", se = FALSE)
  
  return(list(plot=plot))
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
  
  plot <- ggplot(combined, aes(x=Treatment, y=RateDiff)) + geom_boxplot() + 
    ggtitle(paste("Difference of Extended Branch Rates of", protein)) + 
    ylab("Difference in Rate") + xlab("Branches Compared") + 
    theme(axis.text.x=element_text(size=10)) + stat_n_text()
  
  full_vs_20.test <- wilcox.test(full_vs_20$Rate.x, full_vs_20$Rate.y, paired = TRUE)
  full_vs_30.test <- wilcox.test(full_vs_30$Rate.x, full_vs_30$Rate.y, paired = TRUE)
  `20_vs_30.test` <- wilcox.test(`20_vs_30`$Rate.x, `20_vs_30`$Rate.y, paired = TRUE)
  ggsave(paste("figs/", protein, "_btime_extended.png", sep = ""), plot=plot, type="cairo-png")
  
  
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
      MY0.MY30.p=full_vs_30.test
    ))
  }
  
  # FDR Correction
  test_results$MY0.MY20.p.fdr <- p.adjust(test_results$MY0.MY20.p, method = "BH")
  test_results$MY0.MY30.p.fdr <- p.adjust(test_results$MY0.MY30.p, method = "BH")
  test_results$MY20.MY30.p.fdr <- p.adjust(test_results$MY20.MY30.p, method = "BH")
  
  plot1 <- ggplot(test_results, aes(x = MY0.MY20.p.fdr)) + geom_histogram(binwidth = 0.05) +
    stat_bin(aes(y=..count.., label=..count..), geom="text", vjust=-.5, binwidth = 0.05) + 
    ggtitle("Paired Wilcoxon Tests p-values of change in protein rate, 0MY vs 20MY") + 
    xlab("Wilcoxon p-value (FDR Adjusted)") + geom_vline(xintercept=0.025, col='red') + ylab("Number of Proteins")
  plot2 <- ggplot(test_results, aes(x = MY20.MY30.p.fdr)) + geom_histogram(binwidth = 0.05) +
    stat_bin(aes(y=..count.., label=..count..), geom="text", vjust=-.5, binwidth = 0.05) + 
    ggtitle("Paired Wilcoxon Tests p-values of change in protein rate, 20MY vs 30MY") + 
    xlab("Wilcoxon p-value (FDR Adjusted)") + geom_vline(xintercept=0.025, col='red') + ylab("Number of Proteins")
  plot3 <- ggplot(test_results, aes(x = MY0.MY30.p.fdr)) + geom_histogram(binwidth = 0.05) +
    stat_bin(aes(y=..count.., label=..count..), geom="text", vjust=-.5, binwidth = 0.05) + 
    ggtitle("Paired Wilcoxon Tests p-values of change in protein rate, 0MY vs 30MY") + 
    xlab("Wilcoxon p-value (FDR Adjusted)") + geom_vline(xintercept=0.025, col='red') + ylab("Number of Proteins")
  
  return(list(MY0.MY20=plot1, MY20.MY30=plot2, MY0.MY30=plot3, test.data=test_results))
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

PROTS = c("ACE2", "GEN1", "XCR1", "CLU", "IFNAR2", "APOB", "PLA2R1", "CAT", 
          "CERS3")


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

write_xlsx(time_regression_df, "figs/time_regression.xlsx")
write_xlsx(ace2_regression_df, "figs/ace2_regression.xlsx")
write_xlsx(ace2_regression_coefs_df, "figs/ace2_regression_coef.xlsx")
write_xlsx(ace2_regression_no_order_df, "figs/ace2_regression_no_order.xlsx")
write_xlsx(ace2_contrasts_df, "figs/ace2_contrasts.xlsx")
write_xlsx(prot_resid_time_corr_df, "figs/prot_resid_time_corr.xlsx")
write_xlsx(prot_branch_extended_df, "figs/prot_branch_extension.xlsx")

