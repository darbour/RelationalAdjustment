library(ggplot2)
library(stringr)
library(reshape2)

results <- read.csv("../../experiments/results.csv")
configuration <- read.csv("../../experiments/all_configurations.csv")

# reshape the data into a plottable form
melt.results <- melt(results, idvar=c("method", "config", "trial"), 
                     measure.vars=c(paste0("mf_", seq(0, 1, length.out=11)), 
                                    "global_treatment", "global_control", "mt_0", "mt_1"))
melt.results$effect.type <- ifelse(str_detect(melt.results$variable, "mf_"), "Marginal Peer", "")
melt.results$effect.type <- ifelse(str_detect(melt.results$variable, "mt_"), "Marginal Individual", melt.results$effect.type)
melt.results$effect.type <- ifelse(str_detect(melt.results$variable, "global_"), "Global", melt.results$effect.type)

melt.results[melt.results$effect.type == "Marginal Peer", "setting"] <- as.numeric(str_replace(melt.results[melt.results$effect.type == "Marginal Peer", "variable"], "mf_", ""))
melt.results[melt.results$effect.type == "Marginal Individual", "setting"] <- as.numeric(str_replace(melt.results[melt.results$effect.type == "Marginal Individual", "variable"], "mt_", ""))
melt.results[melt.results$effect.type == "Global" & melt.results$variable == "global_treatment", "setting"] <- 1
melt.results[melt.results$effect.type == "Global" & melt.results$variable == "global_control", "setting"] <- 0

results.with.actuals <- merge(melt.results, subset(subset(melt.results, method=="Actual"), select=-method), 
                              by=c("config", "trial", "variable", "effect.type", "setting"), suffixes=c(".est", ".act"))
results.with.config <- merge(results.with.actuals, configuration, by.x="config", by.y="X")

results.summ <- ddply(subset(results.with.config, !is.na(value.est)), .(config, method, setting, exposure.type, effect.type), 
                      summarize, rmse = sqrt(sum((value.est - value.act)^2) / length(value.act)))
# bring the configuration back
results.summ <- merge(results.summ, configuration, by.x="config", by.y="X")

ggplot(subset(results.summ, graph.cluster.randomization == FALSE & effect.type == "Marginal Peer" &
                method != "Actual" & confounding.coeff == 1 & treatment.autocorr.coeff == 0), 
       aes(x=factor(setting), y=rmse, fill=method)) + facet_grid(exposure.type.x~graph.type, scales="free") +
      geom_boxplot()
