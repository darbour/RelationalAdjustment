library(ggplot2)
library(stringr)
library(reshape2)
library(plyr)

results <- read.csv("../../experiments/results_rw.csv")
configuration <- read.csv("../../experiments/all_configurations_rw.csv")

results$global.effect <- results$global_treatment - results$global_control
results$indiv.effect <- results$mt_1 - results$mt_0

# reshape the data into a plottable form
melt.results <- melt(results, idvar=c("method", "config", "trial", "global.effect", "indiv.effect"), 
                     measure.vars=c(paste0("mf_", seq(0, 1, length.out=11))))
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

results.with.config <- subset(results.with.config, 
  (method == "Exp-LM-IND" & graph.cluster.randomization == TRUE) |
    (method %in% c("Obs-GBM-Sufficient", "Exp-GBM") & graph.cluster.randomization == FALSE))

results.with.config$global.effect.err <- results.with.config$global.effect.est - results.with.config$global.effect.act

results.with.config$exposure.type <- revalue(results.with.config$exposure.type, c("linear"="Linear", "rbf-friends"="RBF", "sigmoid"="Sigmoid"))
results.with.config$method <- revalue(results.with.config$method, 
                                        c("Obs-GBM-Sufficient"="ObsGBM", "Exp-LM-IND"="ExpLM", "Exp-GBM"="ObsGBM-U"))
g <- ggplot(subset(results.with.config, method %in% c("ObsGBM", "ExpLM") & (treatment.autocorr.coeff == 2 | graph.cluster.randomization)), 
       aes(x=method, y=global.effect.err, fill=method)) + geom_boxplot(notch=TRUE) + theme_bw() + guides(fill="none") + 
      facet_wrap(~exposure.type) + theme_bw(base_size=15) +
      labs(x="Method", y="Error in Total Effect") +
      scale_fill_manual(values=c("ExpHT"="slateblue", "ExpLM"="gray47", "ObsGBM"="firebrick", "ObsLM"="deepskyblue3", "ObsGBM-U"="indianred2"))
print(g)
png("enron-total-effect.png", width=400, height=300)
print(g)
dev.off()

single.inst.per.trial <- subset(results.with.config, setting == 0)
single.inst.per.trial$indiv.effect.err <- sqrt((single.inst.per.trial$indiv.effect.est - single.inst.per.trial$indiv.effect.act) ** 2)
indiv.effect.perf <- ddply(subset(single.inst.per.trial, (treatment.autocorr.coeff == 2 | graph.cluster.randomization) & method %in% c("ExpLM", "ObsGBM")), 
      .(method, exposure.type), summarize, mean.err = mean(indiv.effect.err), sd.err = sd(indiv.effect.err))
indiv.effect.perf$cell <- paste0(round(indiv.effect.perf$mean.err, 4), " (", round(indiv.effect.perf$sd.err, 4), ")")
indiv.effect.perf <- subset(indiv.effect.perf, select=-c(mean.err, sd.err))
results.table <- reshape(indiv.effect.perf, timevar ="method", direction="wide", idvar="exposure.type")
colnames(results.table) <- c("", "ExpLM", "ObsGBM")
library(xtable)
cat("Individual effect performance:\n")
print(xtable(results.table), row.names=FALSE)


experimental.results <- ddply(single.inst.per.trial, .(method, exposure.type, graph.type), summarize, mean.err = mean(global.effect.err), sd.err=sd(global.effect.err))
#experimental.results$exposure.type <- revalue(experimental.results$exposure.type, c("exponential"="Exponential", "linear"="Linear", "rbf-friends"="RBF", "sigmoid"="Sigmoid"))
#experimental.results$graph.type <- revalue(experimental.results$graph.type, c("barabasi-albert"="Pref.\nAttach.", "small-world"="Small\nWorld"))
g <- ggplot(subset(experimental.results, method != "Actual" & method != "Obs-RKS-Sufficient" & method != "Exp-LM-INT"), 
       aes(x=method, y=mean.err, ymin=mean.err-sd.err, ymax=mean.err+sd.err, color=method)) + 
    facet_grid(graph.type~exposure.type) + geom_errorbar(size=2) + geom_hline(yintercept=0, color="blue") + theme_bw(base_size=20) + 
    labs(x="Method", y="Estimated - Actual") + scale_x_discrete(breaks=c("Exp-HT", "Exp-LM-IND", "Obs-GBM-Sufficient"), labels=c("HT", "LM", "GBM")) + 
    scale_color_manual(values=c("Exp-HT"="slateblue", "Exp-LM-IND"="gray47", "Obs-GBM-Sufficient"="red")) + guides(color="none")
png("experimental-perf.png", width=800, height=300)
print(g)
dev.off()
print(g)

# prune exponential results that exploded
results.summ <- ddply(subset(results.with.config, !is.na(value.est)), .(config, method, setting, exposure.type, effect.type), 
                      summarize, rmse = sqrt(sum((value.est - value.act)^2) / length(value.act)))
# bring the configuration back
results.summ <- merge(results.summ, configuration, by.x="config", by.y="X")
  
table.by.method <- ddply(subset(results.summ, (treatment.autocorr.coeff == 2 | graph.cluster.randomization)),
                              .(exposure.type.x, method, of.beta, ot.beta), 
                            summarize, meanrmse = mean(rmse), sdrmse=sd(rmse))
table.by.method <- subset(subset(table.by.method, method %in% c("ExpLM", "ObsGBM")), select=-c(of.beta, ot.beta))
table.by.method$contents <- with(table.by.method, paste0(round(meanrmse, ifelse(meanrmse > 100, 0, 4)), " (", round(sdrmse, ifelse(meanrmse > 100, 0, 3)), ")"))
table.by.method <- subset(table.by.method, select=-c(meanrmse, sdrmse))
results.table <- reshape(table.by.method, timevar ="method", direction="wide", idvar="exposure.type.x")
colnames(results.table) <- str_replace(colnames(results.table), "contents.", "")
library(xtable)
cat("Peer effect performance:\n")
print(xtable(results.table), row.names=FALSE)

ggplot(subset(results.summ, graph.cluster.randomization == FALSE & effect.type == "Marginal Individual" &
                method != "Actual" & confounding.coeff == 1 & treatment.autocorr.coeff == 2 & 
                graph.type == "small-world" & method %in% c("Exp-LM-IND", "Obs-LM-Simple", "Obs-GBM-Sufficient", "Obs-LM-Sufficient")), 
       aes(x=factor(setting), y=rmse, fill=method)) + facet_grid(exposure.type.x~., scales="free") +
  geom_violin()

ggplot(subset(results.summ, graph.cluster.randomization == FALSE & effect.type == "Marginal Peer" &
                method != "Actual" & confounding.coeff == 1 & treatment.autocorr.coeff == 10 & of.beta > 2), 
       aes(x=factor(setting), y=rmse, fill=method)) + facet_grid(exposure.type.x~graph.type, scales="free") +
      geom_boxplot()
