library(ggplot2)
library(stringr)
library(reshape2)
library(plyr)

results <- read.csv("../../experiments/results.csv")
configuration <- read.csv("../../experiments/all_configurations.csv")

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

results.with.config$global.effect.err <- results.with.config$global.effect.est - results.with.config$global.effect.act
single.inst.per.trial <- subset(results.with.config, setting == 0 & graph.cluster.randomization == TRUE & 
                                  ((of.beta == 1 & ot.beta == 1) | (of.beta == 2 & ot.beta == 2 & exposure.type %in% c("exponential", "sigmoid"))))
experimental.results <- ddply(single.inst.per.trial, .(method, exposure.type, graph.type), summarize, mean.err = mean(global.effect.err), sd.err=sd(global.effect.err))
experimental.results$exposure.type <- revalue(experimental.results$exposure.type, c("exponential"="Exponential", "linear"="Linear", "rbf-friends"="RBF", "sigmoid"="Sigmoid"))
experimental.results$graph.type <- revalue(experimental.results$graph.type, c("barabasi-albert"="Pref.\nAttach.", "small-world"="Small\nWorld"))
g <- ggplot(subset(experimental.results, method != "Actual" & method != "Obs-RKS-Sufficient" & method != "Exp-LM-IND"), 
       aes(x=method, y=mean.err, ymin=mean.err-sd.err, ymax=mean.err+sd.err, color=method)) + 
    facet_grid(graph.type~exposure.type) + geom_errorbar(size=2) + geom_hline(yintercept=0, color="blue") + theme_bw(base_size=20) + 
    labs(x="Method", y="Estimated - Actual") + scale_x_discrete(breaks=c("Exp-HT", "Exp-LM-INT", "Obs-GBM-Sufficient"), labels=c("HT", "LM", "GBM")) + 
    scale_color_manual(values=c("Exp-HT"="slateblue", "Exp-LM-INT"="gray47", "Obs-GBM-Sufficient"="red")) + guides(color="none")
png("experimental-perf.png", width=800, height=300)
print(g)
dev.off()
print(g)

results.summ <- ddply(subset(results.with.config, !is.na(value.est)), .(config, method, setting, exposure.type, effect.type), 
                      summarize, rmse = sqrt(sum((value.est - value.act)^2) / length(value.act)))
# bring the configuration back
results.summ <- merge(results.summ, configuration, by.x="config", by.y="X")

table.by.method <- ddply(subset(results.summ, graph.cluster.randomization == FALSE & 
                                  effect.type == "Marginal Peer" & method != "Actual" & 
                                  treatment.autocorr.coeff == 2), 
                            .(exposure.type.x, method), summarize, meanrmse = mean(rmse))
results.table <- reshape(table.by.method, timevar ="exposure.type.x", direction="wide", idvar="method")
library(xtable)
print(xtable(results.table))

ggplot(subset(results.summ, graph.cluster.randomization == FALSE & effect.type == "Marginal Individual" &
                method != "Actual" & confounding.coeff == 1 & treatment.autocorr.coeff == 2 & 
                graph.type == "small-world" & method %in% c("Exp-LM-IND", "Obs-LM-Simple", "Obs-GBM-Sufficient", "Obs-LM-Sufficient")), 
       aes(x=factor(setting), y=rmse, fill=method)) + facet_grid(exposure.type.x~., scales="free") +
  geom_violin()

ggplot(subset(results.summ, graph.cluster.randomization == FALSE & effect.type == "Marginal Peer" &
                method != "Actual" & confounding.coeff == 1 & treatment.autocorr.coeff == 10 & of.beta > 2), 
       aes(x=factor(setting), y=rmse, fill=method)) + facet_grid(exposure.type.x~graph.type, scales="free") +
      geom_boxplot()

# plot some outcome functions
plotFriends <- function(outcome.fun, file.name=NULL) {
  fp <- seq(0, 1, length.out=2000)
  outcome <- aaply(fp, 1, function(v) outcome.fun(friendt=v))
  g <- qplot(fp, outcome, geom="smooth", se=FALSE) + theme_bw(base_size=20) + labs(x=expression(theta), y="Y")  
  if(!is.null(file.name)) {
    png(file.name, width=250, height=250)
    print(g)
    dev.off()
  } else {
    print(g)
  }
}
gd <- generate.by.index(configuration, 34)
plotFriends(gd$outcome.function, "sigmoid.png")
gd <- generate.by.index(configuration, 2256)
plotFriends(gd$outcome.function, "rbf.png")
gd <- generate.by.index(configuration, 2253)
plotFriends(gd$outcome.function, "linear.png")
gd <- generate.by.index(configuration, 2243)
plotFriends(gd$outcome.function, "exponential.png")

