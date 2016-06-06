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

specific.results <- subset(results.with.config, of.beta == 10 & ot.beta == 10 &
                           ((graph.cluster.randomization == FALSE & treatment.autocorr.coeff == 2 & confounding.coeff == 1 &
                              method %in% c("Obs-GBM-Sufficient", "Actual")) | (method  == "Exp-LM-IND" & graph.cluster.randomization == TRUE)) &
                             exposure.type == "sigmoid" & graph.type == "small-world" & p == 0.1)
g <- ggplot(specific.results, aes(x=setting, y=value.est, color=method)) + 
  geom_smooth(size=1) + theme_bw(base_size=20) + labs(x=expression(theta), y="Outcome") + 
  scale_color_manual(breaks=c("Actual", "Obs-GBM-Sufficient", "Exp-LM-IND"), 
                       labels=c("Actual", "ObsGBM", "ExpLM"), values=c(1, 3, 2), name="Method") +
  theme(legend.position="top") + guides(color=guide_legend(ncol=2))
print(g)
png("sigmoid-example.png", width=400, height=300)
print(g)
dev.off()


observational.configs <- subset(results.with.config, 
                                setting == 0 & graph.cluster.randomization == FALSE & 
                                  ((of.beta == 1 & ot.beta == 1) | (of.beta == 2 & ot.beta == 2 & exposure.type %in% c("sigmoid"))) & 
                                  method %in% c("Obs-GBM-Sufficient", "Obs-LM-Sufficient", "Exp-GBM"))
observational.configs$method <- revalue(observational.configs$method, 
                                        c("Obs-GBM-Sufficient"="ObsGBM", "Obs-LM-Sufficient"="ObsLM", 
                                          "Exp-GBM"="ObsGBM-U"))
experimental.configs <- subset(results.with.config, 
                               setting == 0 & graph.cluster.randomization == TRUE &
                                 ((of.beta == 1 & ot.beta == 1) | (of.beta == 2 & ot.beta == 2 & exposure.type %in% c("sigmoid"))) & 
                                 method %in% c("Exp-GBM", "Exp-HT", "Exp-LM-IND"))
experimental.configs$method <- revalue(experimental.configs$method, 
                                        c("Exp-GBM"="ExpGBM", "Exp-HT"="ExpHT", 
                                          "Exp-LM-IND"="ExpLM"))
compare.exp.obs <- rbind(observational.configs, experimental.configs)
g <- ggplot(subset(compare.exp.obs, exposure.type == "linear" & method != "ObsLM"), 
       aes(x=method, y=global.effect.err, fill=method)) + geom_boxplot(notch=TRUE) + theme_bw() + guides(fill="none") + 
      labs(x="Method", y="Error in Total Effect") + geom_vline(xintercept=2.5, size=1, linetype="longdash", color="gray47") +
      scale_fill_manual(values=c("ExpHT"="slateblue", "ExpLM"="gray47", "ObsGBM"="firebrick", "ObsLM"="deepskyblue3", "ObsGBM-U"="indianred2"))
print(g)
png("exp-obs-compare.png", width=400, height=300)
print(g)
dev.off()

experimental.configs2 <- merge(subset(experimental.configs, select=-confounding.coeff), data.frame(confounding.coeff=c(1, 3)))
compare.exp.obs <- rbind(observational.configs, experimental.configs2)
plot.compare.dat <- subset(compare.exp.obs, exposure.type != "exponential" & method != "ObsLM" & 
                             (treatment.autocorr.coeff == 2 | 
                                (treatment.autocorr.coeff == 0 & graph.cluster.randomization == TRUE)))
plot.compare.dat$exposure.type <- revalue(plot.compare.dat$exposure.type, c("linear"="Linear", "rbf-friends"="RBF", "sigmoid"="Sigmoid"))
plot.compare.dat$confounding.coeff <- paste0("Confounding ", plot.compare.dat$confounding.coeff)
g <- ggplot(subset(plot.compare.dat, confounding.coeff == "Confounding 1" & method != "ExpGBM" & method != "ObsGBM-U"), aes(x=method, y=global.effect.err, fill=method)) + 
    geom_violin(scale="width") + facet_grid(~exposure.type, scales="free") + 
    theme_bw(base_size=16) + guides(fill="none") + labs(x="Method", y="Error in Total Effect") + 
  geom_vline(xintercept=1.5, linetype="longdash", color="gray47")
print(g)
png("cross-exposure-compare.png", width=1000, height=300)
print(g)
dev.off()

single.inst.per.trial <- subset(results.with.config, setting == 0 & method != "Actual" & 
                                  ((graph.cluster.randomization == TRUE & method %in% c("Exp-LM-IND")) | 
                                     (graph.cluster.randomization == FALSE & method %in% c("Obs-GBM-Sufficient") & 
                                        treatment.autocorr.coeff == 2 & confounding.coeff == 1)
                                  ))
single.inst.per.trial$indiv.effect.err <- sqrt((single.inst.per.trial$indiv.effect.est - single.inst.per.trial$indiv.effect.act) ** 2)
indiv.effect.perf <- ddply(subset(single.inst.per.trial, ((of.beta == 5 & ot.beta == 5) | (of.beta == 2 & ot.beta == 2 & exposure.type == "sigmoid")) & exposure.type != "exponential"), 
      .(method, exposure.type), summarize, mean.err = mean(indiv.effect.err), sd.err = sd(indiv.effect.err))
indiv.effect.perf$cell <- paste0(round(indiv.effect.perf$mean.err, 4), " (", round(indiv.effect.perf$sd.err, 4), ")")
indiv.effect.perf <- subset(indiv.effect.perf, select=-c(mean.err, sd.err))
results.table <- reshape(indiv.effect.perf, timevar ="method", direction="wide", idvar="exposure.type")
colnames(results.table) <- c("", "Exp. LM", "Obs. GLM")
library(xtable)
cat("--------------\n")
cat("RMSE individual effect estimation\n")
cat("--------------\n")
print(xtable(results.table), row.names=FALSE)


experimental.results <- ddply(single.inst.per.trial, .(method, exposure.type, graph.type), summarize, mean.err = mean(global.effect.err), sd.err=sd(global.effect.err))
experimental.results$exposure.type <- revalue(experimental.results$exposure.type, c("exponential"="Exponential", "linear"="Linear", "rbf-friends"="RBF", "sigmoid"="Sigmoid"))
experimental.results$graph.type <- revalue(experimental.results$graph.type, c("barabasi-albert"="Pref.\nAttach.", "small-world"="Small\nWorld"))
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
results.prune.exponential <- subset(results.with.config, confounding.coeff < 3 & of.beta < 10 & ot.beta < 10)
results.summ <- ddply(subset(results.prune.exponential, !is.na(value.est)), .(config, method, setting, exposure.type, effect.type), 
                      summarize, rmse = sqrt(sum((value.est - value.act)^2) / length(value.act)))
# bring the configuration back
results.summ <- merge(results.summ, configuration, by.x="config", by.y="X")
  
table.by.method <- ddply(subset(results.summ, effect.type == "Marginal Peer" & method != "Actual" & 
                                ((graph.cluster.randomization == TRUE & method %in% c("Exp-LM-IND")) | 
                                   (graph.cluster.randomization == FALSE & method %in% c("Obs-GBM-Sufficient") & 
                                      treatment.autocorr.coeff == 2 & confounding.coeff == 1)
                                )), .(exposure.type.x, method, of.beta, ot.beta), 
                            summarize, meanrmse = mean(rmse), sdrmse=sd(rmse))
table.by.method <- subset(subset(table.by.method, of.beta >= 2 & ot.beta >= 2), select=-c(of.beta, ot.beta))
table.by.method$contents <- with(table.by.method, paste0(round(meanrmse, ifelse(meanrmse > 100, 0, 4)), " (", round(sdrmse, ifelse(meanrmse > 100, 0, 3)), ")"))
table.by.method$exposure.type.x <- revalue(table.by.method$exposure.type.x, c("exponential"="Exponential", "linear"="Linear", "rbf-friends"="RBF", "sigmoid"="Sigmoid"))
table.by.method$method <- revalue(table.by.method$method, c("Exp-LM-IND"="Exp. LM", "Obs-GBM-Sufficient"="Obs. GBM", "Exp-GBM"="Obs. GBM (unadj)"))
table.by.method <- subset(table.by.method, select=-c(meanrmse, sdrmse))
results.table <- reshape(table.by.method, timevar ="method", direction="wide", idvar="exposure.type.x")
colnames(results.table) <- str_replace(colnames(results.table), "contents.", "")
library(xtable)
cat("-----------\n")
cat("RMSE marginal peer effects\n")
cat("-----------\n")
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

