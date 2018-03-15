library(mvtnorm)
library(doParallel)
library(mcclust.mod)
library(BNPmix)

################################################################################################
# SIMULATE THE DATA
################################################################################################

grid <- matrix(0, ncol = 2, nrow = 2)

ssize <- c(20, 40, 100, 300, 1000)
k <- 1
dataetc <- list()

set.seed(42)
for(i in 1:10){
  for(j in 1:5){
    dataetc[[k]] <- as.matrix(rbind(rmvnorm(ssize[j]/4, mean = c(-2,-2), sigma = diag(1,2)),
                                    rmvnorm(ssize[j]/4, mean = c(2,-2), sigma = diag(1,2)),
                                    rmvnorm(ssize[j]/4, mean = c(-2,2), sigma = diag(1,2)),
                                    rmvnorm(ssize[j]/4, mean = c(2,2), sigma = diag(1,2))))
    k <- k+1 
  }
}

################################################################################################
# ESTIMATE MODELS & PARTITION & CREDIB. BALLS
################################################################################################

detectCores()
no_cores <- 5
registerDoParallel(cores = no_cores)
cl <- makeCluster(no_cores)
results <- foreach(i = 1:50, .packages = c("BNPmixturesVar", "mcclust.ext")) %dopar% {
  mod <-  DPmixMulti(nsim = 25000,
                     nburn = 5000,
                     napprox = 100,
                     nparam = round(2 * nrow(dataetc[[i]])),
                     d = 2,
                     grid_l = nrow(grid),
                     data = as.matrix(dataetc[[i]]),
                     grid = as.matrix(grid),
                     mu_start = c(0,3),
                     Lambda_start = diag(1,2),
                     theta_start = 1,
                     m0 = c(0,3),
                     B0 = diag(1,2),
                     nu0 = 9,
                     sigma = diag(1,2),
                     b1 = 4,
                     B1 = diag(1,2),
                     m1 = c(0,0),
                     M1 = diag(1,2),
                     s1 = 6,
                     S1 = diag(1,2),
                     nupd = 200,
                     plim = 10,
                     fix = TRUE)
  psm <- comp.psm(mod[[2]] + 1)
  VI <- minVI(psm, mod[[2]] + 1, method = "all", include.greedy=TRUE)
  BI <- minbinder.ext(psm, mod[[2]] + 1, method= "all", include.greedy=TRUE)
  ball1 <- apply(VI$cl[2:5,], 1, function(x) credibleball(x,mod[[2]] + 1, c.dist = "VI", alpha = 0.05))
  ball2 <- apply(BI$cl[2:5,], 1, function(x) credibleball(x,mod[[2]] + 1, c.dist = "Binder", alpha = 0.05))
  
  grVI <- greedy2(psm = psm, cls.draw = mod[[2]] + 1, greedy.out = TRUE,
                  suppress.comment = TRUE, start.cl = VI$cl[4,], loss = "VI")
  grBI <- greedy2(psm = psm, cls.draw = mod[[2]] + 1, greedy.out = TRUE,
                  suppress.comment = TRUE, start.cl = BI$cl[4,], loss = "Binder")
  grballVI <- credibleball(grVI$cl,grVI$greed, c.dist = "VI", alpha = 0.05)
  grballBI <- credibleball(grVI$cl,grVI$greed, c.dist = "Binder", alpha = 0.05)
  list(mod, VI, BI, ball1, ball2, grVI, grBI, grballVI, grballBI)
}
stopCluster(cl)

################################################################################################
# AMAZING PLOT
################################################################################################

library(ggplot2)
library(ggpubr)

VImatS <- BImatS <-VImatL <- BImatL <- VImatU <- 
  BImatU <- VImatH <- BImatH <- matrix(0, ncol = 5, nrow = 10)
k <- 1
for(i in 1:10){
  for(j in 1:5){
    VImatS[i,j] <- length(table(results[[k]][[4]][[4]]$c.star))
    VImatU[i,j] <- length(table(results[[k]][[4]][[4]]$c.uppervert))
    VImatL[i,j] <- length(table(results[[k]][[4]][[4]]$c.lowervert))
    VImatH[i,j] <- length(table(results[[k]][[4]][[4]]$c.horiz))
    
    BImatS[i,j] <- length(table(results[[k]][[5]][[4]]$c.star))
    BImatU[i,j] <- length(table(results[[k]][[5]][[4]]$c.uppervert))
    BImatL[i,j] <- length(table(results[[k]][[5]][[4]]$c.lowervert))
    BImatH[i,j] <- length(table(results[[k]][[5]][[4]]$c.horiz))
    
    k <- k+1
  }
}

as.vector(t(cbind(colMeans(VImatS), colMeans(VImatU), colMeans(VImatL), colMeans(VImatH))))
VIend <- as.data.frame(cbind(log(ssize), colMeans(VImatS), colMeans(VImatU), colMeans(VImatL)))
BIend <- as.data.frame(cbind(log(ssize), colMeans(BImatS), colMeans(BImatU), colMeans(BImatL)))
BIVIend <- as.data.frame(cbind(c(rep("VI", 5), rep("BI", 5)), rbind(VIend, BIend)))
names(BIVIend)[1] <- "Gr"

col1 <- "#E69F00"
col2 <- "#009E73"
awesome_plot <- ggplot(BIVIend, aes(x = V1, y = V2, colour = Gr), col = c(col1, col2)) + 
  geom_ribbon(data = VIend, aes(ymin = V3, ymax = V4, colour = col1), fill = col1, alpha = "0.0") +
  geom_ribbon(data = BIend, aes(ymin = V3, ymax = V4, colour = col2), fill = col2, alpha = "0.0") +
  geom_line(data = BIend, aes(x = V1, y = V2), col = col2) +
  geom_line(data = VIend, aes(y = V1), col = col1) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  ylab("partition size") +
  xlab("log sample size") +
  labs(title = "Comparison between VI and Binder", subtitle = "greedy algorithm") +
  theme() +
  scale_color_manual("Loss", limits=c("VI", "Binder"), values = c(col1, col2)) +
  guides(colour = guide_legend(override.aes = list(pch = c(1,1), fill = c(col1, col2), alpha = c(1,1)))) +
  border()
awesome_plot

VImatS <- BImatS <-VImatL <- BImatL <- VImatU <- BImatU <- VImatH <- BImatH <- matrix(0, ncol = 5, nrow = 10)
k <- 1
for(i in 1:10){
  for(j in 1:5){
    VImatS[i,j] <- length(table(results[[k]][[4]][[3]]$c.star))
    VImatU[i,j] <- length(table(results[[k]][[4]][[3]]$c.uppervert))
    VImatL[i,j] <- length(table(results[[k]][[4]][[3]]$c.lowervert))
    VImatH[i,j] <- length(table(results[[k]][[4]][[3]]$c.horiz))
    
    BImatS[i,j] <- length(table(results[[k]][[5]][[3]]$c.star))
    BImatU[i,j] <- length(table(results[[k]][[5]][[3]]$c.uppervert))
    BImatL[i,j] <- length(table(results[[k]][[5]][[3]]$c.lowervert))
    BImatH[i,j] <- length(table(results[[k]][[5]][[3]]$c.horiz))
    
    k <- k+1
  }
}

as.vector(t(cbind(colMeans(VImatS), colMeans(VImatU), colMeans(VImatL), colMeans(VImatH))))
VIend <- as.data.frame(cbind(log(ssize), colMeans(VImatS), colMeans(VImatU), colMeans(VImatL)))
BIend <- as.data.frame(cbind(log(ssize), colMeans(BImatS), colMeans(BImatU), colMeans(BImatL)))
BIVIend <- as.data.frame(cbind(c(rep("VI", 5), rep("BI", 5)), rbind(VIend, BIend)))
names(BIVIend)[1] <- "Gr"
expect <- as.data.frame(cbind(log(ssize), log(ssize)))

awesome_plot2 <- ggplot(BIVIend, aes(x = V1, y = V2, colour = Gr), col = c(col1, col2)) + 
  geom_line(data = BIend, aes(x = V1, y = V2), col = col2) +
  geom_line(data = VIend, aes(y = V1), col = col1) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  ylab("partition size") +
  xlab("log sample size") +
  labs(title = "Comparison between VI and Binder", subtitle = "draw algorithm") +
  scale_color_manual("Loss", limits=c("VI", "Binder"), values = c(col1, col2)) +
  border()

pdf("planC_greedy2.pdf", width = 8, height = 3.5)
ggarrange(awesome_plot2, awesome_plot, ncol = 2, widths = c(1, 1.27))
dev.off()



