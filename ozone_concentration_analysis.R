
rm(list=ls())

library(QRegVCM)


data(airquality)



# Remove rows with NAs
airquality_clean <- na.omit(airquality[, c("Ozone", "Solar.R", "Wind", "Temp", "Month", "Day")])

# Create a new continuous variable for day
airquality_clean$DayInStudy <- 1:nrow(airquality_clean)


# Density estimation for IVs and DV
par(mfrow=c(2,2))
x1 <- airquality_clean$Ozone
fx1 = density(x1, kernel="gaussian")
plot(fx1, col="red", lwd=2,  main="Density Plot of Ozone",
     xlab="Ozone (ppb)")

x2 <- airquality_clean$Solar.R
fx2 = density(x2, kernel="gaussian")
plot(fx2, col="blue", lwd=2,  main="Density Plot of Solar Radiation",
     xlab="Solar Radiation (Ly)")

x3 <- airquality_clean$Wind 
fx3 = density(x3, kernel="gaussian")
plot(fx3, col="green", lwd=2, main="Density Plot of Wind",
     xlab = "Wind (mph)")

x4 <- airquality_clean$Temp
fx4 = density(x4, kernel="gaussian")
plot(fx4, col="purple", lwd=2, main="Density plot of Temp",
     xlab = "Temperature (°F)")

par(mfrow=c(3,1)) # I recommend "shrinking" the x axis when clicking on "Zoom", otherwise it looks 
# too stretched out.


# NW and local polynomial regression 
library("KernSmooth")

# Ozone and Solar Radiation
h1 = sd(x2)*(length(x2))^(-1/5) ## Rule of thumb bandwidth for NW
h2 <- dpill(x2, x1) ## Plug-in method for LP

plot(x2,x1, main='NP regression - Ozone in function of Solar Radiation',
     xlab = 'Solar Radiation (Ly)', ylab='Ozone (ppb)')
lines(ksmooth(x2, x1, kernel="normal", bandwidth=h1), col="red", lwd=3, lty=1)
lines(locpoly(x2, x1, degree=1,  kernel="normal", bandwidth=h2, range.x=range(x2), binned=FALSE),
      col="blue", lwd=3, lty=2)	
legend("topright", c("NW", "LP"),ncol=1, col=c("red","blue"), 
       lwd=c(2,2), lty=c(1,2), cex=0.5)


# Ozone and Wind
h1 = sd(x3)*(length(x3))^(-1/5) ## Rule of thumb bandwidth for NW
h2 <- dpill(x3, x1) ## Plug-in method for LP

plot(x3,x1, main='NP regression - Ozone in function of Wind',
     xlab = 'Wind (mph)', ylab='Ozone (ppb)')
lines(ksmooth(x3, x1, kernel="normal", bandwidth=h1), col="red", lwd=3, lty=1)
lines(locpoly(x3, x1, degree=1,  kernel="normal", bandwidth=h2, range.x=range(x3), binned=FALSE),
      col="blue", lwd=3, lty=2)	
legend("topright", c("NW", "LP"),ncol=1, col=c("red","blue"), 
       lwd=c(2,2), lty=c(1,2), cex=0.5)

# Ozone and Temperature
h1 = sd(x4)*(length(x4))^(-1/5) ## Rule of thumb bandwidth for NW
h2 <- dpill(x4, x1) ## Plug-in method for LP

plot(x4,x1, main='NP regression - Ozone in function of Temperature',
     xlab = 'Temperature (°F)', ylab='Ozone (ppb)')
lines(ksmooth(x4, x1, kernel="normal", bandwidth=h1), col="red", lwd=3, lty=1)
lines(locpoly(x4, x1, degree=1,  kernel="normal", bandwidth=h2, range.x=range(x4), binned=FALSE),
      col="blue", lwd=3, lty=2)	
legend("topright", c("NW", "LP"),ncol=1, col=c("red","blue"), 
       lwd=c(2,2), lty=c(1,2), cex=0.5)




par(mfrow=c(1,1))
# AHeVT model
y <- airquality_clean$Ozone          # Ozone concentration (response)
times <- airquality_clean$DayInStudy  # Day in the study (time variable)
subj <- airquality_clean$Month        # Month (subject/grouping identifier)
dim <- length(y)                      # Number of observations after cleaning

# Define covariates
x0 <- rep(1, dim)                       # Intercept
x1 <- airquality_clean$Solar.R          # Solar Radiation (continuous)
x2 <- airquality_clean$Wind             # Wind Speed (continuous)
x3 <- airquality_clean$Temp             # Temperature (continuous)
X <- cbind(x0, x1, x2, x3)

# Average day
#VecX <- c(1, mean(x1), mean(x2), mean(x3))

# Extreme day
VecX <- c(1, max(x1), min(x2), max(x3))


kn <- c(8, 8, 8, 8)          # Number of knots for each varying coefficient.
degree <- c(3, 3, 3, 3)      # Degree of B-spline basis
taus <- seq(0.1, 0.9, 0.1)   # Quantiles of interest
lambdas <- c(1, 1.5, 1.5, 1.5) # Smoothing parameters
d <- c(1, 1, 1, 1)           # Order of differencing operator
gam <- 1/2                   # Power used in estimating smoothing parameter


AHe_air <- AHeVT(VecX = VecX, times = times, subj = subj, X = X, y = y, d = d,
                 tau = taus, kn = kn, degree = degree, lambda = lambdas, gam = gam)

hat_bt50_air <- AHe_air$hat_bt50 
hat_VT_air <- AHe_air$hat_Vt     
qhat_air <- AHe_air$qhat         

qhat_list_air <- lapply(1:ncol(qhat_air), function(i) qhat_air[,i])
names(qhat_list_air) <- paste0("qhat", 1:9, "_air")
list2env(qhat_list_air, envir = .GlobalEnv) 


hat_bt0_air <- hat_bt50_air[seq(1, dim)]
hat_bt1_air <- hat_bt50_air[seq((dim + 1), (2 * dim))]       
hat_bt2_air <- hat_bt50_air[seq((2 * dim + 1), (3 * dim))]   
hat_bt3_air <- hat_bt50_air[seq((3 * dim + 1), (4 * dim))]  


order_indices <- order(times, hat_VT_air, qhat_air[,1],
                       hat_bt0_air, hat_bt1_air, hat_bt2_air, hat_bt3_air)

times_ordered <- times[order_indices]
hat_VT_ordered <- hat_VT_air[order_indices]
qhat_ordered_matrix <- qhat_air[order_indices, ] 
hat_bt0_ordered <- hat_bt0_air[order_indices]
hat_bt1_ordered <- hat_bt1_air[order_indices]
hat_bt2_ordered <- hat_bt2_air[order_indices]
hat_bt3_ordered <- hat_bt3_air[order_indices]


# Plot coefficient estimators
par(mfrow=c(2,1)) 

plot(hat_bt0_ordered ~ times_ordered, lwd = 2, type = "l",
     xlab = "Day in Study", ylab = "Baseline Ozone",
     main = "(A) Varying Intercept (Baseline Ozone) over Study Period");

plot(hat_bt1_ordered ~ times_ordered, lwd = 2, type = "l",
     xlab = "Day in Study", ylab = "Coefficient of Solar.R",
     main = "(B) Varying Effect of Solar Radiation over Study Period");

plot(hat_bt2_ordered ~ times_ordered, lwd = 2, type = "l",
     xlab = "Day in Study", ylab = "Coefficient of Wind",
     main = "(C) Varying Effect of Wind Speed over Study Period");

plot(hat_bt3_ordered ~ times_ordered, lwd = 2, type = "l",
     xlab = "Day in Study", ylab = "Coefficient of Temperature",
     main = "(D) Varying Effect of Temperature over Study Period");

### Plot variability V(t)
par(mfrow=c(1,1))
ylim_vt <- range(hat_VT_ordered, na.rm = TRUE)
plot(hat_VT_ordered ~ times_ordered, ylim = c(min(hat_VT_air), max(hat_VT_air)),
     xlab = "Day in Study", ylab = "", type = "l", lwd = 2,
     main = "(E) Estimated Variability V(t) of Ozone over Study Period");
mtext(expression(hat(V)(t)), side = 2, cex = 1, line = 3)


### Plot conditional quantiles estimators
library(ggplot2)
library(tidyr) 
library(dplyr) 

plot_data_gg <- as.data.frame(qhat_ordered_matrix)
colnames(plot_data_gg) <- paste0("Tau_", taus) 
plot_data_gg$DayInStudy <- times_ordered      


plot_data_gg_long <- plot_data_gg %>%
  pivot_longer(cols = starts_with("Tau_"),
               names_to = "Quantile_Level",
               names_prefix = "Tau_",
               values_to = "Ozone_Concentration") %>%
  mutate(Quantile_Level = as.factor(as.numeric(Quantile_Level))) 


gg_colors <- c("0.1" = "magenta", "0.2" = "gray50", "0.3" = "blue",
               "0.4" = "brown",   "0.5" = "black",  "0.6" = "orange",
               "0.7" = "darkcyan","0.8" = "green",  "0.9" = "red")
gg_lty <-    c("0.1" = "dashed", "0.2" = "dotted", "0.3" = "longdash", 
               "0.4" = "dotdash","0.5" = "solid",  "0.6" = "dotdash",
               "0.7" = "longdash","0.8" = "dotted", "0.9" = "dashed") 
gg_lwd <-    c("0.1" = 0.8, "0.2" = 0.8, "0.3" = 0.8, 
               "0.4" = 0.8, "0.5" = 1.2, "0.6" = 0.8, #
               "0.7" = 0.8, "0.8" = 0.8, "0.9" = 1.2)


p_legend_right <- ggplot(plot_data_gg_long, aes(x = DayInStudy, y = Ozone_Concentration,
                                                group = Quantile_Level,
                                                color = Quantile_Level,
                                                linetype = Quantile_Level,
                                                linewidth = Quantile_Level)) +
  geom_line() +
  scale_color_manual(values = gg_colors, name = "Tau") +
  scale_linetype_manual(values = gg_lty, name = "Tau") +
  scale_linewidth_manual(values = gg_lwd, name = "Tau") +
  labs(title = "(G) Conditional Quantile Estimators of Ozone (extreme day)",
       x = "Day in Study",
       y = "Ozone Concentration (ppb)") +
  theme_minimal() + 
  theme(
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"), 
    legend.text = element_text(size = 8),  
    legend.title = element_text(size = 9, face = "bold"), 
    legend.background = element_rect(fill="NA", colour = "NA"),
    legend.box.margin = margin(t = 0, r = 5, b = 0, l = 0, unit = "pt"),
    guides(linetype = "none", 
           linewidth = "none", 
           color = guide_legend(title = "Tau")) 
  )
print(p_legend_right)





