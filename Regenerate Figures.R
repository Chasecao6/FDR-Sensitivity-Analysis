library("lattice")
library("zipfR") 
library(ggplot2)
library(cowplot)
library(patchwork)
library(ggpubr)
library(ggtext)



#########################For Figure 1 (step function SPM)
### Bias function for step function SPM 
Biastep=function(power,b,lambda,alpha,pi){
  gamma=log(power)/log(alpha)
  d=pi*(1+b)*alpha/2+(1-pi)*((alpha/2)^gamma+b*(alpha^gamma-(alpha/2)^gamma))
  co=alpha/(alpha-lambda)
  if (lambda>=alpha/2){
    n=co*(pi*b*(alpha-lambda)+(1-pi)*b*(alpha^gamma-lambda^gamma))-pi*(1+b)*alpha/2
    
  }else{
    n=co*(pi*(b*alpha/2+alpha/2-lambda)+(1-pi)*(b*(alpha^gamma-(alpha/2)^gamma)+(alpha/2)^gamma-lambda^gamma))-pi*(1+b)*alpha/2
  }
  bias=n/d
  
  return(bias)
}

# Function to create the plot data
generate_plot_data <- function(lambda, pi) {
  x <- seq(0.1, 1, length.out = 101)  # power
  y <- seq(0, 1, length.out = 101)  # publication bias parameter
  
  data <- expand.grid(X = x, Y = y)
  data$Z <- mapply(Biastep, data$X, data$Y, lambda = lambda, alpha = 0.05, pi = pi)
  
  # Use Lambda and Pi as factors for faceting
  data$Lambda <- as.factor(lambda)
  data$Pi <- as.factor(pi)
  
  return(data)
}

# List of parameter settings
params <- list(
  list(lambda = 0.01, pi = 0.8),
  list(lambda = 0.025, pi = 0.8),
  list(lambda = 0.045, pi = 0.8),
  list(lambda = 0.01, pi = 0.5),
  list(lambda = 0.025, pi = 0.5),
  list(lambda = 0.045, pi = 0.5),
  list(lambda = 0.01, pi = 0.3),
  list(lambda = 0.025, pi = 0.3),
  list(lambda = 0.045, pi = 0.3)
)

# Generate plot data for each parameter set
plot_data_list <- lapply(params, function(p) generate_plot_data(p$lambda, p$pi))

# Combine all data into a single data frame
plot_data <- do.call(rbind, plot_data_list)



# Plotting function with Greek symbols in facet labels
plot_bp_ggplot <- function() {
  ggplot(plot_data, aes(x = X, y = Y, fill = Z)) +
    geom_tile() +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue",
                         midpoint = 0, limits = c(-0.7, 0.8),
                         breaks = seq(-1, 1, by = 0.1)) +
    labs(x = "Power", y = bquote(rho), fill = "Bias") +
    theme_minimal() +
    facet_grid(
      rows = vars(PiLabel),
      cols = vars(LambdaLabel),
      labeller = label_parsed
    ) +
    theme(legend.position = "right")
}

plot_data$PiLabel <- with(plot_data, paste0("pi[0]==", Pi))
plot_data$LambdaLabel <- with(plot_data, paste0("lambda==", Lambda))
# Set correct order for facets
plot_data$PiLabel <- factor(plot_data$PiLabel, levels = c("pi[0]==0.8", "pi[0]==0.5", "pi[0]==0.3"))

combined_plot <- plot_bp_ggplot() +
  plot_annotation(caption = "", theme = theme(plot.caption = element_text(hjust = 0.5)))

print(combined_plot)








########################################### For Figure 2 (Beta density SPM)
### Bias function for Beta density SPM 
Bias=function(power,b,lambda,alpha,pi){
  gamma=log(power)/log(alpha)
  co=alpha/(alpha-lambda)
  n=(pi/b)*((1-lambda)^b-(1-alpha)^b)+(1-pi)*gamma*(Ibeta(alpha,gamma,b)-Ibeta(lambda,gamma,b))
  d=(pi/b)*(1-(1-alpha)^b)+(1-pi)*gamma*Ibeta(alpha,gamma,b)
  
  bias=(co*n-(pi/b)*(1-(1-alpha)^b))/d
  
  return(bias)
}
# Function to create the plot data
generate_plot_data <- function(lambda, pi) {
  x <- seq(0.1, 1, length.out = 101)  # power
  y <- seq(1, 40, length.out = 101)  # publication bias parameter
  
  data <- expand.grid(X = x, Y = y)
  data$Z <- mapply(Bias, data$X, data$Y, lambda = lambda, alpha = 0.05, pi = pi)
  
  # Use Lambda and Pi as factors for faceting
  data$Lambda <- as.factor(lambda)
  data$Pi <- as.factor(pi)
  
  return(data)
}

# List of parameter settings
params <- list(
  list(lambda = 0.01, pi = 0.8),
  list(lambda = 0.025, pi = 0.8),
  list(lambda = 0.045, pi = 0.8),
  list(lambda = 0.01, pi = 0.5),
  list(lambda = 0.025, pi = 0.5),
  list(lambda = 0.045, pi = 0.5),
  list(lambda = 0.01, pi = 0.3),
  list(lambda = 0.025, pi = 0.3),
  list(lambda = 0.045, pi = 0.3)
)

# Generate plot data for each parameter set
plot_data_list <- lapply(params, function(p) generate_plot_data(p$lambda, p$pi))

# Combine all data into a single data frame
plot_data <- do.call(rbind, plot_data_list)



# Plotting function with Greek symbols in facet labels
plot_bp_ggplot <- function() {
  ggplot(plot_data, aes(x = X, y = Y, fill = Z)) +
    geom_tile() +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue",
                         midpoint = 0, limits = c(-0.7, 0.8),
                         breaks = seq(-1, 1, by = 0.1)) +
    labs(x = "Power", y = bquote(eta), fill = "Bias") +
    theme_minimal() +
    facet_grid(
      rows = vars(PiLabel),
      cols = vars(LambdaLabel),
      labeller = label_parsed
    ) +
    theme(legend.position = "right")
}

plot_data$PiLabel <- with(plot_data, paste0("pi[0]==", Pi))
plot_data$LambdaLabel <- with(plot_data, paste0("lambda==", Lambda))
# Set correct order for facets
plot_data$PiLabel <- factor(plot_data$PiLabel, levels = c("pi[0]==0.8", "pi[0]==0.5", "pi[0]==0.3"))

combined_plot <- plot_bp_ggplot() +
  plot_annotation(caption = "", theme = theme(plot.caption = element_text(hjust = 0.5)))

print(combined_plot)
