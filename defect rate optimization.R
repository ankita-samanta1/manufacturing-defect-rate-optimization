install.packages("leaps")
library(leaps)

# Load dataset
df <- read.csv("C:\\Users\\user\\Downloads\\manufacturing_defect_dataset.csv")

# Best subset selection method
best_fit <- regsubsets(DefectRate ~ MaintenanceHours + ProductionVolume + AdditiveProcessTime + 
                         AdditiveMaterialCost+ProductionCost+SupplierQuality+DeliveryDelay+
                         QualityScore+DowntimePercentage+InventoryTurnover+DefectStatus+
                         StockoutRate+WorkerProductivity+SafetyIncidents+EnergyConsumption+EnergyEfficiency
                       ,
                       data = df,
                       nbest = 1)   

# Summary 
summary_best <- summary(best_fit)
summary_best
which.max(summary_best$adjr2)
summary_best$which[which.max(summary_best$adjr2),]

###############################################################
##------------------------------------------------------------
## UNCONSTRAINED OPTIMIZATION USING GOLDEN SECTION METHOD
## Objective: Minimize DefectRate with respect to MaintenanceHours
##------------------------------------------------------------
###############################################################

## Load dataset
df <- read.csv("C:\\Users\\STUDENT\\Downloads\\manufacturing_defect_dataset.csv")

## Step 1: Fit quadratic regression model
model <- lm(DefectRate ~ MaintenanceHours + I(MaintenanceHours^2), data = df)
summary(model)

## Step 2: Define the fitted objective function f(x)
f <- function(x) {
  beta <- coef(model)
  return(beta[1] + beta[2]*x + beta[3]*x^2)
}

## Step 3: Golden Section Search implementation
golden_section <- function(f, a, b, tol = 1e-5, max_iter = 1000) {
  gr <- (sqrt(5) - 1) / 2  # Golden ratio constant
  iter <- 1
  
  # Initialize points
  x1 <- b - gr * (b - a)
  x2 <- a + gr * (b - a)
  
  f1 <- f(x1)
  f2 <- f(x2)
  
  while (abs(b - a) > tol && iter < max_iter) {
    if (f1 > f2) {
      a <- x1
      x1 <- x2
      f1 <- f2
      x2 <- a + gr * (b - a)
      f2 <- f(x2)
    } else {
      b <- x2
      x2 <- x1
      f2 <- f1
      x1 <- b - gr * (b - a)
      f1 <- f(x1)
    }
    iter <- iter + 1
  }
  
  xmin <- (a + b) / 2
  fmin <- f(xmin)
  
  cat("Converged in", iter, "iterations\n")
  cat("Minimum Defect Rate =", round(fmin, 5), "at Maintenance Hours =", round(xmin, 5), "\n")
  
  return(list(xmin = xmin, fmin = fmin))
}

## Step 4: Set search interval (within realistic data range)
a <- min(df$MaintenanceHours)
b <- max(df$MaintenanceHours)

## Step 5: Run Golden Section Search
result <- golden_section(f, a, b)

###############################################################
##------------------------------------------------------------
## UNCONSTRAINED OPTIMIZATION USING STEEPEST DESCENT METHOD
## WITH EXACT LINE SEARCH (via Golden Section Search)
## Objective: Minimize DefectRate 
## w.r.t MaintenanceHours, ProductionVolume, AdditiveProcessTime, QualityScore
##------------------------------------------------------------
###############################################################

## Load dataset
df <- read.csv("C:\\Users\\user\\Downloads\\manufacturing_defect_dataset.csv")

## Step 1: Fit quadratic regression model
model <- lm(DefectRate ~ MaintenanceHours + ProductionVolume + AdditiveProcessTime + QualityScore +
              I(MaintenanceHours^2) + I(ProductionVolume^2) + 
              I(AdditiveProcessTime^2) + I(QualityScore^2) +
              I(MaintenanceHours * ProductionVolume) +
              I(MaintenanceHours * AdditiveProcessTime) +
              I(MaintenanceHours * QualityScore) +
              I(ProductionVolume * AdditiveProcessTime) +
              I(ProductionVolume * QualityScore) +
              I(AdditiveProcessTime * QualityScore),
            data = df)

summary(model)

## Extract coefficients
beta <- coef(model)
b0 <- beta[1]

## Step 2: Define objective function
f <- function(x1, x2, x3, x4) {
  b0 +
    beta["MaintenanceHours"] * x1 +
    beta["ProductionVolume"] * x2 +
    beta["AdditiveProcessTime"] * x3 +
    beta["QualityScore"] * x4 +
    beta["I(MaintenanceHours^2)"] * x1^2 +
    beta["I(ProductionVolume^2)"] * x2^2 +
    beta["I(AdditiveProcessTime^2)"] * x3^2 +
    beta["I(QualityScore^2)"] * x4^2 +
    beta["I(MaintenanceHours * ProductionVolume)"] * x1 * x2 +
    beta["I(MaintenanceHours * AdditiveProcessTime)"] * x1 * x3 +
    beta["I(MaintenanceHours * QualityScore)"] * x1 * x4 +
    beta["I(ProductionVolume * AdditiveProcessTime)"] * x2 * x3 +
    beta["I(ProductionVolume * QualityScore)"] * x2 * x4 +
    beta["I(AdditiveProcessTime * QualityScore)"] * x3 * x4
}

## Step 3: Define gradient
grad_f <- function(x1, x2, x3, x4) {
  df_dx1 <- beta["MaintenanceHours"] + 2 * beta["I(MaintenanceHours^2)"] * x1 +
    beta["I(MaintenanceHours * ProductionVolume)"] * x2 +
    beta["I(MaintenanceHours * AdditiveProcessTime)"] * x3 +
    beta["I(MaintenanceHours * QualityScore)"] * x4
  
  df_dx2 <- beta["ProductionVolume"] + 2 * beta["I(ProductionVolume^2)"] * x2 +
    beta["I(MaintenanceHours * ProductionVolume)"] * x1 +
    beta["I(ProductionVolume * AdditiveProcessTime)"] * x3 +
    beta["I(ProductionVolume * QualityScore)"] * x4
  
  df_dx3 <- beta["AdditiveProcessTime"] + 2 * beta["I(AdditiveProcessTime^2)"] * x3 +
    beta["I(MaintenanceHours * AdditiveProcessTime)"] * x1 +
    beta["I(ProductionVolume * AdditiveProcessTime)"] * x2 +
    beta["I(AdditiveProcessTime * QualityScore)"] * x4
  
  df_dx4 <- beta["QualityScore"] + 2 * beta["I(QualityScore^2)"] * x4 +
    beta["I(MaintenanceHours * QualityScore)"] * x1 +
    beta["I(ProductionVolume * QualityScore)"] * x2 +
    beta["I(AdditiveProcessTime * QualityScore)"] * x3
  
  return(c(df_dx1, df_dx2, df_dx3, df_dx4))
}

## Step 4: Golden Section Search for exact line search
golden_section <- function(phi, a = 0, b = 1, tol = 1e-3) {
  gr <- (sqrt(5) - 1) / 2
  x1 <- b - gr * (b - a)
  x2 <- a + gr * (b - a)
  f1 <- phi(x1)
  f2 <- phi(x2)
  
  while (abs(b - a) > tol) {
    if (f1 > f2) {
      a <- x1
      x1 <- x2
      f1 <- f2
      x2 <- a + gr * (b - a)
      f2 <- phi(x2)
    } else {
      b <- x2
      x2 <- x1
      f2 <- f1
      x1 <- b - gr * (b - a)
      f1 <- phi(x1)
    }
  }
  return((a + b) / 2)
}

## Step 5: Steepest Descent with detailed iteration logging
steepest_descent <- function(f, grad_f, x_init, tol = 1e-3, max_iter = 500) {
  x <- x_init
  history <- data.frame(iter = 0, x1 = x[1], x2 = x[2], x3 = x[3], x4 = x[4], fval = f(x[1], x[2], x[3], x[4]))
  
  for (k in 1:max_iter) {
    g <- grad_f(x[1], x[2], x[3], x[4])
    
    ## Stop if gradient is small
    if (sqrt(sum(g^2)) < tol) {
      cat("\n✅ Converged in", k, "iterations\n")
      break
    }
    
    ## Define phi(alpha)
    phi <- function(alpha) {
      f(x[1] - alpha*g[1], x[2] - alpha*g[2], x[3] - alpha*g[3], x[4] - alpha*g[4])
    }
    
    ## Find optimal alpha using Golden Section Search
    alpha_opt <- golden_section(phi, a = 0, b = 1)
    
    ## Update x
    x_new <- x - alpha_opt * g
    f_new <- f(x_new[1], x_new[2], x_new[3], x_new[4])
    
    ## Record progress
    history <- rbind(history, data.frame(iter = k, x1 = x_new[1], x2 = x_new[2], x3 = x_new[3], x4 = x_new[4], fval = f_new))
    
    ## Print iteration details
    cat("-----------------------------------------------------------\n")
    cat("Iteration:", k, "\n")
    cat("Gradient:", paste(round(g, 6), collapse = ", "), "\n")
    cat("Alpha (step size):", round(alpha_opt, 6), "\n")
    cat("MaintenanceHours:", round(x_new[1], 6), 
        "| ProductionVolume:", round(x_new[2], 6), 
        "| AdditiveProcessTime:", round(x_new[3], 6),
        "| QualityScore:", round(x_new[4], 6), "\n")
    cat("Predicted DefectRate:", round(f_new, 8), "\n")
    cat("-----------------------------------------------------------\n")
    
    ## Check convergence
    if (sqrt(sum((x_new - x)^2)) < tol) {
      cat("\n✅ Converged in", k, "iterations\n")
      break
    }
    
    ## Update for next iteration
    x <- x_new
  }
  
  return(list(optimum = x, min_value = f(x[1], x[2], x[3], x[4]), history = history))
}

## Step 6: Run optimization
x_init <- c(mean(df$MaintenanceHours), mean(df$ProductionVolume),
            mean(df$AdditiveProcessTime), mean(df$QualityScore))

result <- steepest_descent(f, grad_f, x_init)

cat("\n================= FINAL RESULTS =================\n")
cat("Optimal MaintenanceHours:", round(result$optimum[1], 6), "\n")
cat("Optimal ProductionVolume:", round(result$optimum[2], 6), "\n")
cat("Optimal AdditiveProcessTime:", round(result$optimum[3], 6), "\n")
cat("Optimal QualityScore:", round(result$optimum[4], 6), "\n")
cat("Minimum Predicted DefectRate:", round(result$min_value, 8), "\n")


