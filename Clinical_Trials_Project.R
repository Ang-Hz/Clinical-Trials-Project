#### Clinical Trials Project ####
#
# This project evaluates the understanding and application of core Phase I–II methodologies, randomization, and group sequential designs.
#

#### Phase I methodology ####
# The purpose is to study the behavior of the "3+3" algorithm under different toxicity profiles.
##

# Defining a custom function to calculate the operating characteristics of a standard 3+3 design
OC_3plus3 <- function(p, doses = NULL) {
  # If no dose names are provided, generate default labels (Dose1, Dose2, etc.)
  if(is.null(doses)) {
    doses <- paste0("Dose", seq_along(p))
  }
  
  # Pi: Probability of escalation at dose level i.
  # Calculated as: P(0 DLTs in 1st cohort of 3) + [P(1 DLT in 1st cohort of 3) * P(0 DLTs in 2nd cohort of 3)]
  Pi <- dbinom(0, size=3, prob = p) + 
        dbinom(1, size=3, prob = p) * dbinom(0, size=3, prob = p)
        
  # Qi: Cumulative probability of escaping/escalating through all doses up to level i
  Qi <- cumprod (Pi)
  
  # OC: Probability of stopping the trial at or before dose i (1 - Qi)
  OC <- 1 - Qi
  
  # Return results compiled nicely into a data frame
  data.frame(
    dose = doses,
    p = p,
    Pi = Pi,
    Qi = Qi,
    OC = OC
  )
}

# Scenarios of true toxicity probabilities (p) across 4 dose levels
p_A <- c(0.08, 0.18, 0.30, 0.45)
p_B <- c(0.12, 0.22, 0.28, 0.40)


## Transition Probabilities (Pi) calculation for both scenarios
res_A <- OC_3plus3(p_A, doses = paste0("Dose ", 1:4))
res_A
res_B <- OC_3plus3(p_B, doses = paste0("Dose ", 1:4))
res_B


## Probability of early stopping before reaching Dose 4
## early stopping A
early_A <- 1 - res_A[3, "Qi"]
early_A

## early stopping B
early_B <- 1 - res_B[3, "Qi"]
early_B

# Final Dose Recommendation
## Scenario A calculation of MTD (Maximum Tolerated Dose)
# Shift the escalation probabilities to find the next dose's escalation probability
P_next_A <- c(res_A$Pi[-1], NA)

# Si_A: Probability that dose level i is selected as the recommended phase II dose
res_A$Si_A <- res_A$Qi * (1 - P_next_A)
res_A$Si_A[4] <- res_A[4, "Qi"] # 4th Dose boundary condition (cannot escalate further)
res_A$Si_A 

# Scenario B
# Shift the escalation probabilities for Scenario B
P_next_B <- c(res_B$Pi[-1], NA)

# Si_B: Probability of dose recommendation for Scenario B
res_B$Si_B <- res_B$Qi * (1 - P_next_B)
res_B$Si_B[4] <- res_B[4, "Qi"] 
res_B$Si_B

# Explanatory Table
# Compiling and formatting the final recommendation probabilities for both scenarios side-by-side
combined_table <- data.frame(
  "Dosage" = res_A$dose,
  
  # Results for A
  "Probability_A" = round(res_A$Si_A, 4),
  "Percentage_A" = paste0(round(res_A$Si_A * 100, 1), "%"),
  
  # Results for B
  "Probability_B" = round(res_B$Si_B, 4), 
  "Percentage_B" = paste0(round(res_B$Si_B * 100, 1), "%")
)

# Displaying the Table 
View(combined_table, "Final Recommended Dose")


## Consistency Check
# Assessing how a global increase of +0.05 in toxicity probabilities affects the metrics
p_Ac <- p_A + 0.05
p_Bc <- p_B + 0.05
  
res_Ac <- OC_3plus3(p_Ac, doses = paste0("Dose ", 1:4))
res_Bc <- OC_3plus3(p_Bc, doses = paste0("Dose ", 1:4))


#### Phase II methodology ####
##Single-Stage Exact Binomial Design

# Installing and loading the clinical trial design package 'clinfun'
#install.packages("clinfun")
library(clinfun)

# Setting up design parameters for Scenario A
p0_A <- 0.15    # unacceptable/historical response level (Null Hypothesis efficacy)
p1_A <- 0.40    # desirable/clinically meaningful level (Alternative Hypothesis efficacy)
alpha_A <- 0.10 # one-sided Type I error rate
beta_A <- 0.20  # Type II error rate (power = 1 - beta = 0.80)

# Finding exact binomial single-stage designs matching the criteria
out_A <- ph2single(pu = p0_A, pa = p1_A, ep1 = alpha_A, ep2 = beta_A, nsoln = 10)
out_A                  

# Sorting the potential solutions to identify the one with the minimal sample size
out_smallest_n_A <- out_A[order(out_A$n), ][1, ]
out_smallest_n_A # the smallest sample size solution for Scenario A

# Setting up design parameters for Scenario B
p0_B <- 0.20
p1_B <- 0.40
alpha_B <- 0.05
beta_B <- 0.10

# Finding exact binomial single-stage designs for Scenario B
out_B <- ph2single(pu = p0_B, pa = p1_B, ep1 = alpha_B, ep2 = beta_B, nsoln = 10)
out_B  

# Sorting to find the minimal sample size solution for Scenario B
out_smallest_n_B <- out_B[order(out_B$n), ][1, ]
out_smallest_n_B


##Simon's Two-Stage Designs (Optimal and Minimax)

# Generating Simon's 2-stage designs for Scenario A (Optimal minimize EN0, Minimax minimize max N)
sim_design_A <- ph2simon(pu = p0_A, pa = p1_A,
                         ep1 = alpha_A, ep2 = beta_A,
                         nmax = 100)
sim_design_A

# Generating Simon's 2-stage designs for Scenario B
sim_design_B <- ph2simon(pu = p0_B, pa = p1_B,
                           ep1 = alpha_B, ep2 = beta_B,
                           nmax = 100)
sim_design_B


#### Randomization and Imbalance ####
##Theoretical Probability of Perfect Balance
N <- 80
# Calculation of P(N_A = 40) using the exact binomial density formula
dbinom(x = N/2, size = N, p = 0.5)


##Simple Randomization Simulation
set.seed(222) # for reproducibility 
K <- 1000     # number of replications

# Storage of results in a structured dataframe
results <- data.frame(
  replication = 1:K,
  NumberA = NA_integer_
)
  
# Executing Loop to simulate K clinical trials using simple randomization
for (i in 1:K) {
  u <- runif(N)
  results$NumberA[i] <- sum(u < 0.5) # 50% probability for each patient to be assigned to group A
}

# Finding the 2.5th and 97.5th percentiles to build the empirical interval
conf_interval_sim <- quantile(results$NumberA, probs = c(0.025, 0.975))
print(conf_interval_sim) # 95% C.I. of patients who ended up in group A 

# lower limit in proportion
L_CI <- 31 / 80
L_CI

# upper limit in proportion
U_CI <- 49 / 80
U_CI

# Calculating how often the trial results in a severe balance allocation worse than 60/40 or 40/60
imbalance_rate <- mean(results$NumberA > 48 | results$NumberA < 32)
#Imbalance Rate
cat("Imbalance Rate > 60/40:", imbalance_rate * 100, "%\n")


##Permuted Block Randomization
#install.packages("blockrand")
library(blockrand)
set.seed(333)

# Creating a block randomization schedule
block_rand <- blockrand(
  n = 80,
  num.levels = 2,
  levels = c("A", "B"),
  block.sizes = 2 # 2*2 levels = block size 4
)

# How many rows did it return in total? (Accounts for full blocks generated past N=80)
nrow(block_rand)

# Keeping exactly the first 80 patients as required by the protocol
block_rand_80 <- block_rand[1:80, ]

# Distribution of A/B in the first 80 individuals to check for balance
table(block_rand_80$treatment)

# Expressing the absolute distribution as relative proportions
prop.table(table(block_rand_80$treatment))

# Displaying the first 20 lines/blocks of the allocation schedule
head(block_rand, 20)


##Empirical Comparison of Randomization Methods
set.seed(444)
N <- 14
it <- 200
results <- data.frame(
  SR = rep(NA_real_, it),
  BR = rep(NA_real_, it)
)

# Simulating 200 trials of small sample size (N=14) to see the variance in allocation balance
for (i in 1:it) {
  # SR: Simple Randomization allocation
  trt_sr <- ifelse(runif(N) < 0.5, "A", "B")
  results$SR[i] <- mean(trt_sr == "A")
  
  # BR: Permuted Block Randomization allocation (dynamic block size selection up to 4)
  br <- blockrand(n = N, num.levels = 2, levels = c("A", "B"),
                  block.sizes = 1:4)[1:N, ]
  results$BR[i] <- mean(br$treatment == "A")
}
colMeans(results)  # Checking average allocation proportion (aiming close to 0.5)
apply(results, 2, sd) # Evaluating the standard deviation (variability/risk of imbalance)


#### Interim analysis, O’Brien–Fleming and Pocock ####
##Group Sequential Boundaries (Interim Analysis)

install.packages("gsDesign")
library(gsDesign)

# Pocock Boundaries generation (constant critical thresholds across all looks)
xA = gsDesign(k = 3, timing = c(0.30, 0.60, 1.00), test.type = 2, alpha = 0.025, sfu = "Pocock")
gsBoundSummary(xA)

# Extracting and converting Pocock Z-scores and 2-sided p-values per look
cbind(Z = gsBoundSummary(xA)$Eff[gsBoundSummary(xA)$Val=="Z"],
      P = 2*gsBoundSummary(xA)$Eff[gsBoundSummary(xA)$Val=="p (1-sided)"])

# O'Brien-Fleming boundaries generation (conservative initial boundaries, easier to stop later)
xB = gsDesign(k = 3, timing = c(0.30, 0.60, 1.00), test.type = 2, alpha = 0.025, sfu = "OF")
gsBoundSummary(xB)

# Extracting and converting O'Brien-Fleming Z-scores and 2-sided p-values per look
cbind(Z = gsBoundSummary(xB)$Eff[gsBoundSummary(xB)$Val=="Z"],
      P = 2*gsBoundSummary(xB)$Eff[gsBoundSummary(xB)$Val=="p (1-sided)"])
