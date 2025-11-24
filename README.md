# nmfkc: Non-negative Matrix Factorization with Kernel Covariates

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![GitHub version](https://img.shields.io/github/r-package/v/ksatohds/nmfkc)](https://github.com/ksatohds/nmfkc)
**nmfkc** is an R package that extends Non-negative Matrix Factorization (NMF) by incorporating **covariates** using kernel methods. It supports advanced features like rank selection via cross-validation, time-series modeling (NMF-VAR), and supervised classification (NMF-LAB).

# Installation

You can install the **nmfkc** package directly from [GitHub](https://github.com/ksatohds/nmfkc).

## Quick Installation (No Vignettes)

For the fastest installation without building the package vignettes (tutorials), use the following commands.

```r
# If you don't have the 'remotes' package installed, uncomment the line below:
# install.packages("remotes")

# Install the package without building vignettes (default behavior)
remotes::install_github("ksatohds/nmfkc")

# Load the package
library(nmfkc)
```

## Detailed Installation (With Vignettes)

Recommended if you want to view the tutorials locally.

```r
# Install the package and build the vignettes
remotes::install_github("ksatohds/nmfkc", build_vignettes = TRUE)

# Load the package
library(nmfkc)

# You can view the vignettes with:
browseVignettes("nmfkc")
```

## Alternative Download (Source Archive)

The Package Archive File (.tar.gz) is also available on [Dropbox](https://www.dropbox.com/scl/fo/1jxzzfqz9xsl3wlzwm3pj/ACYDNzn-h54VqgAngTKVyc0?rlkey=i5bcup2qxzwqgles27ke0f9ab&st=mdbkongw&dl=0).

## Package Explanation Video

A short video explaining the package can be found here: [YouTube](https://youtu.be/T48pPGReVu4)

# Help and Usage

```r
browseVignettes("nmfkc")
ls("package:nmfkc")
?nmfkc
```

# Citation

Please use the following command to cite the **nmfkc** package:

```r
library(nmfkc)
citation("nmfkc")
```

# Comparison with Standard NMF

| Feature | Standard NMF | nmfkc |
| :--- | :--- | :--- |
| **Handles covariates** | No | **Yes** |
| **Classification** | No | **Yes** (NMF-LAB) |
| **Time series modeling** | No | **Yes** (NMF-VAR) |
| **Nonlinearity** | No | **Yes** (Kernel) |
| **Clustering support** | Limited | **Yes** (Hard/Soft) |
| **Element-wise CV** | No | **Yes** |






# Functions

The **nmfkc** package provides a comprehensive suite of functions for performing, diagnosing, and visualizing Non-negative Matrix Factorization with Kernel Covariates.

## 1. Core Algorithm
The heart of the package is the **nmfkc** function, which estimates the basis matrix $X$ and parameter matrix $C$ (where $B=CA$) based on the observation $Y$ and covariates $A$.

- **nmfkc**: Fits the NMF model ($Y \approx XCA$) using multiplicative update rules. Supports missing values (NA) and observation weights.

## 2. Model Selection & Diagnostics
Tools to determine the optimal number of bases (Rank $Q$) and evaluate model performance.

- **nmfkc.rank**: **(Recommended)** Main diagnostic tool for rank selection. Automatically suggests the optimal rank by comparing:
    - **Elbow Method**: Based on the curvature of $R^2$.
    - **Wold's CV**: Based on Element-wise Cross-Validation (Min RMSE).
    - Also computes stability metrics like Cophenetic Correlation Coefficient (CPCC) and Silhouette scores.
- **nmfkc.ecv**: Performs **Element-wise Cross-Validation** (Wold's CV). Randomly masks elements of the matrix to evaluate structural reconstruction error. Theoretically robust for rank selection.
- **nmfkc.cv**: Performs **Column-wise Cross-Validation**. Useful for evaluating the prediction performance for new (unseen) individuals.
- **nmfkc.residual.plot**: Visualizes the original matrix, fitted matrix, and residual matrix side-by-side to check for randomness in errors.

## 3. Covariate Engineering
Functions to construct the covariate matrix $A$ from raw features $U$.

- **nmfkc.kernel**: Generates a kernel matrix (e.g., Gaussian, Polynomial) from covariates.
- **nmfkc.ar**: Constructs lagged matrices for **Vector Autoregression (NMF-VAR)** models.
- **nmfkc.class**: Converts categorical vectors into one-hot encoded matrices for classification tasks (NMF-LAB).

## 4. Parameter Tuning
Helper functions to optimize hyperparameters other than rank.

- **nmfkc.kernel.beta.cv**: Optimizes the Gaussian kernel width ($\beta$) via cross-validation.
- **nmfkc.kernel.beta.nearest.med**: Heuristically estimates $\beta$ using the median of nearest-neighbor distances.
- **nmfkc.ar.degree.cv**: Selects the optimal lag order ($D$) for NMF-VAR models.

## 5. Prediction & Forecasting
- **predict.nmfkc**: Predicts values or classes for new covariate data ($A_{new}$).
- **nmfkc.ar.predict**: Performs multi-step ahead forecasting for Time Series models.
- **nmfkc.ar.stationarity**: Checks the stationarity (stability) of the estimated VAR system.

## 6. Visualization (Graphviz/DOT)
Generates DOT scripts to visualize relationships between variables.

- **nmfkc.DOT**: Visualizes the structure of standard NMF or Kernel NMF ($A \to X \to Y$).
- **nmfkc.ar.DOT**: Visualizes the Granger causality network in Time Series models.

## 7. Utilities
- **nmfkc.normalize** / **nmfkc.denormalize**: Helper functions to scale data to $[0,1]$ and back.




# Statistical Model

The **nmfkc** package builds upon the standard NMF framework by incorporating external information (covariates). The statistical modeling process involves three steps:

1.  **Standard NMF**:
    Let $Y$ be a $P \times N$ observation matrix. Standard Non-negative Matrix Factorization approximates $Y$ as the product of two non-negative matrices:
    $$Y \approx X B$$
    where $X$ is the **basis matrix** (capturing latent patterns) and $B$ is the **coefficient matrix** (weights for each observation). Unlike ordinary linear models where $X$ is known, NMF optimizes both $X$ and $B$.

2.  **NMF with Covariates**:
    We assume that the coefficient matrix $B$ is not arbitrary but driven by known covariates $A$ (e.g., time, location, or attributes). We model this relationship as:
    $$B \approx C A$$
    where $C$ is a parameter matrix mapping covariates to the latent bases.

3.  **Tri-Factorization Model**:
    Substituting $B$ yields the core model of **nmfkc**:
    $$Y \approx X C A$$
    * **Kernel Extension**: By constructing $A$ using a kernel function (e.g., Gaussian kernel) on raw features, the model can capture non-linear relationships (Satoh, 2024).
    * **Theoretical Roots**: This formulation aligns with Orthogonal Tri-NMF (Ding et al., 2006) and generalizes the Growth Curve Models (Potthoff and Roy, 1964).

# Matrices

The goal of **nmfkc** is to optimize the unknown matrices $X$ and $C$, given the observation $Y$ and covariates $A$, minimizing the reconstruction error:

$$Y(P,N) \approx X(P,Q) \times C(Q,R) \times A(R,N)$$

### Given (Input)
-   **$Y(P,N)$**: Observation matrix ($P$ features $\times$ $N$ samples). Missing values are supported.
-   **$A(R,N)$**: Covariate matrix ($R$ covariates $\times$ $N$ samples).
    -   Can be created using `nmfkc.kernel` (for kernel methods) or `nmfkc.ar` (for time series).
    -   If no covariates are provided, $A$ defaults to the identity matrix (Standard NMF).

### Optimized (Output)
-   **$X(P,Q)$**: Basis matrix ($P$ features $\times$ $Q$ bases).
    -   $Q$ is the rank (number of bases), constrained by $Q \le \min(P,N)$.
-   **$C(Q,R)$**: Parameter matrix ($Q$ bases $\times$ $R$ covariates).
    -   This matrix links the covariates to the latent structure.
    -   Corresponds to $\Theta$ in Satoh (2024).

### Derived
-   **$B(Q,N)$**: Coefficient matrix, calculated as $B = C A$.
    -   Represents the activation of each basis for each sample.
    -   Normalized columns of $B$ (proportions) can be used for **soft clustering** probabilities.
    




# Source

-   Satoh, K. (2024) Applying Non-negative Matrix Factorization with Covariates to the Longitudinal Data as Growth Curve Model. arXiv preprint arXiv:2403.05359. <https://arxiv.org/abs/2403.05359>
-   Satoh, K. (2025) Applying non-negative Matrix Factorization with Covariates to Multivariate Time Series Data as a Vector Autoregression Model, Japanese Journal of Statistics and Data Science, in press. <https://doi.org/10.1007/s42081-025-00314-0>
-   Satoh, K. (2025) Applying non-negative matrix factorization with covariates to label matrix for classification. arXiv preprint arXiv:2510.10375. <https://arxiv.org/abs/2510.10375>


# References

-   Ding, C., Li, T., Peng, W. and Park, H. (2006) Orthogonal Nonnegative Matrix Tri-Factorizations for Clustering, Proceedings of the 12th ACM SIGKDD international conference on Knowledge discovery and data mining, 126-135. <https://doi.org/10.1145/1150402.1150420>
-   Potthoff, R.F., and Roy, S.N. (1964). A generalized multivariate analysis of variance model useful especially for growth curve problems. Biometrika, 51, 313-326. <https://doi.org/10.2307/2334137>




# Examples

Here are practical examples demonstrating the capabilities of **nmfkc**, ranging from simple matrix operations to complex time-series forecasting and classification tasks.

## Note that these scripts have been tested on Windows 11.

| ID | Description | Key Features |
|----|-------------|--------------|
| 0  | Simple matrix operations | Basic usage, Matrix decomposition |
| 1  | Longitudinal data (COVID-19) | Rank selection, Clustering on Map |
| 2  | Spatiotemporal Analysis (Weather) | **Kernel NMF**, Cross-validation for beta |
| 3  | Topic model (Inaugural address) | Text mining, Temporal topic evolution |
| 4  | OD data (Inter-prefecture flow) | **Sankey diagram**, Stability analysis |
| 5  | Kernel ridge regression (Motorcycle) | Non-linear regression, Interpolation |
| 6  | Growth curve model (Orthodont) | Handling dummy variables (Sex) |
| 7  | Autoregression (AirPassengers) | **NMF-VAR**, Forecasting, Stationarity check |
| 8  | Vector Autoregression (Canada) | Multivariate Time Series, **Granger Causality (DOT)** |
| 9 | Classification (Iris) | **NMF-LAB**, Supervised learning |

-----

## 0\. Simple matrix operations

This example demonstrates the basic concept of NMF: decomposing a matrix $Y$ into a basis matrix $X$ and a coefficient matrix $B$.

```r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")

# 1. Prepare synthetic data
# X: Basis matrix (2 bases), B: Coefficient matrix
(X <- cbind(c(1,0,1),c(0,1,0)))
(B <- cbind(c(1,0),c(0,1),c(1,1)))
(Y <- X %*% B) # Y is the observation matrix to be decomposed

# 2. Perform NMF
library(nmfkc)
# Decompose Y into X and B with Rank(Q) = 2
res <- nmfkc(Y, Q=2, epsilon=1e-6)

# 3. Check results
res$X # Estimated Basis
res$B # Estimated Coefficients

# 4. Diagnostics
plot(res)    # Convergence plot of the objective function
summary(res) # Summary statistics
```

## 1\. Longitudinal data: COVID-19 in Japan

An example of analyzing spatiotemporal data (prefecture x time). It includes **Rank Selection** to determine the optimal number of clusters and visualizes the results on a map.

  - Data Source: [https://www.mhlw.go.jp/stf/covid-19/open-data.html](https://www.mhlw.go.jp/stf/covid-19/open-data.html)

<!-- end list -->

```r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
# install.packages("NipponMap")

# 1. Data Preparation
d <- read.csv("https://covid19.mhlw.go.jp/public/opendata/newly_confirmed_cases_daily.csv")
n <- dim(d)
# Log-transform the infection counts (excluding 'Date' and 'ALL' columns)
Y <- log10(1 + d[, 2 + 1:47]) 
rownames(Y) <- unique(d$Date)

# 2. Rank Selection Diagnostics
library(nmfkc)
par(mfrow=c(1,1))
# Evaluate Ranks Q=2 to 8. 
# save.time=FALSE enables Element-wise Cross-Validation (Robust)
result.rank <- nmfkc.rank(Y, Q=2:8, save.time=FALSE)
round(result.rank$criteria, 3)

# 3. Perform NMF with Optimal Rank (e.g., Q=3)
# Q <- result.rank$rank.best
Q <- 3
result <- nmfkc(Y, Q=Q, epsilon=1e-5, prefix="Region")
plot(result, type="l", col=2) # Convergence

# 4. Visualization: Individual Fit
# Compare observed (black) vs fitted (red) for each prefecture
par(mfrow=c(7,7), mar=c(0,0,0,0)+0.1, cex=1)
for(n in 1:ncol(Y)){
  plot(Y[,n], axes=FALSE, type="l") 
  lines(result$XB[,n], col=2) 
  legend("topleft", legend=colnames(Y)[n], x.intersp=-0.5, bty="n")  
  box()
}

# 5. Visualization: Temporal Patterns (Basis X)
# Show the 3 extracted temporal patterns (Waves of infection)
par(mfrow=c(Q,1), mar=c(0,0,0,0), cex=1)
for(q in 1:Q){
  barplot(result$X[,q], col=q+1, border=q+1, las=3,
    ylim=range(result$X), ylab=paste0("topic ", q)) 
  legend("left", fill=q+1, legend=colnames(result$X)[q])
}

# 6. Visualization: Clustering on Map
# Show which prefecture belongs to which pattern
library(NipponMap)
par(mfrow=c(1,1), mar=c(5,4,4,2)+0.1)
jmap <- JapanPrefMap(col="white", axes=TRUE)

# Soft clustering (Pie charts on map)
stars(x=t(result$B.prob), scale=FALSE,
      locations=jmap, key.loc=c(145,34),
      draw.segments=TRUE, len=0.7, labels=NULL,
      col.segments=c(1:Q)+1, add=TRUE)

# Heatmap of coefficients
heatmap(t(result$B.prob))

# Hard clustering (Colored map)
library(NipponMap)
par(mfrow=c(1,1), mar=c(5,4,4,2)+0.1)
jmap <- JapanPrefMap(col=result$B.cluster+1, axes=TRUE)
```

## 2\. Spatiotemporal Analysis: CanadianWeather

This example compares standard NMF with **Kernel NMF**. By using location information as covariates (Kernel), we can predict coefficients for unobserved locations (Interpolation).

```r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
# install.packages("fda")

library(fda)
data(CanadianWeather)
d <- CanadianWeather$dailyAv[,,1] # Temperature data
Y <- d - min(d) # Ensure non-negative
u0 <- CanadianWeather$coordinates[,2:1] # Lat/Lon
u0[,1] <- -u0[,1]

# -----------------------------------
# Method A: Standard NMF (No Covariates)
# -----------------------------------
library(nmfkc)
result <- nmfkc(Y, Q=2, prefix="Trend")
plot(result)

# Visualization: Basis functions (Temperature patterns)
par(mfrow=c(1,1), mar=c(5,4,2,2)+0.1, cex=1)
plot(result$X[,1], type="n", ylim=range(result$X), ylab="basis function")
Q <- ncol(result$X)  
for(q in 1:Q) lines(result$X[,q], col=q+1)
legend("topright", legend=1:Q, fill=1:Q+1)

# Visualization: Spatial Distribution
par(mfrow=c(1,1), mar=c(5,4,2,2)+0.1, cex=1)
plot(u0, pch=19, cex=3, col=result$B.cluster+1, main="Hard Clustering")
text(u0, colnames(Y), pos=1)

# -----------------------------------
# Method B: NMF with Covariates (Linear)
# -----------------------------------
# Using raw coordinates as covariates A
U <- t(nmfkc.normalize(u0)) 
A <- rbind(1, U) # Add intercept
result <- nmfkc(Y, A, Q=2, prefix="Trend")
result$r.squared

# -----------------------------------
# Method C: NMF with Covariates (Gaussian Kernel)
# -----------------------------------
# A is constructed using a Kernel function of locations U.
# K(u,v) = exp{-beta * |u-v|^2}

# 1. Optimize Gaussian kernel parameter (beta) via Cross-Validation
result.beta <- nmfkc.kernel.beta.cv(Y, Q=2, U=U, beta=c(0.5, 1, 2, 5, 10))
(best.beta <- result.beta$beta)

# 2. Perform Kernel NMF
A <- nmfkc.kernel(U, beta=best.beta)
result <- nmfkc(Y, A, Q=2, prefix="Trend")
result$r.squared 

# 3. Prediction / Interpolation on a Mesh
# Predict coefficients B for new locations V
v <- seq(from=0, to=1, length=20)
V <- t(cbind(expand.grid(v,v))) # Grid points
A_new <- nmfkc.kernel(U, V, beta=best.beta)

# B_new = C * A_new
B <- result$C %*% A_new
B.prob <- prop.table(B, 2)

# Visualize estimated probability surface
q <- 1
z <- matrix(B.prob[q,], nrow=length(v)) 
par(mfrow=c(1,1), mar=c(5,4,2,2)+0.1, cex=1)
filled.contour(v, v, z, main=paste0("Predicted Probability of Pattern ",q),
  color.palette = function(n) hcl.colors(n, "Greens3", rev=TRUE),
  plot.axes={
     points(t(U), col=7, pch=19)
     text(t(U), colnames(U), pos=3)
  }
)
```

## 3\. Topic model: US Presidential Inaugural Addresses

Application to text mining. NMF extracts topics ($X$) and their weights ($B$). By using "Year" as a covariate, we can visualize how **topics evolve over time**.

```r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
# install.packages("quanteda")

# 1. Text Preprocessing using quanteda
library(quanteda)
corp <- corpus(data_corpus_inaugural)
tok <- tokens(corp)
tok <- tokens_remove(tok, pattern=stopwords("en", source="snowball"))
df <- dfm(tok)
df <- dfm_select(df, min_nchar=3)
df <- dfm_trim(df, min_termfreq=100)
d <- as.matrix(df)
# Sort by frequency
index <- order(colSums(d), decreasing=TRUE) 
d <- d[,index] 

# -----------------------------------
# Standard Topic Model (NMF)
# -----------------------------------
Y <- t(d) # Term-Document Matrix
Q <- 3
library(nmfkc)
result <- nmfkc(Y, Q=Q, prefix="Topic")

# Interpret Topics (High probability words)
Xp <- result$X.prob
for(q in 1:Q){
  message(paste0("----- Featured words on Topic [", q, "] -----"))
  print(paste0(rownames(Xp), "(", rowSums(Y), ") ", round(100*Xp[,q], 1), "%")[Xp[,q]>=0.5])
}

# Visualize Topic Proportions per President
par(mfrow=c(1,1), mar=c(10,4,4,2)+0.1, cex=1)
barplot(result$B.prob, col=1:Q+1, legend=TRUE, las=3,
  ylab="Probabilities of topics")

# -----------------------------------
# Temporal Topic Model (Kernel NMF)
# -----------------------------------
# Use 'Year' as covariate to smooth topic trends
U <- t(as.matrix(corp$Year))

# Optimize beta
result.beta <- nmfkc.kernel.beta.cv(Y, Q=3, U, beta=c(0.2, 0.5, 1, 2, 5)/10000)
(best.beta <- result.beta$beta)

# Perform Kernel NMF
A <- nmfkc.kernel(U, beta=best.beta)
result <- nmfkc(Y, A, Q, prefix="Topic")

# Visualize Smooth Topic Evolution
colnames(result$B.prob) <- corp$Year
par(mfrow=c(1,1), mar=c(5,4,2,2)+0.1, cex=1)
barplot(result$B.prob, col=1:Q+1, legend=TRUE, las=3, 
        ylab="Probability of topic (Smoothed)")
```

## 4\. Origin-Destination (OD) data

Analyzing flow data between prefectures. This example uses **Sankey diagrams** to visualize the stability of clustering results across different ranks ($Q$).

```r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
# install.packages("httr")
# install.packages("readxl")
# install.packages("alluvial")
# install.packages("NipponMap")
# install.packages("RColorBrewer")

# 1. Data Download & Formatting (Japanese OD data)
library(httr)
url <- "https://www.e-stat.go.jp/stat-search/file-download?statInfId=000040170612&fileKind=0"
GET(url, write_disk(tmp <- tempfile(fileext=".xlsx")))
library(readxl)
d <- as.data.frame(read_xlsx(tmp, sheet=1, skip=2))
# Filter for specific flow type
d <- d[d[,1]=="1" & d[,3]=="02" & d[,5]!="100" & d[,7]!="100", ]
pref <- unique(d[,6])

# Create OD Matrix Y (47x47)
Y <- matrix(NA, nrow=47, ncol=47)
colnames(Y) <- rownames(Y) <- pref
d[,5] <- as.numeric(d[,5]); d[,7] <- as.numeric(d[,7])
for(i in 1:47) for(j in 1:47) Y[i,j] <- d[which(d[,5]==i & d[,7]==j), 9]
Y <- log(1 + Y) # Log transformation

# 2. Rank Selection & Diagnostic Plot
library(nmfkc)
# Check AIC, BIC, RSS for Q=2 to 12
nmfkc.rank(Y, Q=2:12, save.time=FALSE)

# 3. NMF Analysis (Q=7)
Q0 <- 7
res <- nmfkc(Y, Q=Q0, save.time=FALSE, prefix="Region")

# 4. Visualization: Silhouette Plot
si <- res$criterion$silhouette
barplot(si$silhouette, horiz=TRUE, las=1, col=si$cluster+1, 
        cex.names=0.5, xlab="Silhouette width")
abline(v=si$silhouette.mean, lty=3)

# 5. Visualization: Sankey Diagram (Stability Analysis)
# Visualize how clusters change as Q increases from 6 to 8
Q <- 6:8
cluster <- NULL
for(i in 1:length(Q)){
  res_tmp <- nmfkc(Y, Q=Q[i])
  cluster <- cbind(cluster, res_tmp$B.cluster)
}
library(alluvial)
alluvial(cluster, freq=1, axis_labels=paste0("Q=", Q), cex=2,
         col=cluster[,2]+1, border=cluster[,2]+1)
title("Cluster evolution (Sankey Diagram)")

# 6. Visualization: Map
library(NipponMap); library(RColorBrewer)
mypalette <- brewer.pal(12, "Paired")
par(mfrow=c(1,1), mar=c(5,4,4,2)+0.1)
jmap <- JapanPrefMap(col=mypalette[res$B.cluster], axes=TRUE)
text(jmap, pref, cex=0.5)  
title(main="OD Clustering Result")
```

## 5\. Kernel ridge regression

Demonstrates that **nmfkc** can be used for non-linear regression (smoothing) by treating 1D data as a matrix row and using Gaussian kernels.

```r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")

library(MASS)
d <- mcycle
x <- d$times
y <- d$accel
# Treat Y as a 1 x N matrix
Y <- t(as.matrix(y - min(y)))
U <- t(as.matrix(x))

# 1. Linear Regression (NMF with Linear Covariate)
A <- rbind(1, U) # Intercept + Linear term
library(nmfkc)
result <- nmfkc(Y, A, Q=1)
par(mfrow=c(1,1), mar=c(5,4,2,2)+0.1, cex=1)
plot(U, Y, main="Linear Fit")
lines(U, result$XB, col=2)

# 2. Kernel Regression (Gaussian Kernel)
# Optimize beta via Cross-Validation
result.beta <- nmfkc.kernel.beta.cv(Y, Q=1, U, beta=2:5/100)
(beta.best <- result.beta$beta)  

# Create Kernel Matrix
A <- nmfkc.kernel(U, beta=beta.best)
result <- nmfkc(Y, A, Q=1)

# Plot Fitted Curve (Smoothing)
plot(U, Y, main="Kernel Smoothing")
lines(U, as.vector(result$XB), col=2, lwd=2)

# Prediction on new data points
V <- matrix(seq(from=min(U), to=max(U), length=100), ncol=100)
A_new <- nmfkc.kernel(U, V, beta=beta.best)
XB_new <- predict(result, newA=A_new)
lines(V, as.vector(XB_new), col=4, lwd=2)
```

## 6\. Growth curve model

Analysis of repeated measures (Growth Curve). Shows how to incorporate categorical covariates (Sex) into the $A$ matrix to analyze group differences.

```r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
# install.packages("nlme")

library(nlme)
d <- Orthodont
t <- unique(d$age)
# Matrix: Age x Subject
Y <- matrix(d$distance, nrow=length(t))
colnames(Y) <- unique(d$Subject)
rownames(Y) <- t

# 1. Construct Covariate Matrix A
# Include Intercept and Dummy variable for Male
Male <- 1 * (d$Sex == "Male")[d$age == 8]
A <- rbind(rep(1, ncol(Y)), Male)
rownames(A) <- c("Const", "Male")

# 2. Perform NMF with Covariates
library(nmfkc)
result <- nmfkc(Y, A, Q=2, epsilon=1e-8)

# 3. Check Coefficients
# B = C * A. We can see the coefficients for Male vs Female.
(A0 <- t(unique(t(A)))) # Unique patterns (Female, Male)
B <- result$C %*% A0

# 4. Visualization
plot(t, Y[,1], ylim=range(Y), type="n", ylab="Distance", xlab="Age")
mycol <- ifelse(Male==1, 4, 2)
for(n in 1:ncol(Y)){
  lines(t, Y[,n], col=mycol[n]) # Individual trajectories
}
XB <- result$X %*% B
# Mean trajectories
lines(t, XB[,1], col=4, lwd=5) # Male
lines(t, XB[,2], col=2, lwd=5) # Female
legend("topleft", title="Sex", legend=c("Male", "Female"), fill=c(4,2))
```

## 7\. Autoregression: AirPassengers (NMF-VAR)

**nmfkc** can perform Vector Autoregression (NMF-VAR). This example includes **Lag Order Selection** and **Future Forecasting**.

```r
# ===================================================
# NMF-VAR forecasting example using AirPassengers data
# ---------------------------------------------------
# Required packages:
#   remotes::install_github("ksatohds/nmfkc")
#   install.packages("DOT")
# ===================================================

library(nmfkc)
library(DOT)

# ---- 1) Data Preparation ----
d_air <- AirPassengers

# Generate a time vector for labeling
time_air <- time(ts(1:length(d_air), start = c(1949, 1), frequency = 12))
time_vec_air <- round(as.vector(t(time_air)), 2)

# Convert the AirPassengers series to a 1×N log-transformed matrix
Y0_air <- log10(matrix(as.vector(d_air), nrow = 1))
colnames(Y0_air) <- time_vec_air
rownames(Y0_air) <- "Passengers"

# ---- 2) Lag Order Selection (Cross-Validation) ----
# Evaluate lag degrees from 1 to 25
cv_res <- nmfkc.ar.degree.cv(Y0_air, Q = 1, degree = 1:25)

# Print the optimal lag degree
cat("Optimal lag degree:", cv_res$degree, "\n")

# For this demo, use D = 12 (monthly seasonality)
D <- 12

# ---- 3) Construct AR(12) Matrices ----
a_air <- nmfkc.ar(Y0_air, degree = D, intercept = TRUE)
Y_air <- a_air$Y
A_air <- a_air$A

# ---- 4) Fit NMF-AR Model ----
# Rank = 1 for a univariate series
res_air <- nmfkc(Y = Y_air, A = A_air,
                 rank = 1, epsilon = 1e-9, maxit = 500000)

# Display model performance
cat("R-squared:", round(res_air$r.squared, 5), "\n")

# Check dynamic stability (spectral radius < 1 indicates stationarity)
nmfkc.ar.stationarity(res_air)

# ---- 5) Forecast Generation ----
h <- 24  # Forecast horizon (24 months ahead)
pred_res <- nmfkc.ar.predict(x = res_air, Y = Y0_air, n.ahead = h)
pred_val <- 10^as.vector(pred_res$pred)  # Convert from log10 scale back to original

# ---- 6) Generate Future Time Points ----
# Retrieve time settings from the original time series
total_len <- length(d_air)
start_time <- start(d_air)
freq <- frequency(d_air)

# Compute start time for forecasts (immediately after the last observation)
pred_start_time <- start_time[1] + (start_time[2] + total_len) / freq
pred_time <- seq(from = pred_start_time, by = 1 / freq, length.out = h)

# ---- 7) Combine Observed and Forecast Series ----
# Last observed point
last_obs_time <- tail(as.numeric(colnames(Y0_air)), 1)
last_obs_val  <- tail(10^as.vector(Y0_air), 1)

# Concatenate the end of observed data with the forecast start point
plot_pred_time <- c(last_obs_time, pred_time)
plot_pred_val  <- c(last_obs_val, pred_val)

# Define plot axis ranges
xlim_range <- range(c(as.numeric(colnames(Y0_air)), pred_time))
ylim_range <- range(c(10^as.vector(Y0_air), pred_val))

# ---- 8) Visualization ----
# Plot observed, fitted, and forecasted values
plot(as.numeric(colnames(Y0_air)), 10^(as.vector(Y0_air)),
     type = "l", col = "black", lwd = 1,
     xlim = xlim_range, ylim = ylim_range,
     xlab = "Year", ylab = "Air Passengers",
     main = "NMF-VAR Forecast (h = 24)")

# Add fitted values (red)
lines(as.numeric(colnames(Y_air)), 10^as.vector(res_air$XB),
      col = "red", lwd = 2)

# Add forecasts (blue dashed)
lines(plot_pred_time, plot_pred_val, col = "blue", lwd = 2, lty = 2)

# Add legend
legend("topleft",
       legend = c("Observed", "Fitted", "Forecast"),
       col = c("black", "red", "blue"),
       lty = c(1, 1, 2), lwd = 2)
```

## 8\. Vector Autoregression: Canada (Multivariate)

Multivariate NMF-VAR example. It also demonstrates how to visualize Granger causality using **DOT graphs**.

```r
# ===============================================
# NMF with VAR structure (nmfkc) × Canada dataset (vars) demo
# -----------------------------------------------
# Required packages:
#   remotes::install_github("ksatohds/nmfkc")
#   install.packages(c("vars", "DOT"))
# ===============================================

# ---- 0) Packages ----
suppressPackageStartupMessages({
  library(nmfkc)
  library(vars)
  library(DOT)   # dot() renders the Graphviz output (viewable in RStudio or browser)
})

# ---- 1) Data Preparation ----
# Canada: Quarterly ts data from 'vars' (1980Q1–2000Q4, 4 variables)
d0 <- Canada                # ts: multivariate time series (time × variables)
stopifnot(is.ts(d0))

# First difference (to model growth or change rates)
dd <- apply(d0, 2, diff)    # matrix (T-1) × P
# Normalize each column to [0,1] (for stability)
dn <- nmfkc.normalize(dd)   # (T-1) × P
# nmfkc expects variables in rows, samples in columns
Y0 <- t(dn)                 # P × N

# ---- 2) Build matrices for the VAR model (nmfkc.ar) ----
Q <- 2       # Rank (number of latent bases)
D <- 1       # Lag order (VAR(1))
ar_set <- nmfkc.ar(Y0, degree = D, intercept = TRUE)
Y  <- ar_set$Y   # P × (N-D)
A  <- ar_set$A   # (P*D + 1) × (N-D)

# Quick consistency check
if (ncol(Y) != ncol(A)) stop("Y and A must have the same number of columns (samples).")

# ---- 3) Fit NMF model (nmfkc) ----
# Objective: Euclidean distance (default). Use method="KL" for Kullback–Leibler divergence.
res <- nmfkc(Y = Y, A = A, Q = Q, prefix = "Condition", epsilon = 1e-6)
# Key metrics
cat("\n--- Fit summary ---\n",
    "r.squared  :", round(res$r.squared, 4), "\n",
    "method     :", res$method, "\n",
    "rank (Q)   :", res$rank, "\n",
    "runtime    :", res$runtime, "\n", sep="")

# Check dynamic stability (VAR stationarity)
st <- nmfkc.ar.stationarity(res)
cat("\n--- Stationarity (VAR) ---\n",
    "spectral radius :", round(st$spectral.radius, 4), "\n",
    "stationary      :", st$stationary, "\n", sep="")

# ---- 4) Causality visualization (DOT graph) ----
#   X (latent) → Y (observed), T−k → X connections grouped by lag
script <- nmfkc.ar.DOT(res, degree = D, intercept = TRUE, digits = 2, rankdir = "RL")
# Render to screen (display method depends on your environment)
dot(script)

# Optionally save DOT script to file:
# writeLines(script, "nmf_var_graph.dot")

# ---- 5) Transform fitted differences back to levels (denormalization + integration) ----
# res$XB is the fitted values of Y (P × (N-D)), still in diff-normalized scale
xb_diff_norm <- res$XB           # P × (N-D)
stopifnot(all(dim(xb_diff_norm) == dim(t(dn)[, -(1:D), drop=FALSE])))

# Denormalize using original differenced data (dd)
xb_diff <- nmfkc.denormalize(t(xb_diff_norm), dd)  # (N-D) × P (back to diff scale)

# Integrate differences to reconstruct level series:
# rbind(initial level, fitted differences) → cumulative sum
# Drop the first row (the initial point itself)
d1_mat <- apply(rbind(d0[1 + D, ], xb_diff), 2, cumsum)
d1     <- d1_mat[-1, , drop = FALSE]  # (N-D) × P

# ---- 6) Align time axis using ts objects ----
# Convert d1 to a ts with the same frequency as d0, starting D periods later
d1_ts <- ts(d1,
            start     = start(d0) + c(0, D) / frequency(d0),
            frequency = frequency(d0))

# ---- 7) Plot original vs fitted series ----
par(mfrow = c(2, 2))
for (p in 1:nrow(Y0)) {
  plot(d0[, p], col = "gray",
       main = colnames(d0)[p], ylab = colnames(d0)[p], xlab = "")
  lines(d1_ts[, p], col = 2, lwd = 2)
  legend("topleft", bty = "n",
         legend = c("Original", "Fitted (reconstructed)"),
         col = c("gray", "red"), lwd = c(1, 2))
}

# ---- 8) Additional diagnostics (optional) ----
# Convergence plot of the optimization process
par(mfrow = c(1, 1))
plot(res, main = paste0("r.squared = ", round(res$r.squared, 3)))

# Residual heatmaps (Y, XB, E)
# nmfkc.residual.plot(Y = Y, result = res)
```

## 9\. Classification: Iris (NMF-LAB)

**NMF-LAB** (Label-based NMF) example. By treating class labels as the target matrix $Y$ (one-hot encoding) and features as covariates $U$ (kernelized to $A$), NMF can be used for **Supervised Learning**.

```r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
library(nmfkc)

# 1. Data Preparation
label <- iris$Species
# Convert labels to One-Hot Matrix Y
Y <- nmfkc.class(label) 
# Features Matrix U (Normalized)
U <- t(nmfkc.normalize(iris[,-5]))

# 2. Kernel Parameter Optimization
# Heuristic estimation of beta (Median distance)
res.beta <- nmfkc.kernel.beta.nearest.med(U)
# Cross-validation for fine-tuning
res.cv <- nmfkc.kernel.beta.cv(Y, Q=length(unique(label)), U, beta=res.beta$beta_candidates)
best.beta <- res.cv$beta

# 3. Fit NMF-LAB Model
A <- nmfkc.kernel(U, beta=best.beta)
res <- nmfkc(Y=Y, A=A, Q=length(unique(label)), prefix="Class")

# 4. Prediction and Evaluation
# Predict class based on highest probability in XB
fitted.label <- predict(res, type="class")
(f <- table(fitted.label, label))
message("Accuracy: ", round(100 * sum(diag(f)) / sum(f), 2), "%")
```




# Author

-   Kenichi Satoh, [homepage](https://sites.google.com/view/ksatoh/english)
