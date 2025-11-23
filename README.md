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
| 7  | Binary repeated measures (Clinical) | Categorical data analysis |
| 8  | Autoregression (AirPassengers) | **NMF-VAR**, Forecasting, Stationarity check |
| 9  | Vector Autoregression (Canada) | Multivariate Time Series, **Granger Causality (DOT)** |
| 10 | Classification (Iris) | **NMF-LAB**, Supervised learning |
| 11 | Classification (Penguins) | Handling NA, Kernel classification |

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
round(result.rank, 3)

# Plot ICp (Information Criterion)
plot(result.rank$Q, result.rank$ICp, col=2, type="l", main="Model Selection (ICp)")
text(result.rank$Q, result.rank$ICp, result.rank$Q)

# 3. Perform NMF with Optimal Rank (e.g., Q=3)
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
plot(result$X[,1], type="n", ylim=range(result$X[,1]), ylab="basis function")
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
z <- matrix(B.prob[2,], nrow=length(v)) 
par(mfrow=c(1,1), mar=c(5,4,2,2)+0.1, cex=1)
filled.contour(v, v, z, main="Predicted Probability of Pattern 2",
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

## 8\. Autoregression: AirPassengers (NMF-VAR)

**nmfkc** can perform Vector Autoregression (NMF-VAR). This example includes **Lag Order Selection** and **Future Forecasting**.

```r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
# install.packages("DOT")
library(nmfkc)

# 1. Data Preparation
data(AirPassengers)
d <- AirPassengers
time.vec <- round(as.vector(time(d)), 2)
Y0 <- log10(matrix(as.vector(d), nrow=1)) # 1 x N matrix
colnames(Y0) <- time.vec

# 2. Model Selection (Lag Order)
# Find optimal lag D using Cross-Validation
cv_res <- nmfkc.ar.degree.cv(Y0, Q=1, degree=1:16)
D <- cv_res$degree

# 3. Model Fitting (NMF-VAR)
# Create VAR matrices Y (target) and A (lagged covariates)
ar_data <- nmfkc.ar(Y0, degree=D, intercept=TRUE)
res <- nmfkc(Y=ar_data$Y, A=ar_data$A, rank=1, epsilon=1e-9)

# Stationarity Check (Spectral Radius)
print(nmfkc.ar.stationarity(res))

# 4. Forecasting (24 steps ahead)
pred_res <- nmfkc.ar.predict(x=res, Y=Y0, n.ahead=24)
pred_val <- 10^as.vector(pred_res$pred)

# 5. Visualization
plot(as.numeric(colnames(Y0)), 10^as.vector(Y0), type="l", xlim=c(1949, 1963), 
     ylab="Passengers", main="NMF-VAR Forecast")
lines(as.numeric(colnames(ar_data$Y)), 10^as.vector(res$XB), col="red", lwd=2) # Fitted
lines(pred_res$time, pred_val, col="blue", lwd=2) # Forecast
legend("topleft", legend=c("Observed","Fitted","Forecast"), col=c(1,2,4), lty=1, lwd=2)
```

## 9\. Vector Autoregression: Canada (Multivariate)

Multivariate NMF-VAR example. It also demonstrates how to visualize Granger causality using **DOT graphs**.

```r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
# install.packages("vars")
# install.packages("DOT")
library(nmfkc)
library(vars)

# 1. Data Preparation (Canada Economic Data)
d0 <- Canada
dd <- apply(d0, 2, diff) # First difference
dn <- nmfkc.normalize(dd) # Normalize to [0,1]
Y0 <- t(dn) # Variables x Time

# 2. Create VAR Matrices
Q <- 2; D <- 1
a <- nmfkc.ar(Y0, degree=D, intercept=TRUE)
res <- nmfkc(Y=a$Y, A=a$A, Q=Q, prefix="Condition", epsilon=1e-6)

# 3. Causality Visualization (DOT Graph)
# Generates a script for Graphviz/DOT to visualize variable relationships
script <- nmfkc.ar.DOT(res, intercept=TRUE, digits=2)
library(DOT)
dot(script) # Renders the causal graph

# 4. Fitted Curve (with Denormalization)
t1 <- t(res$XB)
dn1 <- nmfkc.denormalize(t1, dd)
# Reconstruct from differences
d1 <- apply(rbind(d0[1+D,], dn1), 2, cumsum) 

# Plot Original vs Fitted
par(mfrow=c(2,2))
for(p in 1:nrow(Y0)){
  plot(d0[,p], type="l", col="gray", ylab=rownames(Y0)[p], main=rownames(Y0)[p])
  lines(time(d0)[-(1:(D+1))], d1[,p], col=2, lwd=2)
}
```

## 10\. Classification: Iris (NMF-LAB)

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







# Examples

## Note that these scripts have been tested on windows 11.

 0.  Simple matrix operations
 1.  Longitudinal data: COVID-19 in Japan
 2.  Spatiotemporal Analysis: CanadianWeather
 3.  Topic model: data_corpus_inaugural
 4.  Origin-Destination (OD) data: Japanese Inter-prefecture flow
 5.  Kernel ridge regression: mcycle
 6.  Growth curve model: Orthodont
 7.  Binary repeated measures: Table 6, Koch et al.(1977)
 8.  Autoregression: AirPassengers
 9.  Vector Autoregression: Canada
10.  Classification: Iris 
11.  Classification: penguins

## 0.  Simple matrix operations

``` r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
(X <- cbind(c(1,0,1),c(0,1,0)))
(B <- cbind(c(1,0),c(0,1),c(1,1)))
(Y <- X %*% B)

library(nmfkc)
res <- nmfkc(Y,Q=2,epsilon=1e-6)
res$X
res$B

plot(res)
summary(res)
``` 

## 1. Longitudinal data

-   COVID-19 in Japan
-   https://www.mhlw.go.jp/stf/covid-19/open-data.html

``` r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
# install.packages("NipponMap")

d <- read.csv("https://covid19.mhlw.go.jp/public/opendata/newly_confirmed_cases_daily.csv")
n <- dim(d)
Y <- log10(1+d[,2+1:47])
rownames(Y) <- unique(d$Date)

# rank selection diagnostics
library(nmfkc)
par(mfrow=c(1,1))
result.rank <- nmfkc.rank(Y,Q=2:8,save.time=F)
round(result.rank,3)

# ICp
plot(result.rank$Q,result.rank$ICp,col=2,type="l")
text(result.rank$Q,result.rank$ICp,result.rank$Q)

# nmf
Q <- 3
result <- nmfkc(Y,Q=Q,epsilon=1e-5,prefix="Region")
plot(result,type="l",col=2)

# individual fit
par(mfrow=c(7,7),mar=c(0,0,0,0)+0.1,cex=1)
for(n in 1:ncol(Y)){
  plot(Y[,n],axes=F,type="l") # observation
  lines(result$XB[,n],col=2) # fitted values
  legend("topleft",legend=colnames(Y)[n],x.intersp=-0.5,bty="n")  
  box()
}

# basis function of which sum is 1
par(mfrow=c(Q,1),mar=c(0,0,0,0),cex=1)
for(q in 1:Q){
  barplot(result$X[,q],col=q+1,border=q+1,las=3,
    ylim=range(result$X),ylab=paste0("topic ",q)) 
  legend("left",fill=q+1,legend=colnames(result$X)[q])
}

# cluster membership probability based on coefficients
n <- 1
result$B[,n]
result$B.prob[,n]
result$B.cluster[n]

# soft clustering based on B.prob
library(NipponMap)
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
jmap <- JapanPrefMap(col="white",axes=TRUE)
stars(x=t(result$B.prob),scale=F,
      locations=jmap,key.loc =c(145,34),
      draw.segments=T,len=0.7,labels=NULL,
      col.segments=c(1:Q)+1,add=T)

# heatmap
heatmap(t(result$B.prob))

# hard clustering based on B.prob
library(NipponMap)
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
jmap <- JapanPrefMap(col=result$B.cluster+1,axes=TRUE)
```

## 2. Spatiotemporal Analysis

-   CanadianWeather

``` r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
# install.packages("fda")

library(fda)
data(CanadianWeather)
d <- CanadianWeather$dailyAv[,,1]
Y <- d-min(d)
u0 <- CanadianWeather$coordinates[,2:1]
u0[,1] <- -u0[,1]

#------------------
# without covariate
#------------------
library(nmfkc)
result <- nmfkc(Y,Q=2,prefix="Trend")
plot(result)

# individual fit
par(mfrow=c(6,6),mar=c(0,0,0,0)+0.1,cex=1)
for(n in 1:ncol(Y)){
  plot(Y[,n],ylim=range(Y),axes=F,type="l") # observation
  lines(result$XB[,n],col=2) # fitted values
  abline(h=-min(d),col="gray",lty=3,lwd=3)
  legend("topleft",legend=colnames(Y)[n],x.intersp=-0.5,bty="n")  
  box()
}

# basis function of which sum is 1
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(result$X[,1],type="n",ylim=range(result$X[,1]),
  ylab="basis function")
Q <- ncol(result$X)  
for(q in 1:Q) lines(result$X[,q],col=q+1)
legend("topright",legend=1:Q,fill=1:Q+1)

# soft clustering based on B.prob
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(u0,type="n")
legend("topright",legend=1:Q,fill=1:Q+1)
stars(t(result$B.prob),
      locations=u0,scale=F,
      draw.segments=TRUE,labels=colnames(Y),
      col.segments=1:Q+1,
      len=max(u0)/30,add=T)

# hard clustering based on B.cluster
plot(u0,pch=19,cex=3,col=result$B.cluster+1)
text(u0,colnames(Y),pos=1)
legend("topright",legend=c("1","2"),fill=c(2,3))

#------------------
# with covariates using location information
#------------------
U <- t(nmfkc.normalize(u0)) # normarization of covariates
A <- rbind(1,U)
result <- nmfkc(Y,A,Q=2,prefix="Trend")
result$r.squared

#------------------
# with covariates using kernel matrix A
#------------------
# U=[u1,...,uN]
# A= |K(u1,u1),...,K(u1,uN)|
#    |                     |
#    |K(uN,u1),...,K(uN,uN)|
# K(u,v)=exp{-beta*|u-v|^2}

# k-fold cross validation for beta
result.beta <- nmfkc.kernel.beta.cv(Y,Q=2,U,beta=c(0.5,1,2,5,10))
(best.beta <- result.beta$beta)

# create kernel with best beta
A <- nmfkc.kernel(U,beta=best.beta)
result <- nmfkc(Y,A,Q=2,prefix="Trend")
result$r.squared # less than nmf without covariates

# prediction of coefficients(b) on mesh point V
# U=[u1,...,uN]
# V=[v1,...,vM]
# A= |K(u1,v1),...,K(u1,vM)|
#    |                     |
#    |K(uN,v1),...,K(uN,vM)|
v <- seq(from=0,to=1,length=20)
V <- t(cbind(expand.grid(v,v)))
plot(t(V))
A <- nmfkc.kernel(U,V,beta=best.beta)
B <- result$C %*% A
B.prob <- prop.table(B,2)
B.prob[,1:6]

# soft clustering based on B.prob (basis function 2) by using covariates
z <- matrix(B.prob[2,],nrow=length(v)) 
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
filled.contour(v,v,z,main="probability of basis function 2",
  color.palette = function(n) hcl.colors(n,"Greens3",rev=TRUE),
               plot.axes={
                 points(t(U),col=7,pch=19)
                 text(t(U),colnames(U),pos=3)
               }
)
```

## 3. Topic model

-   data_corpus_inaugural
-   US presidential inaugural address texts

``` r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
# install.packages("quanteda")

#------------------
# text analysis
#------------------
library(quanteda)
corp <- corpus(data_corpus_inaugural)
head(corp,3)
tail(corp,3)
tok <- tokens(corp)
tok <- tokens_remove(tok,pattern=stopwords("en",source="snowball"))
df <- dfm(tok)
df <- dfm_select(df,min_nchar=3)
df <- dfm_trim(df,min_termfreq=100)
d <- as.matrix(df)
index <- order(colSums(d),decreasing=T) 
d <- d[,index] # document-word matrix
paste0(colnames(d)[1:30],"(",colSums(d)[1:30],")") # Top 30 most frequent words

#------------------
# without covariates
#------------------
Y <- t(d)
Y[1:20,c(1,ncol(Y))]
Q <- 3
library(nmfkc)
result <- nmfkc(Y,Q=Q,prefix="Topic")
result$r.squared # coefficient of determination

# soft clustering
par(mfrow=c(1,1),mar=c(10,4,4,2)+0.1,cex=1)
barplot(result$B.prob,col=1:Q+1,legend=T,las=3,
  ylab="Probabilities of topics")

# hard clustering
par(mfrow=c(1,1),mar=c(10,4,4,2)+0.1,cex=1)
cluster.prob <- 0*1:ncol(Y)
for(n in 1:ncol(Y)) cluster.prob[n] <- result$B.prob[result$B.cluster[n],n]
barplot(cluster.prob,col=result$B.cluster+1,names=colnames(Y),legend=T,
  las=3,ylab="Probability",args.legend=list(bg="white"))

# basis function of which sum is 1
par(mfrow=c(Q,1),mar=c(8,4,1,1),cex=0.6)
for(q in 1:Q){
  barplot(result$X[,q],col=q+1,border=q+1,las=3,
    ylim=range(result$X),ylab=paste0("topic ",q)) 
}

# contribution of words to each topics
Xp <- result$X.prob
par(mfrow=c(1,1),mar=c(8,4,1,1)+0.1,cex=0.6)
barplot(t(Xp),las=3,col=1:Q+1,
  ylab="Proportion of words on each topic")
legend("topleft",fill=1:Q+1,legend=paste0("topic ",1:Q))
for(q in 1:Q){
  message(paste0("----- featured words on topic [",q,"] -----"))
  print(paste0(rownames(Xp),"(",rowSums(Y),")",round(100*Xp[,q],1),"%")[Xp[,q]>=0.5])
}

#------------------
# with covariates using covariate matrix U
#------------------
U <- t(as.matrix(corp$Year))
result.beta <- nmfkc.kernel.beta.cv(Y,Q=3,U,beta=c(0.2,0.5,1,2,5)/10000)
(best.beta <- result.beta$beta)

# create kernel with best beta
A <- nmfkc.kernel(U,beta=best.beta)
result <- nmfkc(Y,A,Q,prefix="Topic")
result$r.squared # less than nmf without covariates

# Topic probability changing over time
colnames(result$B.prob) <- corp$Year
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
barplot(result$B.prob,col=1:Q+1,legend=T,las=3,ylab="Probability of topic")
```

## 4. Origin-Destination (OD) data: Japanese Inter-prefecture flow

-   e-stat run by Japanese Government Statistics

``` r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
# install.packages("httr")
# install.packages("readxl")
# install.packages("alluvial")
# install.packages("NipponMap")
# install.packages("RColorBrewer")

# download data and its formatting
library(httr)
# https://www.e-stat.go.jp/stat-search/files?stat_infid=000040170612
url <- "https://www.e-stat.go.jp/stat-search/file-download?statInfId=000040170612&fileKind=0"
GET(url, write_disk(tmp <- tempfile(fileext=".xlsx")))
library(readxl)
d <- as.data.frame(read_xlsx(tmp,sheet=1,skip=2))
colnames(d)
d <- d[d[,1]=="1",]
d <- d[d[,3]=="02",]
d <- d[d[,5]!="100",]
d <- d[d[,7]!="100",]
pref <- unique(d[,6])

#------------------
# without covariates
#------------------
Y <- matrix(NA,nrow=47,ncol=47)
colnames(Y) <- pref
rownames(Y) <- pref
d[,5] <- as.numeric(d[,5])
d[,7] <- as.numeric(d[,7])
for(i in 1:47)for(j in 1:47){
  Y[i,j] <- d[which(d[,5]==i&d[,7]==j),9]
}
Y <- log(1+Y)
Y[1:6,1:6]

# rank selection diagnostics
library(nmfkc)
nmfkc.rank(Y,Q=2:12,save.time=F)

# nmf
Q0 <- 7
res <- nmfkc(Y,Q=Q0,save.time=F,prefix="Region")
plot(res$objfunc.iter,type="o",
     main=paste0("Q=",Q0,", R^2=",round(res$r.squared,3)))
     
# Silhouette: intra-cluster accumulation and inter-cluster discrepancy
si <- res$criterion$silhouette
barplot(si$silhouette,horiz=T,las=1,col=si$cluster+1,cex.names=0.5,xlab="")
abline(v=si$silhouette.mean,lty=3)
legend("bottomleft",fill=1:Q0,legend=1:Q0)

# Sankey diagram: stability of hard clustering by 
Q <- 6:8
cluster <- NULL
for(i in 1:length(Q)){
  res <- nmfkc(Y,Q=Q[i])
  cluster <- cbind(cluster,res$B.cluster)
}
library(alluvial)
alluvial(cluster,freq=1,axis_labels=paste0("Q=",Q),cex=2,
         col=cluster[,2]+1,border=cluster[,2]+1)

# basis function of which sum is 1
library(NipponMap)
library(RColorBrewer)
mypalette <- brewer.pal(9,"YlOrRd")
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1,cex=0.6)
for(j in 1:Q0){
  cutp <- as.numeric(
    cut(res$B.prob[j,],
        breaks=seq(from=0,to=1,length=10),
        include.lowest=T))
  mycol <- mypalette[cutp]
  par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
  jmap <- JapanPrefMap(col=mycol,axes=TRUE,
                       main=paste0("basis[",j,"]"))
  text(jmap,pref,cex=0.5)  
}

# soft clustering based on B.prob
mypalette <- brewer.pal(12,"Paired")
tp <- t(res$B.prob)
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
jmap <- JapanPrefMap(col="white",axes=TRUE)
stars(x=tp,scale=F,locations=jmap,key.loc =c(145,34),
      draw.segments=TRUE,len=1,labels=NULL,
      col.segments=mypalette[1:Q0],add=T)
title(main="Inter-prefecture flow: weekdays - operations")

# hard clustering based on B.cluster
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
jmap <- JapanPrefMap(col=mypalette[res$B.cluster],axes=TRUE)
text(jmap,pref,cex=0.5)  
legend("topleft",fill=mypalette[1:Q0],
       legend=1:Q0,title="basis")
title(main="Inter-prefecture flow: weekdays - operations")
```

## 5. Kernel ridge regression

-   mcycle

``` r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")

library(MASS)
d <- mcycle
x <- d$times
y <- d$accel
Y <- t(as.matrix(y-min(y)))
U <- t(as.matrix(x))
# scatter plot
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(U,Y)

# linear curve
A <- rbind(1,U)
library(nmfkc)
result <- nmfkc(Y,A,Q=1)
result$r.squared
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(U,Y)
lines(U,result$XB,col=2)

# cv for optimization of beta
result.beta <- nmfkc.kernel.beta.cv(Y,Q=1,U,beta=2:5/100)
(beta.best <- result.beta$beta)  
A <- nmfkc.kernel(U,beta=beta.best)
result <- nmfkc(Y,A,Q=1)
result$r.squared

# fitted curve
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(U,Y)
lines(U,as.vector(result$XB),col=2,lwd=2)
# fitted curve for new data
V <- matrix(seq(from=min(U),to=max(U),length=100),ncol=100)
A <- nmfkc.kernel(U,V,beta=beta.best)
XB <- predict(result,newA=A)
lines(V,as.vector(XB),col=4,lwd=2)
```

## 6. Growth curve model

-   Orthodont

``` r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
# install.packages("nlme")

library(nlme)
d <- Orthodont
head(d)
t <- unique(d$age)
Y <- matrix(d$distance,nrow=length(t))
colnames(Y) <- unique(d$Subject)
rownames(Y) <- t

#------------------
# with covariates
#------------------
Male <- 1*(d$Sex=="Male")[d$age==8]
table(Male)
A <- rbind(rep(1,ncol(Y)),Male)
rownames(A) <- c("Const","Male")

# nmf with covariates
library(nmfkc)
result <- nmfkc(Y,A,Q=2,epsilon=1e-8)
result$r.squared

# basis matrix
print.table(round(result$X,2),zero.print="")

# parameter matrix
print.table(round(result$C,2),zero.print="")

# unique covariates and coefficient matrix
(A0 <- t(unique(t(A))))
B <- result$C %*% A0

# individual fit
plot(t,Y[,1],ylim=range(Y),type="n")
mycol <- ifelse(Male==1,4,2)
for(n in 1:ncol(Y)){
  lines(t,Y[,n],col=mycol[n])
}
XB <- result$X %*% B
lines(t,XB[,1],col=4,lwd=5)
lines(t,XB[,2],col=2,lwd=5)
legend("topleft",title="Sex",legend=c("Male","Female"),fill=c(4,2))
```

## 7. Binary repeated measures

-   Table 6, Koch et al.(1977) Biometrics, 33(1), 133–158. <https://doi.org/10.2307/2529309>
-   Repeated categorical outcome analysis, <https://wnarifin.github.io/medstat.html>

``` r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")

Id <- rep(1:340,each=3)
Mild <- rep(c(1,0),c(150*3,190*3))
NewDrug <- rep(c(0,1,0,1),c(80*3,70*3,100*3,90*3))
Week <- rep(c(1,2,4),times=340)
Normal <- unlist(c(
  rep(list(c(1,1,1),c(1,1,0),c(1,0,1),c(1,0,0),c(0,1,1),c(0,1,0),c(0,0,1),c(0,0,0)),c(16,13,9,3,14,4,15,6)),
  rep(list(c(1,1,1),c(1,1,0),c(1,0,1),c(1,0,0),c(0,1,1),c(0,1,0),c(0,0,1),c(0,0,0)),c(31,0,6,0,22,2,9,0)),
  rep(list(c(1,1,1),c(1,1,0),c(1,0,1),c(1,0,0),c(0,1,1),c(0,1,0),c(0,0,1),c(0,0,0)),c(2,2,8,9,9,15,27,28)),
  rep(list(c(1,1,1),c(1,1,0),c(1,0,1),c(1,0,0),c(0,1,1),c(0,1,0),c(0,0,1),c(0,0,0)),c(7,2,5,2,31,5,32,6))
))
d <- data.frame(Id,Mild,NewDrug,Week,Normal)
ftable(xtabs(~Mild+NewDrug+Week+Normal))

# observation matrix and covariate matrix
Y <- matrix(d$Normal,nrow=3,ncol=340)
A <- matrix(0,nrow=3,ncol=340)
A[1,] <- 1
A[2,] <- d$Mild[d$Week==1]
A[3,] <- d$NewDrug[d$Week==1]
rownames(A) <- c("Const","Mild","NewDrug")

# nmf with covariates
library(nmfkc)
result <- nmfkc(Y,A,Q=2,epsilon=1e-8)
plot(result)

# unique covariates and coefficient matrix
(A0 <- t(unique(t(A))))
B <- result$C %*% A0

# fitted values for coefficient matrix
t <- unique(d$Week)
plot(t,Y[,1],ylim=range(Y),type="n",ylab="Normal",xlab="Week")
XB <- result$X %*% B
mycol <- c(2,7,4,5)
for(j in 1:ncol(XB))lines(t,XB[,j],col=mycol[j],lwd=5)
legend("bottomright",
       title="Diagnosis & Treatment",
       legend=c("Mild & Standard","Mild & NewDrug","Sever & Standard","Sever & NewDrug"),
       fill=mycol)
```

## 8.  Autoregression: AirPassengers

This example demonstrates how to apply the NMF-VAR (Vector Autoregression) model to the famous "AirPassengers" dataset. We will train the model on historical data and forecast future passenger numbers for the next 2 years.

``` r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
# install.packages("DOT")
library(nmfkc)

# --- 1. Data Preparation ---
# Load the AirPassengers dataset
data(AirPassengers)
d <- AirPassengers

# Create a time vector for column names
time <- time(ts(1:length(d), start = c(1949, 1), frequency = 12))
time.vec <- round(as.vector(t(time)), 2)

# Convert to a matrix (1 x N) and apply log10 transformation
Y0 <- log10(matrix(as.vector(d), nrow = 1))
colnames(Y0) <- time.vec
rownames(Y0) <- "Passengers (log10)"

# --- 2. Model Selection (Lag Order) ---
# Determine the optimal lag order (degree) using Cross-Validation.
# We evaluate degrees from 1 to 16 with Q=1 basis.
cv_res <- nmfkc.ar.degree.cv(Y0, Q = 1, degree = 1:16)

# Use the optimal degree found (or manually set D=12 for monthly seasonality)
D <- cv_res$degree
# D <- 12 

# --- 3. Model Fitting (NMF-VAR) ---
# Set parameters
Q <- 1    # Number of latent bases (Rank)

# Construct matrices for the AR model (Y and A)
ar_data <- nmfkc.ar(Y0, degree = D, intercept = TRUE)
Y <- ar_data$Y
A <- ar_data$A

# Fit the NMF model
res <- nmfkc(Y = Y, A = A, rank = Q, prefix = "Factor", epsilon = 1e-9, maxit = 300000)

# Check the goodness of fit (R-squared)
print(paste("R-squared:", round(res$r.squared, 4)))

# Check stationarity (Spectral radius < 1 indicates a stationary process)
print(nmfkc.ar.stationarity(res))

# --- 4. Forecasting ---
# Forecast 24 steps (2 years) ahead
h <- 24
pred_res <- nmfkc.ar.predict(x = res, Y = Y0, n.ahead = h)

# Back-transform predicted values from log10 scale to original scale
pred_val <- 10^as.vector(pred_res$pred)

# --- 5. Visualization ---
# Create the time axis for the forecast period manually
start_year <- 1949
freq <- 12
total_len <- length(d)
pred_time <- seq(from = start_year + total_len/freq, by = 1/freq, length.out = h)

# Prepare data for smooth plotting (connect the lines)
last_obs_time <- tail(as.numeric(colnames(Y0)), 1)
last_obs_val  <- tail(10^as.vector(Y0), 1)

plot_pred_time <- c(last_obs_time, pred_time)
plot_pred_val  <- c(last_obs_val, pred_val)

# Set plot limits
xlim_range <- range(c(as.numeric(colnames(Y0)), pred_time))
ylim_range <- range(c(10^as.vector(Y0), pred_val))

# Draw the plot
# 1. Observed data (Black)
plot(as.numeric(colnames(Y0)), 10^as.vector(Y0), type = "l", col = "black",
     xlim = xlim_range, ylim = ylim_range, lwd = 1,
     xlab = "Time", ylab = "AirPassengers", main = "NMF-VAR Forecast (h=24)")

# 2. Fitted values during training (Red)
lines(as.numeric(colnames(Y)), 10^as.vector(res$XB), col = "red", lwd = 2)

# 3. Future forecast (Blue)
lines(plot_pred_time, plot_pred_val, col = "blue", lwd = 2)

# Add legend
legend("topleft", legend = c("Observed", "Fitted", "Forecast"),
       col = c("black", "red", "blue"), lty = 1, lwd = c(1, 2, 2))
```
 
## 9:  Vector Autoregression: Canada

``` r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
# install.packages("vars")
# install.packages("DOT")
library(nmfkc)

library(vars)
d0 <- Canada
time <- as.vector(time(ts(1:nrow(d0),start=c(1980,1),frequency=4)))
rownames(d0) <- time
dd <- apply(d0,2,diff)
dn <- nmfkc.normalize(dd)

Y0 <- t(dn)
Q <- 2; D <- 1
a <- nmfkc.ar(Y0,degree=D,intercept=T)
Y <- a$Y
A <- a$A
res <- nmfkc(Y=Y,A=A,Q=Q,prefix="Condition",epsilon=1e-6)
res$r.squared

# spectral radius of the companion matrix 
nmfkc.ar.stationarity(res)

# visualize relation between variables
script <- nmfkc.ar.DOT(res,intercept=T,digits=2)
# cat(script)
library(DOT)
dot(script,file="Canada_dot.ps")

# Soft clustering of time trend
barplot(res$B.prob,col=1:Q+1,border=1:Q+1)
legend("topright",legend=colnames(res$X),fill=1:Q+1,bg="white")

# fitted curve
t1 <- t(res$XB)
dn1 <- nmfkc.denormalize(t1,dd)
d1 <- apply(rbind(d0[1+D,],dn1),2,cumsum)
rownames(d1)[1] <- rownames(d0)[1+D]
for(p in 1:nrow(Y)){
  plot(as.numeric(rownames(d0)),d0[,p],type="l",col=8,ylab=rownames(Y)[p])
  lines(as.numeric(rownames(d1)),d1[,p],col=2,lwd=3)
}
```

## 10: Classification: Iris
``` r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
library(nmfkc)

label <- iris$Species
table(label)
Y <- nmfkc.class(label)
U <- t(nmfkc.normalize(iris[,-5]))
dim(U)
range(U)

Q <- length(unique(label))
res.beta <- nmfkc.kernel.beta.nearest.med(U)
beta.med <- res.beta$dist_median
betas <- beta.med*10^(-2:1)

res.cv <- nmfkc.kernel.beta.cv(Y,Q,U)
best.beta <- res.cv$beta
A <- nmfkc.kernel(U,beta=best.beta)
res <- nmfkc(Y=Y,A=A,Q=Q,prefix="Class")
res$r.squared
res$X

fitted.label <- predict(res,type="class")
(f <- table(fitted.label,label))
100*sum(diag(f))/sum(f)
```


## 11: Classification: penguins
``` r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
library(nmfkc)

library(palmerpenguins)
d <- penguins
index <- complete.cases(d)
d <- d[index,]
label <- d$species
table(label)
Y <- nmfkc.class(label)
U <- t(nmfkc.normalize(d[,3:6]))
dim(U)
range(U)

Q <- length(unique(label))
res.beta <- nmfkc.kernel.beta.nearest.med(U)
beta.med <- res.beta$dist_median
betas <- beta.med*10^(-2:1)

res.cv <- nmfkc.kernel.beta.cv(Y,Q,U)
best.beta <- res.cv$beta
A <- nmfkc.kernel(U,beta=best.beta)
res <- nmfkc(Y=Y,A=A,Q=Q,prefix="Class")
res$r.squared
res$X

fitted.label <- predict(res,type="class")
(f <- table(fitted.label,label))
100*sum(diag(f))/sum(f)
```



# Author

-   Kenichi Satoh, [homepage](https://sites.google.com/view/ksatoh/english)
