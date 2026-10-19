# Package index

## All functions

- [`coef(`*`<nmf>`*`)`](https://ksatohds.github.io/nmfkc/reference/coef.nmf.md)
  [`coef(`*`<nmf.ffb>`*`)`](https://ksatohds.github.io/nmfkc/reference/coef.nmf.md)
  : Extract coefficients from NMF models

- [`coef(`*`<nmf.gmm>`*`)`](https://ksatohds.github.io/nmfkc/reference/coef.nmf.gmm.md)
  : Extract the covariate-coefficient matrix from an NMF-GMM fit

- [`fitted(`*`<nmf>`*`)`](https://ksatohds.github.io/nmfkc/reference/fitted.nmf.md)
  [`fitted(`*`<nmfae>`*`)`](https://ksatohds.github.io/nmfkc/reference/fitted.nmf.md)
  [`fitted(`*`<nmfre>`*`)`](https://ksatohds.github.io/nmfkc/reference/fitted.nmf.md)
  [`fitted(`*`<nmf.ffb>`*`)`](https://ksatohds.github.io/nmfkc/reference/fitted.nmf.md)
  : Extract fitted values from NMF models

- [`fitted(`*`<nmf.gmm>`*`)`](https://ksatohds.github.io/nmfkc/reference/fitted.nmf.gmm.md)
  : Fitted (responsibility-averaged) reconstruction of an NMF-GMM fit

- [`nmf.cluster.criteria()`](https://ksatohds.github.io/nmfkc/reference/nmf.cluster.criteria.md)
  : Sample-clustering quality across ranks

- [`nmf.cluster.flow()`](https://ksatohds.github.io/nmfkc/reference/nmf.cluster.flow.md)
  : Cluster-flow (alluvial) diagram across a sequence of fits

- [`nmf.ffb.DOT()`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.DOT.md)
  : Generate a Graphviz DOT Diagram for an NMF-FFB Model

- [`nmf.ffb()`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md) :
  NMF-FFB: non-negative matrix factorization with latent feedback

- [`nmf.ffb.cv()`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.cv.md)
  : Cross-Validation for NMF-FFB

- [`nmf.ffb.diagnostics()`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.diagnostics.md)
  : Structural Diagnostics of a Fitted NMF-FFB Feedback Matrix

- [`nmf.ffb.ecv()`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.ecv.md)
  : Choose the number of factors for NMF-FFB by element-wise
  cross-validation

- [`nmf.ffb.inference()`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.inference.md)
  : Bootstrap inference for NMF-FFB (X held fixed)

- [`nmf.ffb.split()`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.split.md)
  : Heuristic Variable Splitting for NMF-FFB

- [`nmf.ffb.test()`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.test.md)
  : Test the feed-forward null of an NMF-FFB fit

- [`nmf.gmm()`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md) :
  Fit NMF-GMM: a Gaussian-mixture latent-class extension of NMF with
  covariates

- [`nmf.gmm.inference()`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.inference.md)
  : Statistical inference for an NMF-GMM fit (given basis)

- [`nmf.gmm.select()`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.select.md)
  : Choose the number of mixture components K for NMF-GMM

- [`nmf.gmm.twostage()`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.twostage.md)
  : Two-stage (adjust-then-cluster) baseline for NMF-GMM

- [`nmf.rrr.DOT()`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.DOT.md)
  : DOT graph visualization for nmfae objects

- [`nmf.rrr()`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md) :
  Three-Layer Non-negative Matrix Factorization (NMF-AE)

- [`nmf.rrr.cv()`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.cv.md)
  : Sample-wise k-fold Cross-Validation for nmfae

- [`nmf.rrr.ecv()`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.ecv.md)
  : Element-wise Cross-Validation for nmfae (Wold's CV)

- [`nmf.rrr.heatmap()`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.heatmap.md)
  : Heatmap visualization of nmfae factor matrices

- [`nmf.rrr.inference()`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.inference.md)
  : Statistical Inference for NMF-AE Parameter Matrix

- [`nmf.rrr.kernel.beta.cv()`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.kernel.beta.cv.md)
  : Optimize kernel beta for nmfae by cross-validation

- [`nmf.rrr.rank()`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.rank.md)
  : Rank selection for nmfae (paired rank, concise diagnostics)

- [`nmf.rrr.rename()`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.rename.md)
  : Rename decoder and encoder bases

- [`nmf.rrr.signed()`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.md)
  : Signed-Bottleneck NMF-AE: Three-Layer NMF-AE with Signed Bottleneck

- [`nmf.rrr.signed.ecv()`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.ecv.md)
  : Element-wise Cross-Validation for Signed-Bottleneck NMF-AE

- [`nmf.rrr.signed.heatmap()`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.heatmap.md)
  : Heatmap visualization of nmfae.signed factor matrices

- [`nmf.rrr.signed.inference()`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.inference.md)
  : Statistical Inference for Signed-Bottleneck NMF-AE Signed Bottleneck

- [`nmf.rrr.signed.rank()`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.rank.md)
  : Rank selection for nmfae.signed (paired rank, concise diagnostics)

- [`nmf.rrr.signed.rename()`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.rename.md)
  : Rename Resp/Cov labels on nmfae.signed objects

- [`nmfkc.DOT()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.DOT.md)
  : Generate Graphviz DOT Scripts for NMF or NMF-with-Covariates Models

- [`nmfkc()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) :
  Optimize NMF with kernel covariates (Full Support for Missing Values)

- [`nmfkc.ar.DOT()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.DOT.md)
  : Generate a Graphviz DOT Diagram for NMF-AR / NMF-VAR Models

- [`nmfkc.ar()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.md)
  : Construct observation and covariate matrices for a vector
  autoregressive model

- [`nmfkc.ar.degree.cv()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.degree.cv.md)
  : Optimize lag order for the autoregressive model

- [`nmfkc.ar.latent()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.latent.md)
  : Latent transition matrices of an NMF-VAR model

- [`nmfkc.ar.latent.inference()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.latent.inference.md)
  : Bootstrap inference for the latent VAR of an NMF-VAR model

- [`nmfkc.ar.predict()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.predict.md)
  : Forecast future values for NMF-VAR model

- [`nmfkc.ar.stationarity()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.stationarity.md)
  : Stationarity of an NMF-VAR model

- [`nmfkc.ard()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ard.md)
  : Automatic relevance determination for NMF rank (experimental)

- [`nmfkc.bicv()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.bicv.md)
  : Bi-cross-validation for NMF rank selection

- [`nmfkc.class()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.class.md)
  : Create a class (one-hot) matrix from a categorical vector

- [`nmfkc.consensus()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.consensus.md)
  : Consensus-clustering rank selection for NMF (Brunet 2004)

- [`nmfkc.criterion()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.criterion.md)
  : Compute model selection criteria for a fitted nmfkc model

- [`nmfkc.cv()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.cv.md)
  : Perform k-fold cross-validation for NMF with kernel covariates

- [`print(`*`<nmfkc.ecv>`*`)`](https://ksatohds.github.io/nmfkc/reference/nmfkc.cv.methods.md)
  [`plot(`*`<nmfkc.ecv>`*`)`](https://ksatohds.github.io/nmfkc/reference/nmfkc.cv.methods.md)
  [`print(`*`<nmfkc.bicv>`*`)`](https://ksatohds.github.io/nmfkc/reference/nmfkc.cv.methods.md)
  [`plot(`*`<nmfkc.bicv>`*`)`](https://ksatohds.github.io/nmfkc/reference/nmfkc.cv.methods.md)
  [`print(`*`<nmfkc.cv>`*`)`](https://ksatohds.github.io/nmfkc/reference/nmfkc.cv.methods.md)
  [`plot(`*`<nmfkc.cv>`*`)`](https://ksatohds.github.io/nmfkc/reference/nmfkc.cv.methods.md)
  : Print/plot methods for nmfkc cross-validation objects

- [`nmfkc.denormalize()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.denormalize.md)
  : Denormalize a matrix from \\\[0,1\]\\ back to its original scale

- [`nmfkc.ecv()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ecv.md)
  : Perform Element-wise Cross-Validation (Wold's CV)

- [`nmfkc.inference()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.inference.md)
  : Statistical inference for the parameter matrix C (Theta)

- [`nmfkc.kernel()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.md)
  : Create a kernel matrix from covariates

- [`nmfkc.kernel.beta.cv()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.beta.cv.md)
  : Optimize beta of the Gaussian kernel function by cross-validation

- [`nmfkc.kernel.beta.nearest.med()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.beta.nearest.med.md)
  : Estimate Gaussian/RBF kernel parameter beta from covariates
  (supports landmarks)

- [`nmfkc.kernel.gaussian()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.gaussian.md)
  : Create a Gaussian kernel matrix from covariates

- [`nmfkc.kernel.gram()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.gram.md)
  : Block-wise Gram accumulation of Nystroem kernel covariates for large
  N

- [`nmfkc.net.DOT()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.DOT.md)
  : Generate a Graphviz DOT Diagram for a Symmetric NMF Network

- [`nmfkc.net()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.md)
  : Symmetric NMF for networks (tri / bi / signed)

- [`nmfkc.net.ecv()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.ecv.md)
  : Element-wise cross-validation for nmfkc.net (upper-triangle folds)

- [`nmfkc.net.inference()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.inference.md)
  : Statistical Inference for Symmetric NMF Parameters

- [`nmfkc.net.rank()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.rank.md)
  : Rank selection for nmfkc.net (concise diagnostics)

- [`nmfkc.normalize()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.normalize.md)
  : Normalize a matrix to the range \\\[0,1\]\\

- [`nmfkc.rank()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rank.md)
  : Rank selection diagnostics with graphical output

- [`nmfkc.residual.plot()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.residual.plot.md)
  : Plot Diagnostics: Original, Fitted, and Residual Matrices as
  Heatmaps

- [`nmfkc.rff.beta.cv()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rff.beta.cv.md)
  : Select the bandwidth of random-feature covariates by
  cross-validation

- [`nmfkc.rff.positive()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rff.positive.md)
  : Positive (non-negative) random features for the Gaussian kernel

- [`nmfkc.rff.positive.gram()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rff.positive.gram.md)
  : Block-wise Gram accumulation of positive random features for large N

- [`nmfkc.signed()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md)
  : NMF-KC with signed covariate matrix

- [`nmfkc.signed.cv()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.cv.md)
  : Column-wise k-fold cross-validation for nmfkc.signed

- [`nmfkc.signed.ecv()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.ecv.md)
  : Element-wise cross-validation for nmfkc.signed

- [`nmfkc.signed.rank()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.rank.md)
  : Rank selection for nmfkc.signed (concise diagnostics)

- [`nmfkc.signed.rff()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.rff.md)
  : Random Fourier Features for nmfkc.signed()

- [`nmfkc.signed.rff.gram()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.rff.gram.md)
  : Block-wise Gram accumulation of Random Fourier Features for large N

- [`nmfre()`](https://ksatohds.github.io/nmfkc/reference/nmfre.md) :
  Non-negative Matrix Factorization with Random Effects

- [`nmfre.ecv()`](https://ksatohds.github.io/nmfkc/reference/nmfre.ecv.md)
  : NMF-RE element-wise cross-validation for rank selection

- [`nmfre.inference()`](https://ksatohds.github.io/nmfkc/reference/nmfre.inference.md)
  : Statistical inference for the coefficient matrix C from NMF-RE

- [`plot(`*`<nmf.cluster.criteria>`*`)`](https://ksatohds.github.io/nmfkc/reference/plot.nmf.cluster.criteria.md)
  : Plot clustering-quality criteria across a sequence of fits

- [`plot(`*`<nmf.cluster.flow>`*`)`](https://ksatohds.github.io/nmfkc/reference/plot.nmf.cluster.flow.md)
  : Plot a cluster-flow (alluvial) diagram

- [`plot(`*`<nmf.gmm>`*`)`](https://ksatohds.github.io/nmfkc/reference/plot.nmf.gmm.md)
  : Plot method for nmf.gmm objects

- [`plot(`*`<nmf.gmm.select>`*`)`](https://ksatohds.github.io/nmfkc/reference/plot.nmf.gmm.select.md)
  : Plot method for nmf.gmm.select objects

- [`plot(`*`<nmf.rank>`*`)`](https://ksatohds.github.io/nmfkc/reference/plot.nmf.rank.md)
  : Plot a rank-selection (nmf.rank) object

- [`plot(`*`<nmfae>`*`)`](https://ksatohds.github.io/nmfkc/reference/plot.nmfae.md)
  :

  `plot.nmfae` displays the convergence trajectory of the objective
  function across iterations. The title shows the achieved \\R^2\\.

- [`plot(`*`<nmfkc.DOT>`*`)`](https://ksatohds.github.io/nmfkc/reference/plot.nmfkc.DOT.md)
  : Plot method for nmfkc.DOT objects

- [`plot(`*`<nmfkc>`*`)`](https://ksatohds.github.io/nmfkc/reference/plot.nmfkc.md)
  :

  Plot method for objects of class `nmfkc`

- [`plot(`*`<nmfkc.ar.latent>`*`)`](https://ksatohds.github.io/nmfkc/reference/plot.nmfkc.ar.latent.md)
  : Plot the latent dynamics of an NMF-VAR model

- [`plot(`*`<nmfkc.ar.stationarity>`*`)`](https://ksatohds.github.io/nmfkc/reference/plot.nmfkc.ar.stationarity.md)
  : Plot the stationarity bracket of an NMF-VAR model

- [`plot(`*`<nmfkc.ard>`*`)`](https://ksatohds.github.io/nmfkc/reference/plot.nmfkc.ard.md)
  : Plot method for nmfkc.ard objects

- [`plot(`*`<nmfkc.consensus>`*`)`](https://ksatohds.github.io/nmfkc/reference/plot.nmfkc.consensus.md)
  : Plot a consensus rank-selection (nmfkc.consensus) object

- [`plot(`*`<nmfkc.signed>`*`)`](https://ksatohds.github.io/nmfkc/reference/plot.nmfkc.signed.md)
  : Plot method for nmfkc.signed (convergence)

- [`plot(`*`<nmfre>`*`)`](https://ksatohds.github.io/nmfkc/reference/plot.nmfre.md)
  [`plot(`*`<nmf.ffb>`*`)`](https://ksatohds.github.io/nmfkc/reference/plot.nmfre.md)
  : Plot convergence diagnostics for NMF models

- [`plot(`*`<nmfre.ecv>`*`)`](https://ksatohds.github.io/nmfkc/reference/plot.nmfre.ecv.md)
  : Plot method for nmfre.ecv objects

- [`predict(`*`<nmf.gmm>`*`)`](https://ksatohds.github.io/nmfkc/reference/predict.nmf.gmm.md)
  : Class assignments / responsibilities from an NMF-GMM fit

- [`predict(`*`<nmfkc>`*`)`](https://ksatohds.github.io/nmfkc/reference/predict.nmfkc.md)
  :

  Prediction method for objects of class `nmfkc`

- [`predict(`*`<nmfkc.signed>`*`)`](https://ksatohds.github.io/nmfkc/reference/predict.nmfkc.signed.md)
  : Predict method for nmfkc.signed

- [`predict(`*`<nmfre>`*`)`](https://ksatohds.github.io/nmfkc/reference/predict.nmfre.md)
  : Predict method for nmfre objects

- [`print(`*`<nmf>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.nmf.md)
  : Print method for fitted NMF models

- [`print(`*`<nmf.cluster.criteria>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.nmf.cluster.criteria.md)
  : Print method for nmf.cluster.criteria objects

- [`print(`*`<nmf.cluster.flow>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.nmf.cluster.flow.md)
  : Print method for nmf.cluster.flow objects

- [`print(`*`<nmf.ffb.test>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.nmf.ffb.test.md)
  : Print a calibrated test of the feed-forward null

- [`print(`*`<nmf.gmm>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.nmf.gmm.md)
  : Print an NMF-GMM fit

- [`print(`*`<nmf.gmm.select>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.nmf.gmm.select.md)
  : Print method for nmf.gmm.select objects

- [`print(`*`<nmf.inference>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.nmf.inference.md)
  : Print method for NMF inference objects

- [`print(`*`<nmf.rank>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.nmf.rank.md)
  : Print method for rank-selection (nmf.rank) objects

- [`print(`*`<nmfae.ecv>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.nmfae.ecv.md)
  : Print method for nmf.rrr element-wise CV results

- [`print(`*`<nmfkc.DOT>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.nmfkc.DOT.md)
  : Print method for nmfkc.DOT objects

- [`print(`*`<nmfkc.ar.latent>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.nmfkc.ar.latent.md)
  : Print method for nmfkc.ar.latent objects

- [`print(`*`<nmfkc.ar.latent.inference>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.nmfkc.ar.latent.inference.md)
  : Print method for nmfkc.ar.latent.inference objects

- [`print(`*`<nmfkc.ar.stationarity>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.nmfkc.ar.stationarity.md)
  : Print method for nmfkc.ar.stationarity objects

- [`print(`*`<nmfkc.ard>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.nmfkc.ard.md)
  : Print method for nmfkc.ard objects

- [`print(`*`<nmfkc.consensus>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.nmfkc.consensus.md)
  : Print method for nmfkc.consensus objects

- [`print(`*`<nmfkc.gram>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.nmfkc.gram.md)
  : Print method for nmfkc.gram

- [`print(`*`<nmfre.ecv>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.nmfre.ecv.md)
  : Print method for nmfre.ecv objects

- [`print(`*`<summary.nmf.ffb>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.summary.nmf.ffb.md)
  : Print method for summary.nmf.ffb objects

- [`print(`*`<summary.nmf.gmm>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.summary.nmf.gmm.md)
  : Print method for summary.nmf.gmm objects

- [`print(`*`<summary.nmfkc>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.summary.nmfkc.md)
  :

  Print method for `summary.nmfkc` objects

- [`print(`*`<summary.nmfkc.inference>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.summary.nmfkc.inference.md)
  : Print method for summary.nmfkc.inference objects

- [`print(`*`<summary.nmfkc.net>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.summary.nmfkc.net.md)
  : Print method for summary.nmfkc.net objects

- [`print(`*`<summary.nmfkc.net.inference>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.summary.nmfkc.net.inference.md)
  : Print method for summary.nmfkc.net.inference objects

- [`print(`*`<summary.nmfkc.net.signed>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.summary.nmfkc.net.signed.md)
  : Print method for summary.nmfkc.net.signed objects

- [`print(`*`<summary.nmfkc.signed>`*`)`](https://ksatohds.github.io/nmfkc/reference/print.summary.nmfkc.signed.md)
  : Print method for summary.nmfkc.signed

- [`residuals(`*`<nmf>`*`)`](https://ksatohds.github.io/nmfkc/reference/residuals.nmf.md)
  [`residuals(`*`<nmfae>`*`)`](https://ksatohds.github.io/nmfkc/reference/residuals.nmf.md)
  [`residuals(`*`<nmfre>`*`)`](https://ksatohds.github.io/nmfkc/reference/residuals.nmf.md)
  [`residuals(`*`<nmf.ffb>`*`)`](https://ksatohds.github.io/nmfkc/reference/residuals.nmf.md)
  : Extract residuals from NMF models

- [`residuals(`*`<nmf.gmm>`*`)`](https://ksatohds.github.io/nmfkc/reference/residuals.nmf.gmm.md)
  : Residuals from an NMF-GMM fit

- [`summary(`*`<nmf.ffb>`*`)`](https://ksatohds.github.io/nmfkc/reference/summary.nmf.ffb.md)
  : Summary method for nmf.ffb objects

- [`summary(`*`<nmf.gmm>`*`)`](https://ksatohds.github.io/nmfkc/reference/summary.nmf.gmm.md)
  : Summary of an NMF-GMM fit

- [`summary(`*`<nmfkc>`*`)`](https://ksatohds.github.io/nmfkc/reference/summary.nmfkc.md)
  :

  Summary method for objects of class `nmfkc`

- [`summary(`*`<nmfkc.inference>`*`)`](https://ksatohds.github.io/nmfkc/reference/summary.nmfkc.inference.md)
  : Summary method for nmfkc.inference objects

- [`summary(`*`<nmfkc.net>`*`)`](https://ksatohds.github.io/nmfkc/reference/summary.nmfkc.net.md)
  : Summary method for nmfkc.net objects

- [`summary(`*`<nmfkc.net.inference>`*`)`](https://ksatohds.github.io/nmfkc/reference/summary.nmfkc.net.inference.md)
  : Summary method for nmfkc.net.inference objects

- [`summary(`*`<nmfkc.net.signed>`*`)`](https://ksatohds.github.io/nmfkc/reference/summary.nmfkc.net.signed.md)
  : Summary method for nmfkc.net.signed objects

- [`summary(`*`<nmfkc.signed>`*`)`](https://ksatohds.github.io/nmfkc/reference/summary.nmfkc.signed.md)
  : Summary method for nmfkc.signed

- [`summary(`*`<nmfre>`*`)`](https://ksatohds.github.io/nmfkc/reference/summary.nmfre.md)
  [`print(`*`<summary.nmfre>`*`)`](https://ksatohds.github.io/nmfkc/reference/summary.nmfre.md)
  :

  Summary method for objects of class `nmfre`
