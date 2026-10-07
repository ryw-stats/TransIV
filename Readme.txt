TransIV code
============

This folder contains the simulation and real-data analyses for TransIV.
Read this file for setup and the code map. Use Workflow.txt for the exact
run order and data-file handoffs. Study datasets and trained models are not
included.

Choose the workflow you need; they do not all have to run in sequence:
  Linear simulation     Hapgen preprocessing -> five notebooks -> conversion
                        -> IV analyses and figures.
  Nonlinear simulation  Nonlinear notebook -> IV analyses and figures.
  Debiasing simulation  Sim_Debias.R; no external study data required.
  Real-data analysis    Align one target population -> preprocessing
                        -> source-model fitting -> target estimation.

Running the code
----------------
All commands below are R console commands. Start with Code/ as the working
directory; from the project root, run setwd("Code") once. Paths inside an
analysis script are relative to that script's directory. Run it with
chdir = TRUE, for example, after preparing its inputs:

  source("Simulation/Linear-Sim/Linear-sim_Target_Transfer.R", chdir = TRUE)

Alternatively, set the working directory to the script's folder before using
Rscript or running selected lines in RStudio. Scripts no longer depend on
RStudio's active editor to choose the working directory. Keep Code/R/ in its
current location: the entry scripts source the shared functions automatically.

Render an R Markdown notebook from Code/ with:

  rmarkdown::render("Simulation/Linear-Sim/GenoSim-h=0.Rmd")

Rendering executes the notebook, including model training. It is not a
syntax-only check. Most analyses require the input files described in
Workflow.txt before this command will succeed.

Use a fresh R session for each workflow. The linear and nonlinear helper
modules have workflow-specific implementations of some function names and
should not be loaded together into the same environment.

Dependencies
------------
- Linear generation: glmnet, MASS, withr, and base-R parallel.
- Nonlinear generation: keras, tensorflow, glmnet, caret, MASS, randomForest,
  gbm, and ggplot2. Neural networks require a configured TensorFlow backend
  compatible with the existing keras API used in the code.
- IV estimation and figures: MASS, parallel, ggplot2, latex2exp, viridis;
  the linear and debiasing figures also use reshape2 and patchwork.
- Real-data preprocessing: data.table and genio.
- Real-data source fitting and estimation: glmnet, MASS, and data.table.
- Notebook rendering: rmarkdown and knitr, with Pandoc for HTML output.

The linear notebooks use one process on Windows because mclapply does not
support multiple forked workers there. On other systems they retain the
original 20-worker setting. IV analysis scripts use 20 PSOCK workers; adjust
makeCluster(20) to suit available memory and CPU resources. Each PSOCK worker
sources the shared IV functions before evaluating simulation replicates.

Shared functions
----------------
R/model_functions.R
  Cross-validated lasso fitting, target lasso fitting, neural-network model
  construction, correlation-matrix construction, and residual sum of squares.
  train_lasso_model accepts penalty.factor so the LDL-CAD target analysis
  retains its seven unpenalized predictors. include_1se preserves the extra
  diagnostics used by the nonlinear workflow.

R/iv_functions.R
  Instrument-combination, variance, Hessian, bias-corrected IV, GEL, and DEEM
  objectives, plus the trimmed simulation summaries. Objective inputs are
  explicit arguments. abr_rmse(x, beta) takes the true effect explicitly.
  iv_variance_objective preserves each workflow's matrix multiplication
  order: simulations use projection_order = "left", while the real-data
  callers explicitly use "right" to reproduce their original rounding.

R/linear_simulation_functions.R
  Linear data generation, prediction/evaluation functions, one-replicate
  model fitting, and preparation of generated results for the IV analyses.
  The five notebooks pass their original coefficient-mixture settings into
  simulate_linear_data and run_linear_simulation.

R/nonlinear_simulation_functions.R
  Nonlinear data generation, source/target model fitting, hyperparameter
  selection, residual transfer fitting, and prediction/evaluation functions.
  Identical target/residual random-forest fitting uses train_rf_model.

All named functions include Roxygen-style purpose, parameter, and return-value
comments. These are ordinary source annotations; no package build is required.
Shared files define functions without loading datasets or launching analyses.

Analysis files and manuscript outputs
------------------------------------
Simulation/Linear-Sim/GenoSim-h=x.Rmd, x = 0, 0.25, 0.5, 0.75, 1
  Generate the data and predictions for Section 5.1.
Simulation/Linear-Sim/Prepare-Linear-Results.R
  Convert model_results_1.Rdata through model_results_5.Rdata into
  linear_1.Rdata through linear_5.Rdata, including analysis field names.
Simulation/Linear-Sim/Linear-sim_Target_Transfer.R
  Figure 3.
Simulation/Linear-Sim/Linear-sim_MM_CF.R
  Figure 4. Select the comparison setting with mod.
Simulation/Linear-Sim/Illustrative Example.R
  Supplementary Figure S3; reads linear_3.Rdata.

Simulation/Nonlinear-Sim/Hapgen-prePro.R
  Preprocess Hapgen2 genotype files and write SimData_EU.Rdata,
  SimData_AF.Rdata, and SimData_SA.Rdata. The linear notebooks read the first
  two files; the current nonlinear generator uses simulated Gaussian data.
Simulation/Nonlinear-Sim/Gen-Nonlinear.Rmd
  Section 5.2 generation and model evaluation. Runs 200 replicates and saves
  nonlinear_new.Rdata and model_evaluation_results_nonlinear.csv directly.
Simulation/Nonlinear-Sim/Nonlinear-sim_Target_Transfer.R
  Figure 5; reads nonlinear_new.Rdata.
Simulation/Nonlinear-Sim/Nonlinear-sim_MM_CF.R
  Supplementary Figure S2; reads nonlinear_new.Rdata. Set mod as required.
Simulation/Sim_Debias.R
  Supplementary Figure S1; generates its own data and saves Res_debiased.Rdata.

RealData/LDL-CAD/LDL-CAD-Prepro.R
  Preprocess LDL-CAD phenotypes, genotype files, and GWAS selections.
RealData/LDL-CAD/LDL-CAD-SourceModel-LASSO.R
  Train the LDL source models.
RealData/LDL-CAD/LDL-CAD-Transfer&Target.R
  Table 1 analyses.
RealData/BMI-SBP/BMI-SBP-Prepro.R
  Preprocess BMI-SBP phenotypes, genotype files, and GWAS selections.
RealData/BMI-SBP/BMI-SBP-SourceModel-LASSO.R
  Train the BMI source models.
RealData/BMI-SBP/BMI-SBP-Transfer&Target.R
  Table 2 analyses.

The supplied real-data scripts use different example target populations
between stages. They are not ready to run as an unchanged three-script chain.
Align the population and variant selection first, as described in Workflow.txt.