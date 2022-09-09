#!/usr/bin/env Rscript

library(mfittrial)
library(data.table)
library(parallel)
library(bayestestR)
library(optparse)

option_list = list(
  make_option(c("-c", "--cores"), type="integer", default=10,
              help="number of cores to use", metavar="character"),
  make_option(c("-n", "--nsim"), type="integer", default=10,
              help="number of simulations to run under each configuration", metavar="character")
);

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
# To run call, e.g.: Rscript null.R -c 15 -n 2000
opt <- parse_args(OptionParser(option_list = option_list))

num_cores <- opt$cores
num_sims  <- opt$nsim

# ----- Configurations -----
dropout <- 0.8

mean_cfg <- list(
  rep(40, 4),
  c(40, 45, 40, 40),
  c(40, 42.5, 40, 40),
  c(40, 42.5, 42.5, 40),
  c(40, 42.5, 42.5, 42.5)
)

cfg <- data.table(
  sims = num_sims,
  nsubj = list(c(100, 200, 300, 400) * dropout),
  eff_eps = 0.975,
  fut_eps = 0.950,
  sup_eps = 0.975,
  brar = 2,
  allow_stopping = FALSE,
  means = mean_cfg,
  prior = list(c(int_prior_mean = 40, int_prior_sd = 5, b_prior_sd = 5, df = 3, scale = 5))
)
cfg2 <- copy(cfg)[, `:=`(eff_eps = 0.98, sup_eps = 0.98)]
cfg3 <- copy(cfg)[, `:=`(eff_eps = 0.99, sup_eps = 0.99)]
cfg <- rbind(cfg, cfg2, cfg3)

# ----- Which configurations do we want to run? -----

run_row <- seq_len(nrow(cfg))

# ----- Loop over configurations and save results -----

for(z in run_row) {

  start_time <- Sys.time()

  res <- mclapply(1:cfg[z][["sims"]], function(j) {
    # Generate data for trial
    dat <- simulate_continuous_outcome(
      nsubj = max(unlist(cfg[z][["nsubj"]])),
      means = unlist(cfg[z][["means"]])
    )
    simulate_trial_with_control(
      dat = dat,
      n_seq = unlist(cfg[z][["nsubj"]]),
      eff_eps = cfg[z][["eff_eps"]],
      sup_eps = cfg[z][["sup_eps"]],
      brar = cfg[z][["brar"]],
      allow_stopping = cfg[z][["allow_stopping"]],
      prior = unlist(cfg[z][["prior"]]))
  }, mc.cores = num_cores)
  resl <- rbindlist(res, idcol = "trial")
  resl[, analysis := as.numeric(analysis)]

  end_time <- Sys.time()

  saveRDS(list(cfg = cfg[z], res = resl, runtime = end_time - start_time), paste0("~/out_files/mfit_sims/with_control_", formatC(z, width = 2, flag = "0"), ".rds"))
}
