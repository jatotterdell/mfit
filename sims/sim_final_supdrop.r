#!/usr/bin/env Rscript

library(mfittrial)
library(data.table)
library(parallel)
library(bayestestR)
library(optparse)

option_list <- list(
    make_option(c("-c", "--cores"),
        type = "integer",
        default = 10,
        help = "number of cores to use",
        metavar = "character"
    ),
    make_option(c("-n", "--nsim"),
        type = "integer",
        default = 10,
        help = "number of simulations to run under each configuration",
        metavar = "character"
    )
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
# To run call, e.g.: Rscript null.R -c 15 -n 2000
opt <- parse_args(OptionParser(option_list = option_list))

num_cores <- opt$cores
num_sims <- opt$nsim

RNGkind("L'Ecuyer-CMRG")

# ----- Configurations -----

dropout <- 0.8
sigma <- 10
mean <- 35
mean_cfg <- list(
    rep(mean, 4),
    c(mean, mean + sigma, mean, mean),
    c(mean, mean + sigma / 2, mean, mean),
    c(mean, mean + sigma / 2, mean + sigma / 2, mean),
    c(mean, mean + sigma, mean + sigma / 2, mean),
    c(mean, mean + sigma, mean + sigma / 2, mean + sigma / 2),
    c(mean, mean - sigma / 2, mean - sigma / 2, mean - sigma / 2)
)

cfg <- data.table(
    sims = num_sims,
    nsubj = list(c(200, 300, 400) * dropout),
    eff_eps = 0.98,
    fut_eps = 1, # No futility
    sup_eps = 0.98,
    brar = c(
        rep(0, length(mean_cfg)),
        rep(2, length(mean_cfg)),
        rep(2, length(mean_cfg)),
        rep(2, length(mean_cfg))
    ),
    brar_k = 0.5,
    brar_min = c(
        rep(0, length(mean_cfg)),
        rep(0, length(mean_cfg)),
        rep(0.10, length(mean_cfg)),
        rep(0.15, length(mean_cfg))
    ),
    allow_stopping = TRUE,
    drop = "sup",
    means = rep(mean_cfg, times = 4),
    sigma = 10,
    trunc = FALSE,
    prior = list(
        c(
            int_prior_mean = mean,
            int_prior_sd = 5,
            b_prior_sd = 5,
            df = 3,
            scale = 5
        )
    ),
    ctr = contr.treatment
)

# ----- Which configurations do we want to run? -----

run_row <- seq_len(nrow(cfg))

# ----- Loop over configurations and save results -----

for (z in run_row) {
    start_time <- Sys.time()

    res <- mclapply(1:cfg[z][["sims"]], function(j) {
        # Generate data for trial
        dat <- simulate_continuous_outcome(
            nsubj = max(unlist(cfg[z][["nsubj"]])),
            means = unlist(cfg[z][["means"]]),
            sigma = cfg[z][["sigma"]],
            trunc = cfg[z][["trunc"]]
        )
        simulate_trial_with_control3(
            dat = dat,
            n_seq = unlist(cfg[z][["nsubj"]]),
            eff_eps = cfg[z][["eff_eps"]],
            sup_eps = cfg[z][["sup_eps"]],
            brar = cfg[z][["brar"]],
            brar_k = cfg[z][["brar_k"]],
            brar_min = cfg[z][["brar_min"]],
            allow_stopping = cfg[z][["allow_stopping"]],
            drop = cfg[z][["drop"]],
            prior = unlist(cfg[z][["prior"]]),
            ctr = cfg[z][["ctr"]][[1]]
        )
    }, mc.cores = num_cores)
    resl <- rbindlist(res, idcol = "trial")
    resl[, analysis := as.numeric(analysis)]

    end_time <- Sys.time()

    saveRDS(
        list(cfg = cfg[z], res = resl, runtime = end_time - start_time),
        paste0(
            "~/out_files/mfit_sims/sim_final_supdrop_",
            formatC(z, width = 2, flag = "0"),
            ".rds"
        )
    )
}
