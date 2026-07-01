# Generates the bundled example dataset `ivc_example`.
source("R/internal.R"); source("R/dgp.R")
set.seed(2025)
ivc_example <- ivc_simulate(n = 1000, gamma = 1, kappa = 0.2, tau = 0.5, seed = 2025)
ivc_example$U <- NULL  # drop the latent confounder; not observable in practice
save(ivc_example, file = "data/ivc_example.rda", compress = "xz")
