library(cmdstanr)
fp <- file.path("test/special/stan/inv_erf_test.stan")
mod <- cmdstan_model(fp,  
                     include_paths = c("./functions/special",
                                       "./functions/unit_test"))

mod_out <- mod$sample(
  fixed_param = TRUE
)
