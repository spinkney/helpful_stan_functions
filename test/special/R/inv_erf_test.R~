library(cmdstanr)
fp <- file.path("test/special/inv_erf_test.stan")
mod <- cmdstan_model(fp,  
                     include_paths = c("./functions/special",
                                       "./functions/void"))

mod_out <- mod$sample(
  fixed_param = TRUE
)
