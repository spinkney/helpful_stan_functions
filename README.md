<p align="center">
<img src="pretty/helpful_logo2_crop.png" alt="drawing" width="250"/>
 </p>

# Unoffical Stan UDF Repo 

A home for all the **U**ser **D**efined **F**unctions that the Stan community uses or finds useful that haven't made their way into Stan-math. This is Stan with batteries included!

There are three main directories. A UDF must be defined in the `functions` directory. There must be at least an example or a test that goes along with the UDF.

## Documentation 
Documentation for the functions is at hosted on github pages and can be found at 

<a href="https://spinkney.github.io/helpful_stan_functions/"> Helpful Stan Functions Documentation </a>

The repository is setup with github actions to run `doxygen` automatically when a pull-request or push to the `main` branch is triggered. 

### Auto-documentation with `doxygen`

All the `function` files must conform to `javadoc` style documentation. Although `doxygen` will run with a few other styles, this style will be enforced so that all the code is consistently written. See more at:

[How to use doxygen documentation for the project?](https://spinkney.github.io/helpful_stan_functions/doxygen_doc.html)

## Directory Structure

#### `functions`

Contains sub-directories with names indicating the type of UDF. The `functions` folder contains only `.stan` files. For example, all the probability distributions can be found by going into `functions` -> `distribution`. Since UDFs cannot have different signatures or type overloading, each version is contained in the same `.stan` file with a clarifying name (e.g. x_fun and x_fun_vectorized).

If you have a UDF that does not fit in the classification scheme please create an issue for discussion. 

#### `examples`

The examples directory will house an example of each of the functions. The plan is to support many programming languages, however, it is currently all in `R`. 

The directory contains subdirectories for `stan` files and for the programming language that calls the example in Stan. For example, if you are interested in seeing how `lognormal_qf.stan` can be used then go into the examples directory and traverse down the `quantile_function` examples. That is, `examples` -> `quantile_function` -> `R` -> `lognormal_qf_example.R`. 

All the example file names will be the same name as the function file with `_example` attached to the name. 

If you would like to add an example script, a sufficient example will simulate data, pass it to the `stan` example, and have a minimal extraction to compare to the simulated data generating process.

#### `test`

Contains a directory for the `stan` file example that calls the UDF function in the functions directory and the language (e.g. R) to run the test. The naming convention in the test directories follows the same style as in the `examples` directory but with `_test` attached to each file.

Tests are simpler than examples and do not contain data. The program is run with stans `fixed_param` algorithm. The `unit_test` functions such as `expect_equal` will be called to compare to the expected output.


## Stargazers over time
<img src=https://user-images.githubusercontent.com/23704057/143876874-0f6a1f2a-6466-4e4b-a1bd-99d69383a9c3.png width="3" height="3"> 
