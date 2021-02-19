# Unoffical Stan UDF Repo

A home for all the **U**ser **D**efined **F**unctions that the Stan community uses or finds useful that haven't made their way into Stan-math.

There are three main directories. A UDF must be defined in the `functions` directory. There must be at least an example or a test that goes along with the UDF.

#### `functions`

Contains sub-directories by type of UDF. Only `.stan` files are contained in each directory. Since UDFs cannot have different signatures or type overloading, each version is contained in the same `.stan` file with a clarifying name (e.g. x_fun and x_fun_vectorized).

If you have a UDF that does not fit in the classification scheme please create an issue for discussion. 

#### `examples`

Contains a directory for the `stan` file example that calls the UDF function in the functions directory and the language (e.g. R) to run the model. 

A sufficient example will simulate data, pass it to the `stan` example, and have a minimal extraction to compare to the simulated data generating process.

#### `test`

Contains a directory for the `stan` file example that calls the UDF function in the functions directory and the language (e.g. R) to run the test. 

Tests are simpler than examples and do not contain data. The program is run with stans `fixed_param` algorithm. The `void` functions such as `expect_equal` will be called to compare to the expected output.

