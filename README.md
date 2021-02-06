# Unoffical Stan UDF Repo

A home for all the UDFs that the Stan community uses or finds useful that haven't made their way into Stan-math.

There are two main directories
* `functions`
    Contains the UDFs each in their own `.stan` file. Since UDFs cannot have different signatures or type overloading, each version is contained in the same `.stan` file with a clarifying name (e.g. x_fun and x_fun_vectorized).
* `examples`
    Contains two directories, one for `R` (or any other language) to call the example model and one for the example model in `Stan`. 
