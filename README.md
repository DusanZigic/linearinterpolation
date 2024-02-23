linear interpolation class

this is a standalone class that performs linear interpolation for given data,
given on ordered grid; for example a function of 2 variables:
 x | y | f
---+---+---
 1 | 1 | 2
---+---+---
 1 | 2 | 3
---+---+---
 2 | 1 | 3
---+---+---
 2 | 2 | 4

constructors are overloaded so that they can take pointers to arrays or std::vectors
as function arguments;
constructors are overloaded so that they can take different number of inputs given
the number of variables of a function; for example for 1 variable function:
    interpolationF(const double *xData, const double *fData, size_t NofElements) or
    interpolationF(const std::vector<double> &xData, const std::vector<double> &fData)
for multiple variable function, pattern is the same, just with more inputs;
there is also SetData function which can be used as "constructor" if class object needs
to be defined before the variables and function values are set - for example in a
different function;
all data is copied inside private class variables, so there is no need to save it after
constructing class object;
there is an additional constructor for 2 variable function that does not follow pattern
of the rest of them - it takes 3 std::vectors where the first 2 variables grid points,
and the third one is 2D std::vector which are function values;

for acctual interpolation there is "interpolation" class method that takes variable number
inputs of type double and returns interpolated value;

there are additional domain and codomain class methods that return domain of interpolated
function as 2D std::vector and 1D std::vector respectively; domain return 2D std::vector
even is the function is of 1 varible in which case the dimensions of domain vector would
be (1, 2);