<h1><img src="logo/interpolation.png" alt="logo" width='215' align="right"/> Linear Interpolation Class</h1>

**LinearInterpolator** is a high-performance, single-header C++ library for linear interpolation on ordered grids. It is designed for speed-critical applications (such as physics simulations) and supports up to 4D interpolation.

## Features
* **Single-Header**: Just `#include "linearinterpolator.hpp"` and you are ready to go.
* **Performance Optimized**: Internal logic is branchless-friendly, uses precomputed boundaries, and supports compiler vectorization.
* **Hybrid Error Handling**: Fast "poisoned index" crashes in Release mode; detailed diagnostic reports in Debug mode.
* **Fully Templated**: Supports `float`, `double`, and `long double`.

## Data Format

For a function of 2 variables, *f(x,y)=x+y*, data should be in form:
 x | y | f
--- | --- | ---
 1 | 1 | 2
 1 | 2 | 3
 2 | 1 | 3
 2 | 2 | 4

## Usage

### Initialization

The LinearInterpolator uses `std::vector` for data input and includes an optional name parameter. This name is used during **Debug Mode** to identify which specific object is causing an extrapolation error.

### Example: 1D Interpolation

```c++
// Constructor: (x_vector, f_vector, "OptionalName")
LinearInterpolator<float> intFunction(x, f, "intFunctionA");
```

### Example: 4D Interpolation

```c++
// Constructor: (x1, x2, x3, x4, f, "OptionalName")
LinearInterpolator<double> photonFlux(energy, angle, temp, density, flux, "PhotonFlux");
```

### Deferred Initialization

If you need to define the object before the data is ready, you can use the default constructor and the `setData` method later.

```c++
void loadFunction(LinearInterpolator<double> &intFunction) {
     // ... prepare vectors ...
     intFunction.setData(x, y, f, "intFunctionName"); // Name can be set here too
}

int main() {
     LinearInterpolator<double> intFunction;
     loadFunction(intFunction);
}
```

### Interpolation

```c++
// For a 2D-initialized object:
auto val = intFunction.interpolation(2.5, 1.2);
```

## Error Handling & Performance Modes

The library uses a preprocessor-driven safety strategy. Release mode is the default.

### 1. Release Mode (Default)
In this mode, the library is optimized for maximum throughput (30x speedup compared to standard exception-based code).

- Behavior: If an input is out of bounds (extrapolation), the function returns a "poisoned index" (size_t -1), which triggers a hardware-level SegFault upon memory access.

- Compilation flags: -O3 -march=native

### 2. Debug Mode
If you encounter a crash and need to find the specific variable causing the issue, recompile with the DEBUG flag.

- Behavior: Prints a detailed report to stderr (including Object Name, Dimension, Value, and Allowed Range) and throws a std::out_of_range exception.

- Compilation flags: -O3 -march=native -DDEBUG

[!IMPORTANT]
Extrapolation is not supported. Ensure your input values fall within the domain provided during initialization.

## Technical Note
This class performs internal copies of all data. Once the object is constructed or `setData` is called, you may safely delete or modify your original data structures.