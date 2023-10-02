# Description

This is a fork of the original cgmlst-dists project by Torsten Seemann (https://github.com/tseemann/cgmlst-dists).

This version allows to manage large input files and bypasses the Integer Overflow issue when the calculation of the final memory size of the distance vector is greater than MAX_INT.

# Compile
```console
make
```

# Check
```console
make check
```
# Example of usage

```console
./cgmlst-dists test/100.tab > 100output.tab
```
