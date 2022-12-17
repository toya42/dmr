# dmr

A numerical simulation code for a double mach reflection problem.

[![MIT License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)

dmr is a code of finite difference method for computational fluid dynamics, written in FORTRAN77. The code can simulate double mach reflection problem.

<p align="center">
     <img src="https://github.com/toya42/garage/blob/master/dmr/density_contour.jpeg"
width="954" height="401"
alt="double mach reflection"
title="density contour of double mach reflection problem (t=0.2)">
</p>

density contour (t=0.2) of double mach reflection problem

## Numerical method

### spatial discretization and interpolation

finite difference method with third order MUSCL interpolation

### Riemann solver

SLAU

### time integration

third order TVD Runge Kutta method

## Double Mach Reflection problem

The double mach reflection problem is firstly proposed by Woodward and Collela (JCP, 1984)[1].
This problem is an important test case for the assessment of the resolution of Euler codes[2].

<p align="center">
     <img src="https://github.com/toya42/garage/blob/master/dmr/density_contour.gif"
width="954" height="401"
alt="double mach reflection"
title="density contour movie of double mach reflection problem">
</p>

## Dependency

Only a fortran compiler is required to compile this code.

## Compile & Execute

### Compile

```shell
$ gfortran main.f
```

### Execute

```shell
$./a.out
```

## Visualize

### Output file format

`grid.xyz` is a two dimensional grid file.

`flowfield_?????.q` are solution files.

File formats of these files are [**plot3d**](https://www.grc.nasa.gov/www/wind/valid/plot3d.html).

### Visualize

Any visualization softwares can be used that supports plot3d format. For example, Paraview, Fieldveiw, etc...

## License

This software is release under the MIT Licese, see [LICENSE](https://github.com/toya42/dmr/blob/main/LICENSE)

## Contact

Please report bugs and other issues through the issue tracker at:

https://github.com/toya42/dmr/issues

## References

[1] Woodward and Collela, Journal of Computational Physics, 1984.
[2] Vevek et. al., Jornal of Scientific Computing, 2019.