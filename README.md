
## Title: inside-crack

Modified from Ellad B. Tadmor and Ronald E. Miller' friction example

Brief Description: a crack in a FCC Ni single crystal

## Description: 

This example want to insert an inner crack or internal boundary conditon in QC. The Ni is modeled used the Universal 3
form of the Foiles, Baskes and Daw EAM potential. The crytal is intended to be a rectangule, turned out to be an abnormal geometry with an error `Point outside all grains` when running the code. Though the code insert the inner crack successfully, the error and its frustrating shape make the code futile. 

[A similar question](https://groups.google.com/forum/#!topic/qcmethod/QNEJ3tRutHM) asked by another person in QC forum is also wating for answer. I wonder if anybody can help me solve this problem.

## Directory contents:

```
inside-crack
      |
      |-- Makefile..........Make file for building
      |
      |-- README.md............This file
      |
      |-- user_fric.f.......User file for this inside-crack problem 
      |
      |--/Shear
            |
            |-- fric_shear.in.....Input file for the problem
            |
            |-- friction.geo......Grain definition file for this problem
```

## Installation and execution instructions:

1. Place the inside-crack directory under $QC (where $QC is the location where the QC distribution is located).

2. To make the executable:
```bash
cd $QC/Friction-example
make QCCOMPILER=COMP
```

 (Here `COMP` is the desired Fortran compiler. See QC Tutorial guide.)

3. To run inside-crack:
```bash
cd $QC/Friction-example/Shear
../fric < fric_shear.in > fric_shear.out
```
