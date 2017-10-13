
## Title: inside-crack

Modified from Ellad B. Tadmor and Ronald E. Miller' friction example

Brief Description: a crack in FCC Ni single crystal

## Description: 

This example is designed to insert an inner crack or internal boundary conditon in QC. The Ni is modeled used the Universal 3
form of the Foiles, Baskes and Daw EAM potential. To insert an inner crack,

1. we should apply different algorithms to outer boundaries and inner boundaries.

      - Dealing with the outer first, iterate each vertex, get outer boundary nodes with linear interpolation in counter-clockwise. This procedure is similar to finding boundaries for a simple rectangle.

      - Then with inner, iterate each vertex, get inner boundary nodes with linear interpolation in clockwise, and make the last node in `elist` point to the start point of the iteration. 

      - Most importantly, in `.geo` file we should list the vertexes in a special sequence, by which we can draw the outline of our geometry without intersections.  The start and end point of outer or inner should be the same, namely, there should be no less than 5 vertexes if the outer is a rectangular.

[A similar question](https://groups.google.com/forum/#!topic/qcmethod/QNEJ3tRutHM) asked by another person in QC forum is also waiting for answer. I think my code is a proper way to answer that question.

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
      `--/Shear
            |
            |-- fric_shear.in.....Input file for this problem
            |
            `-- friction.geo......Grain definition file for this problem
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
