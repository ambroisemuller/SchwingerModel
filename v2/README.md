# SchwingerModel

Lattice simulation for QED (quantum electrodynamics) in two dimensions (a.k.a. the Schwinger model).

*Based on the <u>Lattice Practices 2021</u> code by Mattia Dalla Brida* (link).<br> **Rewritten in C++ by Ambroise Muller.**

## To download

```
git clone ${link-to-repository}
cd SchwingerModel
```

## To setup for studies

Requirements: at least Python 3.0

```
mkdir results results/observables
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install numpy pandas matplotlib scipy
```

## To build

Requirements: C++ 17, CMake 3.10

On `euler.ethz.ch` use e.g.
```
env2lmod
set_software_stack.sh new
module load gcc/8.2.0 cmake/3.26.3
```
Create makefile and compile using
```
mkdir build
cd build
cmake ..
make
```

## To run

In `SchwingerModel/v2/build/` use

```
./main
```

## Output

*Coming soon.*