C++ T-Matrix code, based on the T-matrix code by Mishchenko 
(https://www.giss.nasa.gov/staff/mmishchenko/t_matrix.html).

The code can be configured using CMake:

```
cmake <PATH TO SOURCE FOLDER>
```

we recommend an out-of-source configuration (where the build folder is 
either a subdirectory of the source folder or a folder that is outside 
the source folder). Once the configuration has succeeded, the code can 
be compiled using Make:

```
make -j <NUMBER OF THREADS TO USE>
```

The CosTuuM Python library can be installed for the current user using

```
make install
```

This will make it available from within Python, simply

```
import CosTuuM
```
