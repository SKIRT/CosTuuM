#! /bin/sh -l

# Action script. This is the shell script that gets executed inside the
# Docker container.
# The Docker container has all the relevant libraries and compilers installed.

# Go into the directory where the repository was checked out. The action
# bot mounts this directory as a volume.
cd "$GITHUB_WORKSPACE"

# Create a build folder and enter it
mkdir build
cd build

# Create a flag variable to store the exit code
exitcode=0

# Configure the code using CMake
cmake ..
if [ $? -ne 0 ] ; then
  echo "Configuration failed!"
  exitcode=1
fi

# Now compile the code...
make
if [ $? -ne 0 ] ; then
  echo "Compilation failed!"
  exitcode=1
fi
# ...compile and run the unit tests...
make check
if [ $? -ne 0 ] ; then
  echo "Unit tests failed!"
  exitcode=1
fi
# ...install the Python library...
make install
if [ $? -ne 0 ] ; then
  echo "Installation failed!"
  exitcode=1
fi
# ...and check that we can actually load it from Python
python3 -c "import CosTuuM"
if [ $? -ne 0 ] ; then
  echo "Python module cannot be imported!"
  exitcode=1
fi

# Check if all tests passed.
if [ $exitcode -ne 0 ] ; then
  echo "One or more tests failed!"
  exit 1
else
  echo "All tests passed!"
fi
