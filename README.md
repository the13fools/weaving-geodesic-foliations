Code for the geodesic field design and integration steps.

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `relax-field` binary.

## Run

From within the build directory just issue:

    ./relax-field

Once inside the program, load an .obj and then press "Whole Pipeline." We apologize, we have not spend much time trying to optimize the code, so it may take 10-15 mins per example depending on your machine specs. "Advanced mode" contains *many* (user-unfriendly) additional options.
