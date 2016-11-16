# cs430-Raytracer

The project was to create a raytracer that can eat in a .json file and turn it into a ppm.
I completed this by creating a raycaster designed to iterate through every pixel of a proposed picture,
and together with camera height and pixel dimentions, translate the real world position in the .json file
and make a 2d picture representing those coordinates.

The raycaster now accounts for lights as well as spheres and planes. With these equations, light is implemented as well as shadow.

The evolution to raytracer happened when I added recursive functionality to calculate reflection of light off other objects.

***Refraction doesn't currently have functionality, as it has been bugging out on me. I am still working on this***

Usage:
Compile with "make"
Execution:
./raytrace <width> <height> <input.json> <output.ppm>
