# cs430-illumination

The project was to create a raycaster that can eat in a .json file and turn it into a ppm.
I completed this by creating a raycaster designed to iterate through every pixel of a proposed picture, 
and together with camera height and pixel dimentions, translate the real world position in the .json file 
and make a 2d picture representing those coordinates.

The raycaster now accounts for lights as well as spheres and planes. With these equations, light is implemented as well as shadow.

Usage:
Compile with "make"
./raycast <width> <height> <input.json> <output.json>