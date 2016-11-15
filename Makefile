all: raytracer.c
	gcc raytracer.c -o raytrace -lm

clean:
	rm -rf raytrace *~
	rm *.ppm
