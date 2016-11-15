#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NS 20
#define PI 3.14159265
#define CAMPOS {0, 0, 0}
#define EPS .01

typedef struct RGBpixel {
  unsigned char r, g, b;
} RGBpixel;

typedef struct Image{
  int width, height, format, range;
  RGBpixel *buffer;
}Image;

typedef struct {
  int type;//1 - sphere, 2 - plane 3 - light
  double* color;
  double radius;
  double* position;
  double* normal;
  double* direction;
  double radialA0;//default 0
  double radialA1;//default 0
  double radialA2;//default 1
  double angularA0;
  double* diffuse;
  double reflectivity;
  double refractivity;
  double ior;
  double* specular;
  double theta;
}Shape;

Shape shapes[180];
Image image;

int height;
int width;
int oCount=0;
double cameraheight;
double camerawidth;

int line = 1;

int write_p6(char* input){
  FILE* fp = fopen(input, "wb");
  //create and print to a p6 file.
  fprintf(fp, "P6\n");
  fprintf(fp, "#This document was converted from json to P6 by my converter\n");
  fprintf(fp, "%d %d\n%d\n", image.width, image.height, image.range);
  fwrite(image.buffer, sizeof(RGBpixel), image.height * image.width, fp);
  fclose(fp);
  return 0;
}

static inline double sqr(double v){
  return v*v;
}

static inline double distance(double* v){
  return sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
}

static inline void normalize(double* v){
  double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

static inline double dot(double* v, double* b){
  return (v[0]*b[0] + v[1]*b[1] + v[2]*b[2]);
}

static inline void vadd(double* v, double* b){
  v[0] += b[0];
  v[1] += b[1];
  v[2] += b[2];
}

void print_scene(){
    for(int i=0; i<oCount; i++){
        if(shapes[i].type != 3) {
            printf("----------\n");
            printf("Type: %d\n", shapes[i].type);
            printf("Color: [%f %f %f]\n", shapes[i].diffuse[0], shapes[i].diffuse[1], shapes[i].diffuse[2]);
            printf("Position: [%f %f %f]\n", shapes[i].position[0], shapes[i].position[1], shapes[i].position[2]);
        }
    }
}


//k - intersecting object
//j - light

void diffuse_color(double* v, int k, int j, double alpha){
  //find the Diffuse color
  ///if Alpha is negative or zero, everything gets set to 0,
  ///else multiply diffuse color of the object, the color of the light and the alpha
  if(alpha <= 0 || shapes[k].diffuse == NULL){
    v[0] = 0;
    v[1] = 0;
    v[2] = 0;
  }else{
    v[0] = shapes[k].diffuse[0] * shapes[j].color[0] * alpha;
    v[1] = shapes[k].diffuse[1] * shapes[j].color[1] * alpha;
    v[2] = shapes[k].diffuse[2] * shapes[j].color[2] * alpha;
  }
}

void specular_color(double* v, int k, int j, double alpha){
  //Find the specular color
  ///if alpha is negative or zero, specular is zero,
  ///else multiply the specular color, light color, and alpha sqared to the NS (constant)
  if(alpha <= 0 || shapes[k].specular == NULL){
    v[0] = 0;
    v[1] = 0;
    v[2] = 0;
  }else{
    alpha = pow(alpha, NS);
    v[0] = shapes[k].specular[0] * shapes[j].color[0] * alpha;
    v[1] = shapes[k].specular[1] * shapes[j].color[1] * alpha;
    v[2] = shapes[k].specular[2] * shapes[j].color[2] * alpha;
  }
}

//run the variables through the quadratic formula and checks the t values.
double quadratic(double b, double c){
  if((sqr(b) - 4 * c) < 0){
    return -1;
  }
  double t0 = (-b-sqrt((pow(b, 2))-4*c))/2;
  double t1 = (-b+sqrt((pow(b, 2))-4*c))/2;

  //printf("Quad a - %f, b - %f, c - %f\n", a, b, c);
  //if(t0 > 0 || t1 > 0) printf("t0 - %f, t1 - %f\n", t0, t1);
  if(t0<0) return t1;
  if(t0<t1) return t0;
  return -1;
}
double sphere_intersection(double* Rn, double* Rd, int k){
  double a, b, c;

  b=2 * (Rd[0] * (Rn[0]-shapes[k].position[0]) + Rd[1] * (Rn[1]-shapes[k].position[1]) + Rd[2] * (Rn[2]-shapes[k].position[2]));
  c=(sqr(Rn[0]-shapes[k].position[0]) + sqr(Rn[1]-shapes[k].position[1]) + sqr(Rn[2]-shapes[k].position[2])-sqr(shapes[k].radius));

  return quadratic(b, c);
}

//Plane intersection formula as well as setting up the variables.
///returns the t value so long as the denominator isn't zero.
double plane_intersection(double* Rd, int k){
  double d, t;
  double* nPos;

  nPos = shapes[k].normal;
  normalize(nPos);

  d = -(nPos[0]*shapes[k].position[0] + nPos[1]*shapes[k].position[1] + nPos[2]*shapes[k].position[2]);
  if(nPos[0]*Rd[0] + nPos[1]*Rd[1] + nPos[2]*Rd[2] == 0){
    return -1;
  }
  t = -d/(nPos[0]*Rd[0] + nPos[1]*Rd[1] + nPos[2]*Rd[2]);
  return t;
}

double frad(int k, double d){
  //inputs distance to light
  ///if distance is infinity, return 1, else return 1/(a2*dl^2 + a1*dl + a0)
  return 1/(shapes[k].radialA2 * sqr(d) + shapes[k].radialA1 * d + shapes[k].radialA0);
}

double fang(double* rDn, int j){
  //1 if not spotlight(theta = 0r no direction)
  //0 if rLt dot shapes[k].direction = cos(alpha) and is less than cos(theta)
  //else rLt dot shapes[k].direction^a0
  if(shapes[j].theta == 0 || shapes[j].direction == NULL){
    return 1;
  }
  double v = dot(rDn, shapes[j].direction);
  if(v < cos(shapes[j].theta)){
    return 0;
  }else{
    return pow(v, shapes[j].angularA0);
  }
}


//Inputs - index position, object index, unit vector pointing towards point of intersection and length of intersection.
///caculates illumination due to lights.
double* shader(double* Ro, double* Rd, int level){
    double light_d, t=0, alpha, fRad, fAng;
    double t_min = -1;
    int i, j, k;
    int closest_object = -1;
    double N[3];
    double rOn[3];

    //set my color offset
    double color[3];
    color[0] = .1;//ambient_color[0]
    color[1] = .1;//ambient_color[1]
    color[2] = .1;//ambient_color[2]

      if(level > 6){
        return color;
    }

      for(k=0; k<oCount; k++){
        //If its a sphere, sphere intersection formula
        if(shapes[k].type == 1){

          //Get some help with quadratic functions
          t=sphere_intersection(Ro, Rd, k);
        //If t is positive and closer than t_min, lets paint.
        if(t>0){
            if(t_min == -1 || t<t_min){
                t_min = t;
                //set t_min and closest_object for shading
                closest_object = k;
                }
            }
        //If its a plane..
        }else if(shapes[k].type == 2){
         //Send pertinent to helper function and get the t-value
            t = plane_intersection(Rd, k);
            //If t is positive and closer than t_min, lets paint.

            if(t>0){
                if(t_min == -1 || t<t_min){
                    t_min = t;
                    closest_object = k;
                    //set t_min and closest object for shading
                }
            }
        }else if(shapes[k].type == 3){

        }else{
            fprintf(stderr, "I'm not sure what shape that is!");
        }
      }

        //printf("%d\n", k);
    if(closest_object != -1 && t_min != -1){



        rOn[0] = t_min * Rd[0] + Ro[0];
        rOn[1] = t_min * Rd[1] + Ro[1];
        rOn[2] = t_min * Rd[2] + Ro[2];

        if(shapes[closest_object].type == 1){
      //Sphere
            N[0] = rOn[0] - shapes[closest_object].position[0];
            N[1] = rOn[1] - shapes[closest_object].position[1];
            N[2] = rOn[2] - shapes[closest_object].position[2];
        }else if(shapes[closest_object].type == 2){
      //Plane
            N[0] = shapes[closest_object].normal[0];
            N[1] = shapes[closest_object].normal[1];
            N[2] = shapes[closest_object].normal[2];
        }else{
            printf("What object am I hitting again?\n");
        }
        //loop through lights
        for(j=0; j<oCount; j++){
            if(shapes[j].type != 3){
                //Only looking for lights
            }else{
                //found a light, calculate the point on the object and a vector towards the light.
                double rDn [3];
                rDn[0] = shapes[j].position[0] - rOn[0];
                rDn[1] = shapes[j].position[1] - rOn[1];
                rDn[2] = shapes[j].position[2] - rOn[2];

                //store distance to light for t checking
                light_d = distance(rDn);

                normalize(rDn);

                t_min = light_d;

                 //loop through objects looking for one in between point and light.
                for(i=0; i<oCount; i++) {
                    //don't want to consider the object the point is on for shading.
                    t = -1;
                    if (i == closest_object) continue;

                    if (shapes[i].type == 1) {
                        //find the sphere intersection.
                        t = sphere_intersection(rOn, rDn, i);
                        //if t is positive and lower than the distance to light, its in the way.
                        if (t > 0) {
                            if (t < t_min) {
                                t_min = t;
                            }
                        }
                    } else if (shapes[i].type == 2) {
                        //Send pertinent to helper function and get the t-value
                        t = plane_intersection(rDn, i);
                        //If t is positive and closer than distance to light, its in the way.
                        if (t > 0) {
                            if (t < t_min) {
                                t_min = t;
                            }
                        }
                        // not looking for lights
                    } else if (shapes[i].type == 3) {

                    } else {
                        fprintf(stderr, "I'm not sure what shape that is!");
                    }
                }
  	          //No object in the way
                if(fabs(t_min - light_d) < 0.00000000001){
  	         //Calculate and update color
  	          //closest_object = intersecting object, j = light object

  	                    //Store L = unit vector from obj to light

                    double L[3];
                    L[0] = rDn[0];
                    L[1] = rDn[1];
                    L[2] = rDn[2];

                  //Store R the reflected vector from the light vector off the object
                    double R[3];
                    alpha = dot(rDn, N);
                    R[0] = rDn[0] - (2 * alpha * N[0]);
                    R[1] = rDn[1] - (2 * alpha * N[1]);
                    R[2] = rDn[2] - (2 * alpha * N[2]);

                    normalize(R);

                    //Store V, the vector from object to camera
                    double V[3];
                    V[0] = Rd[0];
                    V[1] = Rd[1];
                    V[2] = Rd[2];

                    //Not normal and Obj to light vector, for diffuse calculation
                    alpha = dot(N, L);


                    //Initialize diffuse color vector
                    double diff[3];
                    diff[0] = 0;
                    diff[1] = 0;
                    diff[2] = 0;

                    //calculate diffuse color
                    diffuse_color(diff, closest_object, j, alpha);


                    //printf("[%f %f %f]--[%f %f %f]\n", shapes[closest_object].color[0], shapes[closest_object].color[1], shapes[closest_object].color[2], diff[0], diff[1], diff[2]);
                    alpha = dot(V, R);
                    double spec[3];
                    spec[0] = 0;
                    spec[1] = 0;
                    spec[2] = 0;

                    //Calculate specular color
                    specular_color(spec, closest_object, j, alpha);

                    //Calculate radial light
                    fRad = frad(j, light_d);

                    //Calculate angular light
                    fAng = fang(rDn, j);



                    //Add to color, frad(j, light_d) * fang(rDn, j) * (spec + diff)
                    color[0] += fRad * fAng * (spec[0] + diff[0]);
                    color[1] += fRad * fAng * (spec[1] + diff[1]);
                    color[2] += fRad * fAng * (spec[2] + diff[2]);


                }else if(t_min < light_d){

                }else {
                    printf("Not sure what happened\n");
                }
            }
        }
        double nextColor[3];
        nextColor[0] = 0;
        nextColor[1] = 0;
        nextColor[2] = 0;

        if(shapes[closest_object].reflectivity != 0) {
            alpha = dot(Rd, N);
            double reflDirection[3];
            reflDirection[0] = Rd[0] - (2 * alpha * N[0]);
            reflDirection[1] = Rd[1] - (2 * alpha * N[1]);
            reflDirection[2] = Rd[2] - (2 * alpha * N[2]);
            normalize(reflDirection);
            double reflPoint[3];
            reflPoint[0] = rOn[0] + EPS * reflDirection[0];
            reflPoint[1] = rOn[1] + EPS * reflDirection[1];
            reflPoint[2] = rOn[2] + EPS * reflDirection[2];

            vadd(nextColor, shader(reflPoint, reflDirection, level + 1));
            nextColor[0] *= shapes[closest_object].reflectivity;
            nextColor[1] *= shapes[closest_object].reflectivity;
            nextColor[2] *= shapes[closest_object].reflectivity;
        }

        vadd(color, nextColor);
    }
    return color;
}


int caster(){
  //Instantiate the image object, and the buffer inside of it.
  image.buffer = malloc(sizeof(RGBpixel)*width*height);
  image.height = height;
  image.width = width;
  image.range = 255;
  //All of the math variables
  int i, j;

  //Camera Position -- Set to 0.
  double Rn[3];
  Rn[0] = 0;
  Rn[1] = 0;
  Rn[2] = 0;
  //Nested for loop to iterate through all the pixels
  for (i=0;i<height;i++){
    for(j=0;j<width;j++){
      //Pixel center
      double Rd [3];
      Rd[0] = Rn[0] - camerawidth/2+camerawidth/width*(j+0.5);
      Rd[1] = Rn[1] - cameraheight/2+cameraheight/height*(i+0.5);
      Rd[2] = 1;

      //Normalize direction vector
      normalize(Rd);
      //Loop through all objects
      double color[3];
      color[0] = 0;
      color[1] = 0;
      color[2] = 0;
      vadd(color, shader(Rn, Rd, 1));
      //printf("[%f %f %f]\n", color[0], color[1], color[2]);
      //Clamp color
      if(color[0] > 1){
        color[0] = 1;
      }
      if(color[1] > 1){
        color[1] = 1;
      }
      if(color[2] > 1){
        color[2] = 1;
      }
      image.buffer[i*width+j].r=(unsigned char)(color[0]*image.range);
      image.buffer[i*width+j].g=(unsigned char)(color[1]*image.range);
      image.buffer[i*width+j].b=(unsigned char)(color[2]*image.range);
    }
  }

  return 0;
}

// next_c() wraps the getc() function and provides error checking and line
// number maintenance
int next_c(FILE* json) {
  int c = fgetc(json);
#ifdef DEBUG
  printf("next_c: '%c'\n", c);
#endif
  if (c == '\n') {
    line += 1;
  }
  if (c == EOF) {
    fprintf(stderr, "Error: Unexpected end of file on line number %d.\n", line);
    exit(1);
  }
  return c;
}


// expect_c() checks that the next character is d.  If it is not it emits
// an error.
void expect_c(FILE* json, int d) {
  int c = next_c(json);
  if (c == d) return;
  fprintf(stderr, "Error: Expected '%c' on line %d.\n", d, line);
  exit(1);
}


// skip_ws() skips white space in the file.
void skip_ws(FILE* json) {
  int c = next_c(json);
  while (isspace(c)) {
    c = next_c(json);
  }
  ungetc(c, json);
}


// next_string() gets the next string from the file handle and emits an error
// if a string can not be obtained.
char* next_string(FILE* json) {
  char buffer[129];
  int c = next_c(json);
  if (c != '"') {
    fprintf(stderr, "Error: Expected string on line %d.\n", line);
    exit(1);
  }
  c = next_c(json);
  int i = 0;
  while (c != '"') {
    if (i >= 128) {
      fprintf(stderr, "Error: Strings longer than 128 characters in length are not supported.\n");
      exit(1);
    }
    if (c == '\\') {
      fprintf(stderr, "Error: Strings with escape codes are not supported.\n");
      exit(1);
    }
    if (c < 32 || c > 126) {
      fprintf(stderr, "Error: Strings may contain only ascii characters.\n");
      exit(1);
    }
    buffer[i] = c;
    i += 1;
    c = next_c(json);
  }
  buffer[i] = 0;
  return strdup(buffer);
}

double next_number(FILE* json) {
  double value;
  fscanf(json, "%lf", &value);
  // Error check this..
  return value;
}

double* next_vector(FILE* json) {
  double* v = malloc(3*sizeof(double));
  expect_c(json, '[');
  skip_ws(json);
  v[0] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[1] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[2] = next_number(json);
  skip_ws(json);
  expect_c(json, ']');
  return v;
}


void read_scene(char* filename) {
  int c;
  FILE* json = fopen(filename, "r");
  if (json == NULL) {
    fprintf(stderr, "Error: Could not open file \"%s\"\n", filename);
    exit(1);
  }

  skip_ws(json);

  // Find the beginning of the list
  expect_c(json, '[');

  skip_ws(json);

  // Find the objects
  while (1) {
    c = fgetc(json);
    if (c == ']') {
      fclose(json);
      return;
    }
    if (c == '{') {
      skip_ws(json);

      // Parse the object
      char* key = next_string(json);
      if (strcmp(key, "type") != 0) {
	fprintf(stderr, "Error: Expected \"type\" key on line number %d.\n", line);
	exit(1);
      }

      skip_ws(json);

      expect_c(json, ':');

      skip_ws(json);

      char* value = next_string(json);

      if (strcmp(value, "camera") == 0) {   //Set type, decrease object count for camera
	oCount = oCount - 1;
      } else if (strcmp(value, "sphere") == 0) {
	shapes[oCount].type = 1;
    shapes[oCount].ior = 1;
    shapes[oCount].reflectivity = 0;
    shapes[oCount].refractivity = 0;
      } else if (strcmp(value, "plane") == 0) {
	shapes[oCount].type = 2;
    shapes[oCount].ior = 1;
    shapes[oCount].reflectivity = 0;
    shapes[oCount].refractivity = 0;
      } else if (strcmp(value, "light") == 0){
	shapes[oCount].type = 3;
	shapes[oCount].angularA0 = 1;
      } else {
	fprintf(stderr, "Error: Unknown type, \"%s\", on line number %d.\n", value, line);
	exit(1);
      }

      skip_ws(json);

      while (1) {
	// , }
	c = next_c(json);
	if (c == '}') {
	  // stop parsing this object
	  oCount = oCount + 1;
	  break;
	} else if (c == ',') {
	  // read another field
	  skip_ws(json);
	  char* key = next_string(json);
	  skip_ws(json);
	  expect_c(json, ':');
	  skip_ws(json);
	  if (strcmp(key, "width") == 0){//store appropriate value into the shapes struct.
	    camerawidth = next_number(json);
	  }else if (strcmp(key, "height") == 0){
	    cameraheight = next_number(json);
	  }else if (strcmp(key, "radius") == 0){
	    shapes[oCount].radius = next_number(json);
	  }else if (strcmp(key, "color") == 0){
            shapes[oCount].color = next_vector(json);
	  }else if (strcmp(key, "position") == 0){ //For position and normal, flip the y coordinate so that positive Y means up.
	    shapes[oCount].position = next_vector(json);
	    shapes[oCount].position[1] = -1*shapes[oCount].position[1];
	  }else if (strcmp(key, "normal") == 0){
	    shapes[oCount].normal = next_vector(json);
	    shapes[oCount].normal[1] = -1*shapes[oCount].normal[1];
	  }else if (strcmp(key, "direction") == 0){
	    shapes[oCount].direction = next_vector(json);
	    shapes[oCount].direction[1] += -1;
	  }else if (strcmp(key, "radial-a2") == 0){
	    shapes[oCount].radialA2 = next_number(json);
	  }else if (strcmp(key, "radial-a1") == 0){
	    shapes[oCount].radialA1 = next_number(json);
	  }else if (strcmp(key, "radial-a0") == 0){
	    shapes[oCount].radialA0 = next_number(json);
	  }else if (strcmp(key, "angular-a0") == 0){
	    shapes[oCount].angularA0 = next_number(json);
	  }else if (strcmp(key, "diffuse_color") == 0){
	    shapes[oCount].diffuse = next_vector(json);
	  }else if (strcmp(key, "specular_color") == 0){
	    shapes[oCount].specular = next_vector(json);
	  }else if (strcmp(key, "theta") == 0){
	    shapes[oCount].theta = next_number(json);
    }else if (strcmp(key, "refractivity") == 0){
      shapes[oCount].refractivity = next_number(json);
    }else if (strcmp(key, "reflectivity") == 0){
      shapes[oCount].reflectivity = next_number(json);
    }else if (strcmp(key, "ior") == 0){
      shapes[oCount].ior = next_number(json);
	  }else{
	    fprintf(stderr, "Error: Unknown property, \"%s\", on line %d\n",
		    key, line);
	  }
	  skip_ws(json);
	} else {
	  fprintf(stderr, "Error: Unexpected value on line %d\n", line);
	  exit(1);
	}
      }
      skip_ws(json);
      c = next_c(json);
      if (c == ',') {
	// noop
	skip_ws(json);
      } else if (c == ']') {
	fclose(json);
	return;
      } else {
	fprintf(stderr, "Error: Expecting ',' or ']' on line %d.\n", line);
	exit(1);
      }
    }
  }
}


int main(int argc, char** argv) {

  if(argc != 5){
    fprintf(stderr, "Not the correct arguments: please add <width> <height> <input(.json)> <output file name(.ppm)>");
    exit(1);
  }
  width = atoi(argv[1]);
  height = atoi(argv[2]);
  read_scene(argv[3]);
    print_scene();
  caster();
  write_p6(argv[4]);
  free(image.buffer);
  return 0;
}
