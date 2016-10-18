#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


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
  double* specular;
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
  fprintf(fp, "#This document was converted from json to P6 by my converter\n", image.format);
  fprintf(fp, "%d %d\n%d\n", image.width, image.height, image.range);
  fwrite(image.buffer, sizeof(RGBpixel), image.height * image.width, fp);
  fclose(fp);
  return 0;
}

//run the variables through the quadratic formula and checks the t values.
double quadratic(double a, double b, double c){
   double t0 = (-b-sqrt((pow(b, 2))-4*a*c))/2*a;
   double t1 = (-b+sqrt((pow(b, 2))-4*a*c))/2*a;

   if(t0<0) return t1;
   if(t0<t1) return t0;
   return -1;
}

double sphere_intersection(double x_pos, double y_pos, double z_pos, int k){
  double a, b, c;
  a=pow(x_pos, 2)+pow(y_pos, 2)+pow(z_pos, 2);
  b=2 * (x_pos * (0-shapes[k].position[0]) + y_pos * (0-shapes[k].position[1]) + z_pos * (0-shapes[k].position[2]));
  c=(pow(0-shapes[k].position[0], 2) + pow(0-shapes[k].position[1], 2) + pow(0-shapes[k].position[2], 2)-pow(shapes[k].radius, 2));

  return quadratic(a, b, c);
}

//Plane intersection formula as well as setting up the variables.
///returns the t value so long as the denominator isn't zero.
double plane_intersection(double x_pos, double y_pos, double z_pos, int k){
  double a, b, c, d, mag, t;

  a = shapes[k].normal[0];
  b = shapes[k].normal[1];
  c = shapes[k].normal[2];
  mag = sqrt(pow(a, 2)+pow(b, 2)+pow(c, 2));
  a = a/mag;
  b = b/mag;
  c = c/mag;
  d = -(a*shapes[k].position[0] + b*shapes[k].position[1] + c*shapes[k].position[2]);
  if(a*x_pos + b*y_pos + c*z_pos == 0){
    return -1;
  }	    
  t = -d/(a*x_pos + b*y_pos + c*z_pos);
  return t;
}

double distance(double x, double y, double z){
  return sqrt(pow(x, 2)+pow(y, 2)+pow(z, 2));
}

int shader(int k, double x_pos, double y_pos, double z_pos, double t_min){
  double rOn_x, rOn_y, rOn_z, rDn_x, rDn_y, rDn_z, light_d, t, min_t;
  
  double* color = malloc(sizeof(double)*3);
  color[0] = 0;//ambient_color[0]
  color[1] = 0;//ambient_color[1]
  color[2] = 0;//ambient_color[2]

  for(int j=0, shapes[j] != NULL; j+=1){
    if(shapes[j].type != 3){
      printf("Only looking for lights");
    }else{
      rOn_x = min_t * x_pos;
      rOn_y = min_t * y_pos;
      rOn_z = min_t * z_pos;

      rDn_x = shapes[j].position[0] - rOn_x;
      rDn_y = shapes[j].position[1] - rOn_y;
      rDn_z = shapes[j].position[2] - rOn_z;

      light_d = distance(rDn_x, rDn_y, rDn_z);

      mag = sqrt(pow(rDn_x, 2)+pow(rDn_y, 2)+pow(rDn_z, 2));
      rDn_x = rDn_x/mag;
      rDn_y = rDn_y/mag;
      rDn_z = rDn_z/mag;

      min_t = light_d;

      for(int i=0; shapes[i] != NULL; i+=1){
	if(i == k) continue;
	
	if(shapes[i].type == 1){	  
	  
	  //Get some help with quadratic functions
	  t=sphere_intersection(rDn_x, rDn_y, rDn_z,  i);
	  //If t is positive and closer than t_min, lets paint.
	  if(t>0){
	    if(t_min == -1 || t<t_min){
	      t_min = t;
	      //Take the double value. Mult by max color value(255)
	      //Cast to unsigned character. (u_char) number
	      
	      /*
		image.buffer[i*width+j].r=(int)(shapes[k].color[0]*image.range);
		image.buffer[i*width+j].g=(int)(shapes[k].color[1]*image.range);
		image.buffer[i*width+j].b=(int)(shapes[k].color[2]*image.range);
	      */
	    }
	  }
	  //If its a plane..
	}else if(shapes[i].type == 2){
	  //Send pertinent to helper function and get the t-value
	  t = plane_intersection(rDn_x, rDn_y, rDn_z,  i);
	  //If t is positive and closer than t_min, lets paint.
	  if(t>0){
	    if(t_min == -1 || t<t_min){ 
	      t_min = t;
	      //Take the double value. Mult by max color value(255)
	      //Cast to unsigned character. (u_char) number

	      /*
	      image.buffer[i*width+j].r=(int)(shapes[k].color[0]*image.range);
	      image.buffer[i*width+j].g=(int)(shapes[k].color[1]*image.range);
	      image.buffer[i*width+j].b=(int)(shapes[k].color[2]*image.range);
	      */
	    }
	  }
	}else if(shapes[i].type == "light"){
	  printf("light\n");
	}else{
	  fprintf(stderr, "I'm not sure what shape that is!");
	}
	//Check for t, if its less than distance to light
	//Mark for shadow
      } 
    }   
  }
}
int caster(){
  //Instantiate the image object, and the buffer inside of it.
  image.buffer = malloc(sizeof(RGBpixel)*width*height);
  image.height = height;
  image.width = width;
  image.range = 255;
  //All of the math variables
  double x_pos, y_pos, z_pos, x_dist, y_dist, z_dist, distance, a, b, c, d, t0, t1;
  double mag,t=0;
  double t_min=-1.0;
  int i, j, k;
  //Nested for loop to iterate through all the pixels
  for (i=0;i<height;i++){
    for(j=0;j<width;j++){
      t_min = -1;
      //Pixel center
      x_pos = 0 - camerawidth/2+camerawidth/width*(j+0.5);
      y_pos = 0 - cameraheight/2+cameraheight/height*(i+0.5);
      //Normalize direction vector
      mag = sqrt(pow(x_pos, 2)+pow(y_pos, 2)+1);
      x_pos = x_pos/mag;
      y_pos = y_pos/mag;
      z_pos = 1/mag;
      //Loop through all objects
      for(k=0; k<oCount; k++){
	//If its a sphere, sphere intersection formula
 	if(shapes[k].type == 1){	  
	 
	  //Get some help with quadratic functions
	  t=sphere_intersection(x_pos, y_pos, z_pos, k);
	  //If t is positive and closer than t_min, lets paint.
	  if(t>0){
	    if(t_min == -1 || t<t_min){
	      t_min = t;
	      //Take the double value. Mult by max color value(255)
	      //Cast to unsigned character. (u_char) number

	      shading(k, x_pos, y_pos, z_pos, t_min);

	      /*
	      image.buffer[i*width+j].r=(int)(shapes[k].color[0]*image.range);
	      image.buffer[i*width+j].g=(int)(shapes[k].color[1]*image.range);
	      image.buffer[i*width+j].b=(int)(shapes[k].color[2]*image.range);
	      */
	    }
	  }
	  //If its a plane..
	}else if(shapes[k].type == 2){
	  //Send pertinent to helper function and get the t-value
	  t = plane_intersection(x_pos, y_pos, z_pos, k);
	  //If t is positive and closer than t_min, lets paint.
	  if(t>0){
	    if(t_min == -1 || t<t_min){ 
	      t_min = t;
	      //Take the double value. Mult by max color value(255)
	      //Cast to unsigned character. (u_char) number
	      image.buffer[i*width+j].r=(int)(shapes[k].color[0]*image.range);
	      image.buffer[i*width+j].g=(int)(shapes[k].color[1]*image.range);
	      image.buffer[i*width+j].b=(int)(shapes[k].color[2]*image.range);
	    }
	  }
	}else if(shapes[k].type == "light"){
	  printf("light\n");
	}else{
	  fprintf(stderr, "I'm not sure what shape that is!");
	}
      }
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
      shapes[oCount+1] = NULL;
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
      } else if (strcmp(value, "plane") == 0) {
	shapes[oCount].type = 2;
      } else if (strcmp(value, "light") == 0){
	shapes[oCount].type = 3;
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
  caster();
  write_p6(argv[4]);
  free(image.buffer);
  return 0;
}
