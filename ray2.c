// #include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct vector{
  double x, y, z;
} vector;

typedef struct ray{
  vector o, d;
} ray;

typedef struct sphere{
  vector c;
  double r;
  int material;
} sphere;

typedef struct material{
	colour diffuse;
	float reflection;
}material;

typedef struct colour {
  double red, green, blue;
}colour;

typedef struct light{
  vector pos;
  colour intensity;
} light;


void add_colour(colour *c, double r, double g, double b){
  c->red = r;
  c->green = g;
  c->blue = b;
}
// Atribiui is valores para o vetor passado.
void add_vec(vector *v, double x, double y, double z){
  v->x = x;
  v->y = y;
  v->z = z;
}
//Soma os dois vetores
vector * sum_vec(vector *v, vector *w){
  vector *u;
  u = malloc(sizeof(vector));
  u->x = v->x + w->x;
  u->y = v->y + w->y;
  u->z = v->z + w->z;
  return u;
}

//Subtrai os dois vetores
vector * sub_vec(vector *v, vector *w){
  vector *u;
  u = malloc(sizeof(vector));
  u->x = v->x - w->x;
  u->y = v->y - w->y;
  u->z = v->z - w->z;
  return u;
}

//Multiplica o vetor por um real
vector * mult_vec(vector *v, double w){
  vector *u;
  u = malloc(sizeof(vector));
  u->x = v->x * w;
  u->y = v->y * w;
  u->z = v->z * w;
  return u;
}

//Divide o vetor por um real
vector * div_vec(vector *v, double w){
  vector *u;
  u = malloc(sizeof(vector));
  u->x = v->x / w;
  u->y = v->y / w;
  u->z = v->z / w;
  return u;
}

void normalize(vector *v) {
  double mg = sqrt((v->x)*(v->x) + (v->y)*(v->y) + (v->z)*(v->z));
  v->x = v->x/mg;
  v->y = v->y/mg;
  v->z = v->z/mg;
}

void copy_vec(vector *orig, vector *dest){
  dest->x = orig->x;
  dest->y = orig->y;
  dest->z = orig->z;
}



double dot(vector *v, vector *w) {
  return (v->x*w->x + v->y*w->y + v->z*w->z);
}


void add_ray(ray *r, double o_x, double o_y, double o_z, double d_x, double d_y, double d_z){
  add_vec(&(r->o), o_x, o_y, o_z);
  add_vec(&(r->d), d_x, d_y, d_z);
}



void add_sphere(sphere *s, double x, double y, double z, double r){
  add_vec(&(s->c), x, y, z);
  s->r = r;
}

vector * get_normal_sphere(sphere *s, vector *pi) {
  vector *aux;
  aux = sub_vec(pi, &(s->c));
  vector *aux2;
  aux2 = div_vec(aux, s->r);
  free (aux);
  return aux2;
}

int intersect(sphere *s, ray *r, double *t) {
  vector *o;
  o = malloc(sizeof(vector));
  add_vec(o, (r->o).x, (r->o).y, (r->o).z);
  vector *d;
  d = malloc(sizeof(vector));
  add_vec(d, r->d.x, r->d.y, r->d.z);
  vector *oc;
  oc = sub_vec(o, &(s->c));
  double b = 2 * dot(oc, d);
  double c = dot(oc, oc) - (s->r)*(s->r);
  double disc = b*b - 4 * c;
  free (o);
  free(d);
  if (disc < 1e-4) return 0;
  disc = sqrt(disc);
  double t0 = -b - disc;
  double t1 = -b + disc;
  if(t0 > t1)
    t0 = t1;

  if((t0 > 0.001f) && (t0 < *t)){
    *t = t0;
    retval = 1;
  }else
    retval = 0;
}

void clamp255(vector *v) {
  v->x = (v->x > 255) ? 255 : (v->x < 0) ? 0 : v->x;
  v->y = (v->y > 255) ? 255 : (v->y < 0) ? 0 : v->y;
  v->z = (v->z > 255) ? 255 : (v->z < 0) ? 0 : v->z;
}

void saveppm(char *filename, unsigned char *img, int width, int height){
	/* FILE pointer */
	FILE *f;

	/* Open file for writing */
	f = fopen(filename, "wb");

	/* PPM header info, including the size of the image */
	fprintf(f, "P6 %d %d %d\n", width, height, 255);

	/* Write the image data to the file - remember 3 byte per pixel */
	fwrite(img, 3, width*height, f);

	/* Make sure you close the file */
	fclose(f);
}

int main() {

  const int HEIGHT = 500;
  const int WIDTH = 500;
  unsigned char img[3*WIDTH*HEIGHT];
  int hit;

  // vector *white;
  // white = malloc(sizeof(vector));
  // add_vec(white, 255, 255, 255);
  //
  // vector *black;
  // black = malloc(sizeof(vector));
  // add_vec(black, 0, 0, 0);
  //
  // vector *red;
  // red = malloc(sizeof(vector));
  // add_vec(red, 255, 0, 0);

  sphere s;
  add_sphere(&s, WIDTH*0.5, HEIGHT*0.5, 50, 50);

  // sphere light;
  // add_sphere(&light, 0, 0, 50, 1);

  // vector *pix_col;
  // pix_col = malloc(sizeof(vector));
  colour *rgb;
  rgb = malloc(sizeof(colour));
  ray *r;
  r = malloc(sizeof(ray));
  for (int y = 0; y < HEIGHT; ++y) {
    for (int x = 0; x < WIDTH; ++x) {
      // copy_vec(black, pix_col);
      add_colour(rgb, 0, 0, 0);
      int level = 0;
      double coef = 1.0;

      add_ray(r,x,y,0,0,0,1);

      do{
        double t = 20000.0f;
        int currentSphere = -1;
        //procura a esfera mais proxima desse pixel
        for(int i = 0; i < 3; i++){
          if(intersect(&s, r, &t))
            currentSphere = i;
        }
        if (currentSphere == -1) break;



      }while




      hit = intersect(&s, r, &t);
      if(hit) {
        img[(x + y*WIDTH)*3 + 0] = 255;
				img[(x + y*WIDTH)*3 + 1] = 0;
				img[(x + y*WIDTH)*3 + 2] = 0;
			}else{
				img[(x + y*WIDTH)*3 + 0] = 0;
				img[(x + y*WIDTH)*3 + 1] = 0;
				img[(x + y*WIDTH)*3 + 2] = 0;
			}
      // if (intersect(&s, r, &t)) {
      //   vector *pi;
      //   vector *aux;
      //   aux = mult_vec(&(r->d), t);
      //   pi = sum_vec(&(r->o), aux);
      //   vector *L;
      //   L = sub_vec(&(light.c), pi);
      //   vector *N;
      //   N = get_normal_sphere(&s, pi);
      //   normalize(L);
      //   normalize(N);
      //   double dt = dot(L, N);
      //
      //   // pix_col = (red + white*dt) * 0.5;
      //   vector *aux2;
      //   aux2 = mult_vec(white, dt);
      //   vector *aux3;
      //   aux3 = sum_vec(red, aux2);
      //   vector *aux4;
      //   aux4 = mult_vec(aux3, 0.5);
      //   clamp255(aux4);
      //   copy_vec(aux4, pix_col);
      //   free(pi);
      //   free(aux);
      //   free(L);
      //   free(N);
      //   free(aux2);
      //   free(aux3);
      //   free(aux4);
      // }
      // fprintf(fp, "%d %d %d\n", (int)pix_col->x, (int)pix_col->y, (int)pix_col->z);
    }
  }
  saveppm("image.ppm", img, WIDTH, HEIGHT);
}
