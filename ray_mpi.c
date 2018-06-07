// #include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#ifndef NMAQUINAS
#define NMAQUINAS 16
#endif



#define min(a,b) (((a) < (b)) ? (a) : (b))

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

typedef struct colour {
  double red, green, blue;
}colour;

typedef struct material{
	colour diffuse;
	double reflection;
}material;


typedef struct light{
  vector pos;
  colour intensity;
} light;




void add_light(light *l, double x, double y, double z, double r, double g, double b){
  l->pos.x = x;
  l->pos.y = y;
  l->pos.z = z;
  l->intensity.red = r;
  l->intensity.green = g;
  l->intensity.blue = b;

}

void add_material(material *m, double r, double g, double b, double reflec){
  m->diffuse.red = r;
  m->diffuse.green = g;
  m->diffuse.blue = b;
  m->reflection = reflec;
}

void add_colour(colour *c, double r, double g, double b){
  c->red = r;
  c->green = g;
  c->blue = b;
}
// Atribiui is valores para o vetor passado.
void inic_vec(vector *v, double x, double y, double z){
  v->x = x;
  v->y = y;
  v->z = z;
}
//Soma os dois vetores
vector soma_vec(vector *v, vector *w){
  vector u;

  u.x = v->x + w->x;
  u.y = v->y + w->y;
  u.z = v->z + w->z;
  return u;
}

//Subtrai os dois vetores
vector sub_vec(vector *v, vector *w){
  vector u;
  u.x = v->x - w->x;
  u.y = v->y - w->y;
  u.z = v->z - w->z;
  return u;
}

//Multiplica o vetor por um real
vector mult_vec(vector *v, double w){
  vector u;
  u.x = v->x * w;
  u.y = v->y * w;
  u.z = v->z * w;
  return u;
}

//Divide o vetor por um real
vector div_vec(vector *v, double w){
  vector u;
  u.x = v->x / w;
  u.y = v->y / w;
  u.z = v->z / w;
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
  inic_vec(&(r->o), o_x, o_y, o_z);
  inic_vec(&(r->d), d_x, d_y, d_z);
}



void add_sphere(sphere *s, double x, double y, double z, double r, int m){
  inic_vec(&(s->c), x, y, z);
  s->r = r;
  s->material = m;
}

vector get_normal_sphere(sphere *s, vector *pi) {
  vector aux;
  aux = sub_vec(pi, &(s->c));
  vector aux2;
  aux2 = div_vec(&aux, s->r);
  return aux2;
}

int intersect(sphere *s, ray *r, double *t) {

  vector o;
  inic_vec(&o, (r->o).x, (r->o).y, (r->o).z);
  vector d;
  inic_vec(&d, r->d.x, r->d.y, r->d.z);
  vector oc;
  oc = sub_vec(&o, &(s->c));
  double b = 2 * dot(&oc, &d);
  double c = dot(&oc, &oc) - (s->r)*(s->r);
  double disc = b*b - 4 * c;
  if (disc < 1e-4) return 0;
  disc = sqrt(disc);
  double t0 = (-b - disc)/2;
  double t1 = (-b + disc)/2;
  if(t0 > t1)
    t0 = t1;

  if((t0 > 0.001f) && (t0 < *t)){
    *t = t0;
    return 1;
  }else
    return 0;
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

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  if (argc < 3){
    printf("./ray <variacao (1 ou 2)> <tamanho (1, 2, 4)>\n");
    return 0;
  }
  int HEIGHT = 4500; // 4500 6188 9000
  int WIDTH = 8000;
  if(argv[2][0] == '1'){
       HEIGHT = 4480; // 4500 6188 9000
       WIDTH = 8000; // 8000 11000 16000
      // HEIGHT = 30; // 4500 6188 9000
      // WIDTH = 30; // 8000 11000 16000
  } else if(argv[2][0] == '2'){
      HEIGHT = 6176; // 4500 6188 9000
      WIDTH = 11000; // 8000 11000 16000
  } else if(argv[2][0] == '4'){
      HEIGHT = 8992; // 4500 6188 9000
      WIDTH = 16000; // 8000 11000 16000
  } else {
      printf("./ray <variacao (1 ou 2)> <tamanho (1, 2, 4)>\n");
      return 0;
  }
  
  int N_SPHERES = 18;
  unsigned long int tamanho_bloco = (3*WIDTH*HEIGHT) / NMAQUINAS;
  
  int N_LIGHTS = 6;
  int N_MATERIALS = 9;
  // unsigned char img[3*WIDTH*HEIGHT];

  int my_rank, n_procs;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

  unsigned char *img;
  img = malloc(sizeof(unsigned char)*tamanho_bloco);

  

  sphere *spheres;
  spheres = malloc(sizeof(sphere)*N_SPHERES);

  if(argv[1][0] == '1'){
        add_sphere(&spheres[0], 200, 300, 0, 100, 0);
        add_sphere(&spheres[1], 400, 400, 0, 100, 1);
        add_sphere(&spheres[2], 500, 140, 0, 100, 2);
        add_sphere(&spheres[3], 600, 100, 0, 100, 3);
        add_sphere(&spheres[4], 100, 800, 0, 100, 4);
        add_sphere(&spheres[5], 750, 900, 0, 100, 5);
        add_sphere(&spheres[6], 750, 100, 0, 100, 6);
        add_sphere(&spheres[7], 200, 900, 0, 100, 0);
        add_sphere(&spheres[8], 250, 550, 0, 100, 4);
        add_sphere(&spheres[10], 1500, 600, -30, 80, 6);
        add_sphere(&spheres[11], 1700, 800, 10, 130, 5);
        add_sphere(&spheres[12], 6000, 3000, 10, 500, 6);
        add_sphere(&spheres[13], 5000, 3000, 10, 500, 8);
        add_sphere(&spheres[14], 3000, 3000, 10, 500, 3);
        add_sphere(&spheres[15], 2400, 1300, 10, 200, 1);
        add_sphere(&spheres[16], 2800, 1000, 20, 300, 3);
        add_sphere(&spheres[17], 9000, 1000, 0, 1000, 7);
    } else if(argv[1][0] == '2'){
        N_SPHERES = 15;
        add_sphere(&spheres[0], 500, 2300, 0, 100, 0);
        add_sphere(&spheres[1], 700, 2400, 0, 100, 1);
        add_sphere(&spheres[2], 1400, 2140, 0, 100, 2);
        add_sphere(&spheres[3], 900, 2100, 0, 100, 3);
        add_sphere(&spheres[4], 9000, 1000, 0, 1000, 7);
        add_sphere(&spheres[5], 1050, 2900, 0, 100, 5);
        add_sphere(&spheres[6], 4000, 6100, 3000, 3000, 0); // grande
        add_sphere(&spheres[7], 500, 2900, 0, 100, 0);
        add_sphere(&spheres[8], 550, 2550, 0, 100, 4);
        add_sphere(&spheres[9], 2800, 1000, 20, 300, 3);
        add_sphere(&spheres[10], 1500, 2800, -30, 80, 6);
        add_sphere(&spheres[11], 1700, 3000, 10, 130, 5);
        add_sphere(&spheres[12], 6000, 3000, 10, 500, 6);
        add_sphere(&spheres[13], 5000, 3000, 10, 500, 8);
        add_sphere(&spheres[14], 2400, 1300, 10, 200, 1);
    } else {
      printf("./ray <variacao (1 ou 2)> <tamanho (1, 2, 4)>\n");
      return 0;
  }

  material *materials;
  materials = malloc(sizeof(material)*N_MATERIALS);
  add_material(&materials[0], 1, 0, 0, 0.2); //vermelho
  add_material(&materials[1], 0, 1, 0, 0.5); //verde
  add_material(&materials[2], 0, 0, 1, 0.9); //azul
  add_material(&materials[3], 1, 1, 0, 0.9); //amarelo
  add_material(&materials[4], 1, 0, 1, 0.7); //rosa
  add_material(&materials[5], 0, 1, 1, 0.5); //ciano
  add_material(&materials[6], 1, 1, 1, 1); //branco
  add_material(&materials[7], 0.18, 0.02, 0.23, 0.65); //roxo
  add_material(&materials[8], 0.1, 0.1, 0.1, 0.1); // preto


  light *lights;
  lights = malloc(sizeof(light)*N_LIGHTS);
  add_light(&lights[0], 0, 3240, -100, 1, 1, 1);
  add_light(&lights[1], 3200, 3000, -1000, 0.6, 0.7, 1);
  add_light(&lights[2], 600, 0, -100, 0.3, 0.5, 1);
  add_light(&lights[3], 1200, 0, -100, 1, 0.5, 0.2);
  add_light(&lights[4], -300, 1200, 100, 1, 1, 1);
  add_light(&lights[5], 6000, 10, -100, 1, 1, 1);


  ray r;
  //printf("%d, %d\n", (HEIGHT*my_rank)/NMAQUINAS, (HEIGHT*(my_rank+1))/NMAQUINAS);
  for (int y = (HEIGHT*my_rank)/NMAQUINAS; y < (HEIGHT*(my_rank+1))/NMAQUINAS; ++y) {
    for (int x = 0; x < WIDTH; ++x) {
      // copy_vec(black, pix_col);
      double red = 0;
      double green = 0;
      double blue = 0;
      int level = 0;
      double coef = 1.0;

      add_ray(&r,x,y,-2000,0,0,1);
      

      do{
        double t = 20000.0f;
        int currentSphere = -1;
        //procura a esfera mais proxima desse pixel
        for(int i = 0; i < N_SPHERES; i++){
          if(intersect(&spheres[i], &r, &t))
            currentSphere = i;
        }
        if (currentSphere == -1) break;

        //criando o ray atual do pixel
        vector pi;
        vector aux;
        aux = mult_vec(&(r.d), t);
        pi = soma_vec(&(r.o), &aux);

        /* Find the normal for this new vector at the point of intersection */
        vector L;
        L = sub_vec(&pi, &(spheres[currentSphere].c));
        double temp = dot(&L, &L);
        if(temp == 0) break;

        temp = 1.0f / sqrt(temp);

        vector N;
        N = mult_vec(&L, temp);
        //
        material currentMat = materials[spheres[currentSphere].material];

        //calcula o valor da luz nesse pixel
        for(int j = 0; j < N_LIGHTS; j++){
          light currentLight = lights[j];
          vector dist;
          dist = sub_vec(&currentLight.pos, &pi);
          if(dot(&N, &dist) <= 0.0f) continue;
          double t = sqrt(dot(&dist,&dist));
          if(t <= 0.0f) continue;

          ray lightRay;
          lightRay.o = pi;
          lightRay.d = mult_vec(&dist, (1/t));
          /* Lambert diffusion */
          double lambert = dot(&lightRay.d, &N) * coef;
          // double lambert = 0.0;
          red += lambert * currentLight.intensity.red * currentMat.diffuse.red;
          green += lambert * currentLight.intensity.green * currentMat.diffuse.green;
          blue += lambert * currentLight.intensity.blue * currentMat.diffuse.blue;
        }

        /* Iterate over the reflection */
        coef *= currentMat.reflection;

        /* The reflected ray start and direction */
        r.o = pi;
        double reflect = 2.0f * dot(&r.d, &N);
        vector tmp;
        tmp = mult_vec(&N, reflect);
        r.d = sub_vec(&r.d, &tmp);

        level++;

      }while((coef > 0.0f) && (level < 30));
      // printf("%d\n", (x + y*WIDTH));
      //printf("rank: %d x: %d, y: %d\n", my_rank, x, (x + (y-(HEIGHT*(my_rank))/NMAQUINAS)*WIDTH)*3 + 0);
      img[(x + (y-(HEIGHT*(my_rank))/NMAQUINAS)*WIDTH)*3 + 0] = (unsigned char)min(red*255.0f, 255.0f);
      img[(x + (y-(HEIGHT*(my_rank))/NMAQUINAS)*WIDTH)*3 + 1] = (unsigned char)min(green*255.0f, 255.0f);
      img[(x + (y-(HEIGHT*(my_rank))/NMAQUINAS)*WIDTH)*3 + 2] = (unsigned char)min(blue*255.0f, 255.0f);
    }

  }
  

// alocar matriz para cada vetor das maquinas
// rank 0 tem todas, e trabalha no vetor 0
// cada maquina vai trabalhar do y = tambloco * rank até tambloco * (rank + 1)
// formula para gravar na imagem para cada maquina vai ser (x + (y-(tambloco*rank))*WIDTH)*3 + rgb]
// no fim vamos receber de cada maquina um seu vetor e colocar ele na matriz do rank 0, usando seu rank como posicao na matriz
// vamos escrever no arquivo com for da matriz e for do vetor
  unsigned char *img_t; 
  if(my_rank == 0)
    img_t = malloc(sizeof(unsigned char)*tamanho_bloco*NMAQUINAS);

  MPI_Gather(img, tamanho_bloco, MPI_UNSIGNED_CHAR, img_t, tamanho_bloco, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
  if(!img) free(img);


  if (my_rank == 0) saveppm("image.ppm", img_t, WIDTH, HEIGHT);
  if(!img_t) free(img_t);
  MPI_Finalize();

}
