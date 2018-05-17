// #include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#ifndef NTHREADS
#define NTHREADS 4
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

// Atribui valores para a luz
void add_light(light *l, double x, double y, double z, double r, double g, double b){
  l->pos.x = x;
  l->pos.y = y;
  l->pos.z = z;
  l->intensity.red = r;
  l->intensity.green = g;
  l->intensity.blue = b;

}

// Atribui valores para o material
void add_material(material *m, double r, double g, double b, double reflec){
  m->diffuse.red = r;
  m->diffuse.green = g;
  m->diffuse.blue = b;
  m->reflection = reflec;
}

// Atribui Valores para a cor
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

// Normaliza o Vetor
void normalize(vector *v) {
  double mg = sqrt((v->x)*(v->x) + (v->y)*(v->y) + (v->z)*(v->z));
  v->x = v->x/mg;
  v->y = v->y/mg;
  v->z = v->z/mg;
}

// Copia o conteúdo de um vetor orig para o dest
void copy_vec(vector *orig, vector *dest){
  dest->x = orig->x;
  dest->y = orig->y;
  dest->z = orig->z;
}


// Cálculo do Produto escalar
double dot(vector *v, vector *w) {
  return (v->x*w->x + v->y*w->y + v->z*w->z);
}

// Atribui valores aos vetores do Ray
void add_ray(ray *r, double o_x, double o_y, double o_z, double d_x, double d_y, double d_z){
  inic_vec(&(r->o), o_x, o_y, o_z);
  inic_vec(&(r->d), d_x, d_y, d_z);
}


// Atribui valores para a Esfera
void add_sphere(sphere *s, double x, double y, double z, double r, int m){
  inic_vec(&(s->c), x, y, z);
  s->r = r;
  s->material = m;
}


// Retorna se um ray intercepta uma esfera ou não. 
// Se sim, o valor double *t, é atualizado para a distância da esfera
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


// Função que escreve todo o conteũdo renderizado no arquivo
void saveppm(char *filename, unsigned char *img, int width, int height){
	FILE *f;

	f = fopen(filename, "wb");

	fprintf(f, "P6 %d %d %d\n", width, height, 255);

	fwrite(img, 3, width*height, f);

	fclose(f);
}

int main(int argc, char const *argv[])
{
  if (argc < 3){
    printf("./ray <variacao (1 ou 2)> <tamanho (1, 2, 4)>\n");
    return 0;
  }
  int HEIGHT = 4500; 
  int WIDTH = 8000;

  // Baseado na entrada define o número de pixels
  if(argv[2][0] == '1'){
      HEIGHT = 4500; 
      WIDTH = 8000; 
  } else if(argv[2][0] == '2'){
      HEIGHT = 6188; 
      WIDTH = 11000; 
  } else if(argv[2][0] == '4'){
      HEIGHT = 9000; 
      WIDTH = 16000; 
  } else {
      printf("./ray <variacao (1 ou 2)> <tamanho (1, 2, 4)>\n");
      return 0;
  }
    int N_SPHERES = 18;

  
  int N_LIGHTS = 6;
  int N_MATERIALS = 9;
  // unsigned char img[3*WIDTH*HEIGHT];
  unsigned char *img;
  img = malloc(sizeof(unsigned char)*3*WIDTH*HEIGHT);

  sphere *spheres;
  spheres = malloc(sizeof(sphere)*N_SPHERES);

  // Com base na entrada define o setup de esferas da imagem
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

  // Predefinição de materias
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
  // Craiação das luzes
  add_light(&lights[0], 0, 3240, -100, 1, 1, 1);
  add_light(&lights[1], 3200, 3000, -1000, 0.6, 0.7, 1);
  add_light(&lights[2], 600, 0, -100, 0.3, 0.5, 1);
  add_light(&lights[3], 1200, 0, -100, 1, 0.5, 0.2);
  add_light(&lights[4], -300, 1200, 100, 1, 1, 1);
  add_light(&lights[5], 6000, 10, -100, 1, 1, 1);


  ray r;
  #pragma omp parallel for num_threads(NTHREADS) private(r) schedule(dynamic)
  for (int y = 0; y < HEIGHT; ++y) {
    for (int x = 0; x < WIDTH; ++x) {

      double red = 0;
      double green = 0;
      double blue = 0;
      int level = 0;
      double coef = 1.0;
      // Como x, y são as coordenadas da imagem, faz com que o raio passe por um pixel.
      add_ray(&r,x,y,-2000,0,0,1);

      do{
        double t = 20000.0f;
        int currentSphere = -1;
        // Procura a esfera mais proxima desse pixel
        for(int i = 0; i < N_SPHERES; i++){
          if(intersect(&spheres[i], &r, &t))
            currentSphere = i;
        }
        if (currentSphere == -1) break;

        // Ponto de Intercecção
        vector pi;
        vector aux;
        aux = mult_vec(&(r.d), t);
        pi = soma_vec(&(r.o), &aux);

        /* Find the normal for this new vector at the point of intersection */
        // Encontra a normal desse novo vetor no ponto de intercecção
        vector L;
        L = sub_vec(&pi, &(spheres[currentSphere].c));
        double temp = dot(&L, &L);
        if(temp == 0) break;

        temp = 1.0f / sqrt(temp);

        vector N;
        N = mult_vec(&L, temp);
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
          
          // Difusão de Lambert
          double lambert = dot(&lightRay.d, &N) * coef;
          red += lambert * currentLight.intensity.red * currentMat.diffuse.red;
          green += lambert * currentLight.intensity.green * currentMat.diffuse.green;
          blue += lambert * currentLight.intensity.blue * currentMat.diffuse.blue;
        }

        // Iterar pelo numero de reflexões possíveis
        coef *= currentMat.reflection;

        // Direção do inicio e direção do Ray refletido
        r.o = pi;
        double reflect = 2.0f * dot(&r.d, &N);
        vector tmp;
        tmp = mult_vec(&N, reflect);
        r.d = sub_vec(&r.d, &tmp);

        level++;

      }while((coef > 0.0f) && (level < 30));
      // Insere cor no Pixel
      img[(x + y*WIDTH)*3 + 0] = (unsigned char)min(red*255.0f, 255.0f);
      img[(x + y*WIDTH)*3 + 1] = (unsigned char)min(green*255.0f, 255.0f);
      img[(x + y*WIDTH)*3 + 2] = (unsigned char)min(blue*255.0f, 255.0f);
    }

  }

  saveppm("image.ppm", img, WIDTH, HEIGHT);
}
