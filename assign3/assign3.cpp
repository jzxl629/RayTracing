/*
CSCI 420
Assignment 3 Raytracer

Name: Zhouqian Lu
*/

#include <string.h>
#include "pic.h"

#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include <cstring>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits>

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

char *filename = 0;

// different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode = MODE_DISPLAY;

// DEBUG
bool debug = true;

// you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

// the field of view of the camera
#define fov 60.0
#define PI 3.14159265
#define smallValue 0.000001

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

typedef struct _Triangle
{
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b);

// Variable Declaration
double bottomLeft[2];
double topRight[2];

void normalize(double p[3])
{
  double len = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
  p[0] /= len;
  p[1] /= len;
  p[2] /= len;
}

void subtraction(double a[3], double b[3], double c[3])
{
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2];
}

void cross_product(double a[3], double b[3], double c[3])
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

void subtractPoints(double v1[3], double v2[3], double v3[3])
{
  v3[0] = v2[0] - v1[0];
  v3[1] = v2[1] - v1[1];
  v3[2] = v2[2] - v1[2];
}

// Calculate the distance between two points
double calDistance(double direction[3])
{
  return sqrt(pow(direction[0], 2) + pow(direction[1], 2) + pow(direction[2], 2));
}

// Determines if the shadow ray hits a triangle
bool shadowIntersectTriangle(double direction[3], double origin[3], double lightDist)
{
  for (int i = 0; i < num_triangles; i++)
  {
    double v1[3] = {0, 0, 0};
    subtractPoints(triangles[i].v[0].position, triangles[i].v[1].position, v1);
    double v2[3] = {0, 0, 0};
    subtractPoints(triangles[i].v[0].position, triangles[i].v[2].position, v2);
    double norm[3] = {0, 0, 0};
    cross_product(v1, v2, norm);
    normalize(norm);
    double D = -(triangles[i].v[0].position[0] * norm[0] + triangles[i].v[0].position[1] * norm[1] + triangles[i].v[0].position[2] * norm[2]);

    double denominator = norm[0] * direction[0] + norm[1] * direction[1] + norm[2] * direction[2];
    // ray parallel to plane, no intersection
    if (denominator == 0)
    {
      continue;
    }
    double t = -(norm[0] * origin[0] + norm[1] * origin[1] + norm[2] * origin[2] + D) / denominator;
    // the intersection is behind ray origin, no intersection
    if (t < smallValue || t > lightDist)
    {
      continue;
    }

    // Check if the point is inside triangle or not
    double intersectionPoint[3] = {origin[0] + t * direction[0], origin[1] + t * direction[1], origin[2] + t * direction[2]};

    double C[3] = {0, 0, 0};
    double edge[3] = {0, 0, 0};
    subtractPoints(triangles[i].v[0].position, triangles[i].v[1].position, edge);
    double vp[3] = {0, 0, 0};
    subtractPoints(triangles[i].v[0].position, intersectionPoint, vp);
    cross_product(edge, vp, C);
    bool testA = norm[0] * C[0] + norm[1] * C[1] + norm[2] * C[2] >= 0;

    subtractPoints(triangles[i].v[1].position, triangles[i].v[2].position, edge);
    subtractPoints(triangles[i].v[1].position, intersectionPoint, vp);
    cross_product(edge, vp, C);
    bool testB = norm[0] * C[0] + norm[1] * C[1] + norm[2] * C[2] >= 0;

    subtractPoints(triangles[i].v[2].position, triangles[i].v[0].position, edge);
    subtractPoints(triangles[i].v[2].position, intersectionPoint, vp);
    cross_product(edge, vp, C);
    bool testC = norm[0] * C[0] + norm[1] * C[1] + norm[2] * C[2] >= 0;

    if (testA && testB && testC)
    {
      return true;
    }
  }
  return false;
}

// Determines if the shadow ray hits a sphere
bool shadowIntersectSphere(double direction[3], double origin[3], double lightDist)
{
  for (int i = 0; i < num_spheres; i++)
  {
    double diff[3] = {0, 0, 0};
    subtraction(origin, spheres[i].position, diff);
    double a = 1;
    double b = 2 * (direction[0] * diff[0] + direction[1] * diff[1] + direction[2] * diff[2]);
    double c = pow(diff[0], 2) + pow(diff[1], 2) + pow(diff[2], 2) - pow(spheres[i].radius, 2);
    double k = pow(b, 2) - 4 * a * c;
    if (k < 0)
    {
      continue;
    }
    else
    {
      double t0 = (-b + sqrt(k)) / 2;
      double t1 = (-b - sqrt(k)) / 2;
      double closerT = t0 < t1 ? t0 : t1;
      if (closerT > smallValue && closerT < lightDist)
      {
        return true;
      }
    }
  }
  return false;
}

// Determines if a point is in shadow
bool isInShadow(double point[3], Light l)
{
  double direction[3] = {l.position[0] - point[0], l.position[1] - point[1], l.position[2] - point[2]};
  double lightDist = calDistance(direction);
  normalize(direction);
  // Ray = p(0) + dt
  bool hitSphere = shadowIntersectSphere(direction, point, lightDist);
  bool hitTriangle = shadowIntersectTriangle(direction, point, lightDist);
  return hitSphere || hitTriangle;
}

// Returns the closest sphere this ray hits, if any
int raySphere(double direction[3], double origin[3], double &t)
{
  int id = -1;
  for (int i = 0; i < num_spheres; i++)
  {
    double diff[3] = {0, 0, 0};
    subtraction(origin, spheres[i].position, diff);
    double a = 1;
    double b = 2 * (direction[0] * diff[0] + direction[1] * diff[1] + direction[2] * diff[2]);
    double c = pow(diff[0], 2) + pow(diff[1], 2) + pow(diff[2], 2) - pow(spheres[i].radius, 2);
    double k = pow(b, 2) - 4 * c;

    if (k < 0)
    {
      continue;
    }
    else
    {
      double t0 = (-b + sqrt(k)) / 2;
      double t1 = (-b - sqrt(k)) / 2;
      double closerT = t0 < t1 ? t0 : t1;
      if (closerT > smallValue && closerT < t)
      {
        t = closerT;
        id = i;
      }
    }
  }
  return id;
}

// Returns the closest triangle this ray hits, if any
int rayTriangle(double direction[3], double origin[3], double &t)
{
  int id = -1;
  for (int i = 0; i < num_triangles; i++)
  {
    double v1[3] = {0, 0, 0};
    subtractPoints(triangles[i].v[0].position, triangles[i].v[1].position, v1);
    double v2[3] = {0, 0, 0};
    subtractPoints(triangles[i].v[0].position, triangles[i].v[2].position, v2);
    double norm[3] = {0, 0, 0};
    cross_product(v1, v2, norm);
    normalize(norm);
    double D = -(triangles[i].v[0].position[0] * norm[0] + triangles[i].v[0].position[1] * norm[1] + triangles[i].v[0].position[2] * norm[2]);

    double denominator = norm[0] * direction[0] + norm[1] * direction[1] + norm[2] * direction[2];
    // ray parallel to plane, no intersection
    if (denominator == 0)
    {
      continue;
    }

    double tt = -(norm[0] * origin[0] + norm[1] * origin[1] + norm[2] * origin[2] + D) / denominator;
    // the intersection is behind ray origin, no intersection
    if (tt <= smallValue)
    {
      continue;
    }

    // Check if the point is inside triangle or not
    double intersectionPoint[3] = {origin[0] + tt * direction[0], origin[1] + tt * direction[1], origin[2] + tt * direction[2]};

    double C[3] = {0, 0, 0};
    double edge[3] = {0, 0, 0};
    subtractPoints(triangles[i].v[0].position, triangles[i].v[1].position, edge);
    double vp[3] = {0, 0, 0};
    subtractPoints(triangles[i].v[0].position, intersectionPoint, vp);
    cross_product(edge, vp, C);
    bool testA = norm[0] * C[0] + norm[1] * C[1] + norm[2] * C[2] >= 0;

    subtractPoints(triangles[i].v[1].position, triangles[i].v[2].position, edge);
    subtractPoints(triangles[i].v[1].position, intersectionPoint, vp);
    cross_product(edge, vp, C);
    bool testB = norm[0] * C[0] + norm[1] * C[1] + norm[2] * C[2] >= 0;

    subtractPoints(triangles[i].v[2].position, triangles[i].v[0].position, edge);
    subtractPoints(triangles[i].v[2].position, intersectionPoint, vp);
    cross_product(edge, vp, C);
    bool testC = norm[0] * C[0] + norm[1] * C[1] + norm[2] * C[2] >= 0;

    if (testA && testB && testC)
    {
      if (tt < t)
      {
        t = tt;
        id = i;
      }
    }
  }
  return id;
}

// Get color for each ray
void getColor(double direction[3], double origin[3], double color[3])
{
  // Check if it hits anything
  double tSphere = std::numeric_limits<double>::max();
  double tTriangle = std::numeric_limits<double>::max();
  int idSphere = raySphere(direction, origin, tSphere);
  int idTriangle = rayTriangle(direction, origin, tTriangle);

  if (idSphere == -1 && idTriangle == -1)
  {
    color[0] = 1.0;
    color[1] = 1.0;
    color[2] = 1.0;
    return;
  }

  // Hits triangle
  if (idSphere == -1 || tTriangle < tSphere)
  {
    double intersectionPoint[3] = {origin[0] + tTriangle * direction[0], origin[1] + tTriangle * direction[1], origin[2] + tTriangle * direction[2]};

    // Does not work when ABC is 0 when projected to AB-plane
    double ABC = (triangles[idTriangle].v[1].position[0] - triangles[idTriangle].v[0].position[0]) *
                     (triangles[idTriangle].v[2].position[1] - triangles[idTriangle].v[0].position[1]) -
                 (triangles[idTriangle].v[2].position[0] - triangles[idTriangle].v[0].position[0]) *
                     (triangles[idTriangle].v[1].position[1] - triangles[idTriangle].v[0].position[1]);

    double PBC = (triangles[idTriangle].v[1].position[0] - intersectionPoint[0]) *
                     (triangles[idTriangle].v[2].position[1] - intersectionPoint[1]) -
                 (triangles[idTriangle].v[2].position[0] - intersectionPoint[0]) *
                     (triangles[idTriangle].v[1].position[1] - intersectionPoint[1]);
    double APC = (intersectionPoint[0] - triangles[idTriangle].v[0].position[0]) *
                     (triangles[idTriangle].v[2].position[1] - triangles[idTriangle].v[0].position[1]) -
                 (triangles[idTriangle].v[2].position[0] - triangles[idTriangle].v[0].position[0]) *
                     (intersectionPoint[1] - triangles[idTriangle].v[0].position[1]);
    double alpha = 0.33;
    double beta = 0.33;
    double gamma = 0.33;
    if (ABC != 0)
    {
      alpha = PBC / ABC;
      beta = APC / ABC;
      gamma = 1 - alpha - beta;
    }
    double n[3] = {alpha * triangles[idTriangle].v[0].normal[0] + beta * triangles[idTriangle].v[1].normal[0] +
                       gamma * triangles[idTriangle].v[2].normal[0],
                   alpha * triangles[idTriangle].v[0].normal[1] + beta * triangles[idTriangle].v[1].normal[1] +
                       gamma * triangles[idTriangle].v[2].normal[1],
                   alpha * triangles[idTriangle].v[0].normal[2] + beta * triangles[idTriangle].v[1].normal[2] +
                       gamma * triangles[idTriangle].v[2].normal[2]};
    normalize(n);

    for (int i = 0; i < num_lights; i++)
    {
      double l[3] = {lights[i].position[0] - intersectionPoint[0], lights[i].position[1] - intersectionPoint[1], lights[i].position[2] - intersectionPoint[2]};
      normalize(l);

      if (!isInShadow(intersectionPoint, lights[i]))
      {
        // Diffuse
        double ln = n[0] * l[0] + n[1] * l[1] + n[2] * l[2];
        // Clamps ln to zero if it is negative
        if (ln < 0)
        {
          ln = 0;
        }

        // specular
        double v[3] = {-intersectionPoint[0], -intersectionPoint[1], -intersectionPoint[2]};
        normalize(v);
        double r[3] = {0, 0, 0};
        if (ln == 0)
        {
          r[0] = -l[0];
          r[1] = -l[1];
          r[2] = -l[2];
        }
        else
        {
          r[0] = 2 * ln * n[0] - l[0];
          r[1] = 2 * ln * n[1] - l[1];
          r[2] = 2 * ln * n[2] - l[2];
        }
        double rv = r[0] * v[0] + r[1] * v[1] + r[2] * v[2];
        if (rv < 0)
        {
          rv = 0;
        }

        color[0] += lights[i].color[0] * (alpha * triangles[idTriangle].v[0].color_diffuse[0] + beta * triangles[idTriangle].v[1].color_diffuse[0] + gamma * triangles[idTriangle].v[2].color_diffuse[0]) * ln;
        color[1] += lights[i].color[1] * (alpha * triangles[idTriangle].v[0].color_diffuse[1] + beta * triangles[idTriangle].v[1].color_diffuse[1] + gamma * triangles[idTriangle].v[2].color_diffuse[1]) * ln;
        color[2] += lights[i].color[2] * (alpha * triangles[idTriangle].v[0].color_diffuse[2] + beta * triangles[idTriangle].v[1].color_diffuse[2] + gamma * triangles[idTriangle].v[2].color_diffuse[2]) * ln;

        color[0] += lights[i].color[0] * (alpha * triangles[idTriangle].v[0].color_specular[0] + beta * triangles[idTriangle].v[1].color_specular[0] + gamma * triangles[idTriangle].v[2].color_specular[0]) * pow(rv, triangles[idTriangle].v[0].shininess);
        color[1] += lights[i].color[1] * (alpha * triangles[idTriangle].v[0].color_specular[1] + beta * triangles[idTriangle].v[1].color_specular[1] + gamma * triangles[idTriangle].v[2].color_specular[1]) * pow(rv, triangles[idTriangle].v[1].shininess);
        color[2] += lights[i].color[2] * (alpha * triangles[idTriangle].v[0].color_specular[2] + beta * triangles[idTriangle].v[1].color_specular[2] + gamma * triangles[idTriangle].v[2].color_specular[2]) * pow(rv, triangles[idTriangle].v[2].shininess);
      }
    }
  }
  else
  {
    // Hits sphere
    double intersectionPoint[3] = {origin[0] + tSphere * direction[0], origin[1] + tSphere * direction[1], origin[2] + tSphere * direction[2]};
    for (int i = 0; i < num_lights; i++)
    {
      double l[3] = {lights[i].position[0] - intersectionPoint[0], lights[i].position[1] - intersectionPoint[1], lights[i].position[2] - intersectionPoint[2]};
      normalize(l);
      double n[3] = {intersectionPoint[0] - spheres[idSphere].position[0], intersectionPoint[1] - spheres[idSphere].position[1], intersectionPoint[2] - spheres[idSphere].position[2]};
      normalize(n);
      // If shadow ray is not blocked by an object
      if (!isInShadow(intersectionPoint, lights[i]))
      {
        // Diffuse

        double ln = n[0] * l[0] + n[1] * l[1] + n[2] * l[2];
        // Clamps ln to zero if it is negative
        if (ln < 0)
        {
          ln = 0;
        }

        // Specular
        double v[3] = {-intersectionPoint[0], -intersectionPoint[1], -intersectionPoint[2]};
        normalize(v);
        // reflection r = 2(l dot n)n - l
        double r[3] = {0, 0, 0};
        if (ln == 0)
        {
          r[0] = -l[0];
          r[1] = -l[1];
          r[2] = -l[2];
        }
        else
        {
          r[0] = 2 * ln * n[0] - l[0];
          r[1] = 2 * ln * n[1] - l[1];
          r[2] = 2 * ln * n[2] - l[2];
        }
        double rv = r[0] * v[0] + r[1] * v[1] + r[2] * v[2];
        if (rv < 0)
        {
          rv = 0;
        }
        color[0] += lights[i].color[0] * spheres[idSphere].color_diffuse[0] * ln;
        color[1] += lights[i].color[1] * spheres[idSphere].color_diffuse[1] * ln;
        color[2] += lights[i].color[2] * spheres[idSphere].color_diffuse[2] * ln;

        color[0] += lights[i].color[0] * spheres[idSphere].color_specular[0] * pow(rv, spheres[idSphere].shininess);
        color[1] += lights[i].color[1] * spheres[idSphere].color_specular[1] * pow(rv, spheres[idSphere].shininess);
        color[2] += lights[i].color[2] * spheres[idSphere].color_specular[2] * pow(rv, spheres[idSphere].shininess);
      }
    }
  }
}

// Initializes image plane
void initImagePlane()
{
  double aspectRatio = (double)WIDTH / HEIGHT;
  double angle = fov * PI / 180;

  // Calculates vertices of the bottom left corner of the image plane
  bottomLeft[0] = -tan(angle / 2) * aspectRatio;
  bottomLeft[1] = -tan(angle / 2);

  topRight[0] = tan(angle / 2) * aspectRatio;
  topRight[1] = tan(angle / 2);
}

// MODIFY THIS FUNCTION
void draw_scene()
{
  initImagePlane();
  double iw = 2 * topRight[0];
  double ih = 2 * topRight[1];

  // printf("Max double: %f", std::numeric_limits<double>::max());
  //  simple output
  for (unsigned int x = 0; x < WIDTH; x++)
  {
    glPointSize(2.0);
    glBegin(GL_POINTS);
    for (unsigned int y = 0; y < HEIGHT; y++)
    {
      double origin[3] = {0, 0, 0};
      // Fire a ray to the middle of the pixel
      double dir[3] = {bottomLeft[0] + (1.0 * x + 0.5) * iw / WIDTH, bottomLeft[1] + (1.0 * y + 0.5) * ih / HEIGHT, -1};
      normalize(dir);
      double color[3] = {ambient_light[0], ambient_light[1], ambient_light[2]};
      // Get the color for the ray
      getColor(dir, origin, color);
      if (color[0] > 1.0)
      {
        color[0] = 1.0;
      }
      if (color[1] > 1.0)
      {
        color[1] = 1.0;
      }
      if (color[2] > 1.0)
      {
        color[2] = 1.0;
      }
      // Draw pixel
      plot_pixel(x, y, (int)(color[0] * 255), (int)(color[1] * 255), (int)(color[2] * 255));
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n");
  fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((double)r) / 256.f, ((double)g) / 256.f, ((double)b) / 256.f);
  glVertex2i(x, y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[HEIGHT - y - 1][x][0] = r;
  buffer[HEIGHT - y - 1][x][1] = g;
  buffer[HEIGHT - y - 1][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x, y, r, g, b);
  if (mode == MODE_JPEG)
    plot_pixel_jpeg(x, y, r, g, b);
}

void save_jpg()
{
  Pic *in = NULL;

  in = pic_alloc(640, 480, 3, NULL);
  printf("Saving JPEG file: %s\n", filename);

  memcpy(in->pix, buffer, 3 * WIDTH * HEIGHT);
  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);
}

void parse_check(char *expected, char *found)
{
  if (strcasecmp(expected, found))
  {
    char error[100];
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE *file, char *check, double p[3])
{
  char str[100];
  fscanf(file, "%s", str);
  parse_check(check, str);
  fscanf(file, "%lf %lf %lf", &p[0], &p[1], &p[2]);
  printf("%s %lf %lf %lf\n", check, p[0], p[1], p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file, "%s", str);
  parse_check("rad:", str);
  fscanf(file, "%lf", r);
  printf("rad: %f\n", *r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file, "%s", s);
  parse_check("shi:", s);
  fscanf(file, "%lf", shi);
  printf("shi: %f\n", *shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv, "r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file, "%i", &number_of_objects);

  printf("number of objects: %i\n", number_of_objects);
  char str[200];

  parse_doubles(file, "amb:", ambient_light);

  for (i = 0; i < number_of_objects; i++)
  {
    fscanf(file, "%s\n", type);
    printf("%s\n", type);
    if (strcasecmp(type, "triangle") == 0)
    {

      printf("found triangle\n");
      int j;

      for (j = 0; j < 3; j++)
      {
        parse_doubles(file, "pos:", t.v[j].position);
        parse_doubles(file, "nor:", t.v[j].normal);
        parse_doubles(file, "dif:", t.v[j].color_diffuse);
        parse_doubles(file, "spe:", t.v[j].color_specular);
        parse_shi(file, &t.v[j].shininess);
      }

      if (num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if (strcasecmp(type, "sphere") == 0)
    {
      printf("found sphere\n");

      parse_doubles(file, "pos:", s.position);
      parse_rad(file, &s.radius);
      parse_doubles(file, "dif:", s.color_diffuse);
      parse_doubles(file, "spe:", s.color_specular);
      parse_shi(file, &s.shininess);

      if (num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if (strcasecmp(type, "light") == 0)
    {
      printf("found light\n");
      parse_doubles(file, "pos:", l.position);
      parse_doubles(file, "col:", l.color);

      if (num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n", type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0, WIDTH, 0, HEIGHT, 1, -1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0, 0, 0, 0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  // hack to make it only once
  static int once = 0;
  if (!once)
  {
    draw_scene();
    if (mode == MODE_JPEG)
    {
      save_jpg();
    }
  }
  once = 1;
}

int main(int argc, char **argv)
{
  if (argc < 2 || argc > 3)
  {
    printf("usage: %s <scenefile> [jpegname]\n", argv[0]);
    exit(0);
  }
  if (argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
    printf("generate a file");
  }
  else if (argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc, argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0, 0);
  glutInitWindowSize(WIDTH, HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
