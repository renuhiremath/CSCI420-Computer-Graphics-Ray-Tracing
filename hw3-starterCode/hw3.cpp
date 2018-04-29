/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Renu Hiremath
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <imageIO.h>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif


#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;
const double epsilon = 1e-15;
bool antialiasing = false;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

struct Vector
{
  double x;
  double y;
  double z;
};

struct Point
{
  double x;
  double y;
  double z;
};

struct Ray
{
  Point origin;
  Vector direction;
};

struct Color
{
  double r;
  double g;
  double b;
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;
int x, y;
void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

void clampValueToZero(double &value)
{
  if (value<0)
    value = 0.0;
}

void clampValueToOne(double &value)
{
  if (value>1)
    value = 1.0;
}

void normalize(Vector &v)
{
  double mag = (v.x*v.x) + (v.y*v.y) + (v.z*v.z);
  if (mag!=0)
  {
    mag = sqrt(mag);
    v.x = v.x/mag;
    v.y = v.y/mag;
    v.z = v.z/mag;
  }
  else
  {
    v.x=0;
    v.y=0;
    v.z=0;
  }
}

void computeCrossProduct(Vector a, Vector b, Vector &c)
{
  c.x = a.y*b.z - b.y*a.z;
  c.y = a.z*b.x - b.z*a.x;
  c.z = a.x*b.y - b.x*a.y;
}

double computeMagnitudeOfVector(Vector a)
{
  return (pow(a.x,2) + pow(a.y,2) + pow(a.z,2));
}

void computeVector(Point a, Point b, Vector& c)
{
  c.x = a.x - b.x;
  c.y = a.y - b.y;
  c.z = a.z - b.z;
}

double computeAreaOfTriangle(Point a, Point b, Point c)
{
  Vector b_a,c_a;
  Vector area_vec;
  computeVector(b,a,b_a);//B-A
  computeVector(c,a,c_a);//C-A

  computeCrossProduct(b_a,c_a,area_vec);

  double mag = pow(area_vec.x,2) + pow(area_vec.y,2) + pow(area_vec.z,2);
  if (mag>=0)
    return (0.5 * sqrt(mag));
  else
    return 0.0;
}

void computeDotProduct(Vector a, Vector b, double &value)
{
  value = (a.x*b.x) + (a.y*b.y) + (a.z*b.z);
}

void computeRayFromCamera(int x, int y, Ray &ray, float difX, float difY)
{
  double ratio = (double)(WIDTH)/(double)HEIGHT;
  double angle = tan((fov/2)* ((22.0/7) / 180.0));
  double width, height;
  width = 2 * angle * ratio;
  height = 2 * angle;

  Point bottomLeftCorner;
  bottomLeftCorner.x = -1.0 * (width) / 2.0;
  bottomLeftCorner.y = -1.0 * (height) / 2.0;
  bottomLeftCorner.z = -1;

  Vector camera;
  camera.x = 0.0;
  camera.y = 0.0;
  camera.z = 0.0;

  ray.origin.x = 0.0;
  ray.origin.y = 0.0;
  ray.origin.z = 0.0;
  ray.direction.x = (bottomLeftCorner.x + ((x+0.5-difX) * width/(double)WIDTH)) - camera.x;
  ray.direction.y = (bottomLeftCorner.y + ((y+0.5-difY) * height/(double)HEIGHT)) - camera.y;
  ray.direction.z = bottomLeftCorner.z - camera.z;
  normalize(ray.direction);
}

bool pointInTriangleTest(Point c, Triangle triangle)
{
    Point c1, c2, c0;
    c0.x = triangle.v[0].position[0];
    c0.y = triangle.v[0].position[1];
    c0.z = triangle.v[0].position[2];

    c1.x = triangle.v[1].position[0];
    c1.y = triangle.v[1].position[1];
    c1.z = triangle.v[1].position[2];

    c2.x = triangle.v[2].position[0];
    c2.y = triangle.v[2].position[1];
    c2.z = triangle.v[2].position[2];

    //calculate the barycentric co efficients
    double areaTri = computeAreaOfTriangle(c0,c1,c2);
    double alpha = computeAreaOfTriangle(c,c1,c2)/areaTri;
    double beta = computeAreaOfTriangle(c0,c,c2)/areaTri;
    double gamma = computeAreaOfTriangle(c,c0,c1)/areaTri;

    double error = 1e-15 / 2.0;
    //consider a reasonable error due to the float division

    if ((alpha + beta + gamma)>=(1.0-error) && (alpha + beta + gamma)<=(1.0+error) )
    {
      //point lies inside the triangle
      return true;
    }
    return false;
}

bool getIntersectionWithTriangle(Triangle triangle, Ray ray, Point &intersectionPoint)
{
  //get the normal for the plane
  Vector normal,b_a,c_a;
  //B-A
  b_a.x = triangle.v[1].position[0] - triangle.v[0].position[0];
  b_a.y = triangle.v[1].position[1] - triangle.v[0].position[1];
  b_a.z = triangle.v[1].position[2] - triangle.v[0].position[2];

  //C-A
  c_a.x = triangle.v[2].position[0] - triangle.v[0].position[0];
  c_a.y = triangle.v[2].position[1] - triangle.v[0].position[1];
  c_a.z = triangle.v[2].position[2] - triangle.v[0].position[2];

  computeCrossProduct(b_a,c_a,normal);
  normalize(normal);

  double nDOTp = 0, nDOTd = 0;
  computeDotProduct(normal, ray.direction, nDOTd);
  if (abs(nDOTd) < epsilon)
  //no intersection - the plane and the ray are parallel
  {
    return false;
  }

  //compute d
  double d = 0.0;
  //substituing values in the eq of the plane : ax+by+cz+d=0
  d = - ((normal.x * triangle.v[0].position[0])
        + (normal.y * triangle.v[0].position[1])
        + (normal.z * triangle.v[0].position[2]));

  nDOTp = (normal.x * ray.origin.x) + (normal.y * ray.origin.y) + (normal.z * ray.origin.z);
  double t = (-(nDOTp + d))/nDOTd;
  if (t<0)
  {
    //intersection behind the origin
    return false;
  }
  else
  {
    intersectionPoint.x = ray.origin.x + ray.direction.x*t;
    intersectionPoint.y = ray.origin.y + ray.direction.y*t;
    intersectionPoint.z = ray.origin.z + ray.direction.z*t;

    if (pointInTriangleTest(intersectionPoint,triangle))
    {
      return true;
    }
    else
    {
      return false;
    }
  }
}

bool getIntersectionWithSphere(Sphere sphere, Ray ray, Point &intersectionPoint)
{
  double b,c;
  Vector d;

  d.x = ray.origin.x - sphere.position[0];
  d.y = ray.origin.y - sphere.position[1];
  d.z = ray.origin.z - sphere.position[2];

  b = 2 * (ray.direction.x * d.x + ray.direction.y * d.y + ray.direction.z * d.z);

  c = pow(d.x,2) + pow(d.y,2) + pow(d.z,2) - pow(sphere.radius,2);

  double val = pow(b,2)-(4*c);
  if (val<0)
  //no roots - no intersection
  {
    return false;
  }
  val= sqrt(val);

  double t0,t1;
  t0 = (-b - val)/2;
  t1 = (-b + val)/2;

  if (t0>=0 and t1>=0)
  {
    t0 = fmin(t0,t1);
    intersectionPoint.x = ray.origin.x + ray.direction.x*t0;
    intersectionPoint.y = ray.origin.y + ray.direction.y*t0;
    intersectionPoint.z = ray.origin.z + ray.direction.z*t0;
    return true;
  }
  return false;
}

void computeColorForSphere(int index, Ray ray, Point intersectionPoint, Color &color)
{
  Ray shadowRay;
  bool objectLit = true;
  Light light;
  for (int l=0; l<num_lights; l++)
  {
    light.position[0] = lights[l].position[0];
    light.position[1] = lights[l].position[1];
    light.position[2] = lights[l].position[2];
    light.color[0] = lights[l].color[0];
    light.color[1] = lights[l].color[1];
    light.color[2] = lights[l].color[2];

    //get shadow ray
    shadowRay.origin.x = intersectionPoint.x;
    shadowRay.origin.y = intersectionPoint.y;
    shadowRay.origin.z = intersectionPoint.z;

    shadowRay.direction.x = light.position[0] - intersectionPoint.x;
    shadowRay.direction.y = light.position[1] - intersectionPoint.y;
    shadowRay.direction.z = light.position[2] - intersectionPoint.z;
    normalize(shadowRay.direction);

    Point shadowInteractionPoint;
    objectLit = true;
    //compute collision between other objects
    for (int i=0; i<num_spheres; i++)
    {
      if (i == index)
      {
        continue;
      }
      if (getIntersectionWithSphere(spheres[i],shadowRay, shadowInteractionPoint))
      {
        //make sure the intersection point is on the incident ray
        Vector x,y;
        Point l;
        l.x = light.position[0];
        l.y = light.position[1];
        l.z = light.position[2];
        computeVector(shadowInteractionPoint, intersectionPoint, x);
        computeVector(l, intersectionPoint, y);
        if (computeMagnitudeOfVector(x) < computeMagnitudeOfVector(y))
        {
          objectLit = false;
          break;
        }
      }
    }

    if (objectLit)
    {
      for (int j=0; j<num_triangles; j++)
      {
        if (getIntersectionWithTriangle(triangles[j],shadowRay, shadowInteractionPoint))
        {
          //make sure the intersection point is on the incident ray
          Vector x,y;
          Point l;
          l.x = light.position[0];
          l.y = light.position[1];
          l.z = light.position[2];
          computeVector(shadowInteractionPoint, intersectionPoint, x);
          computeVector(l, intersectionPoint, y);
          if (computeMagnitudeOfVector(x) < computeMagnitudeOfVector(y))
          {
            objectLit = false;
            break;
          }
        }
      }
    }
    //if the object is still lit, compute color
    if (objectLit)
    {
      Vector normal,l,r,v;
      double lDOTn;
      double rDOTv;

      //compute normal
      normal.x = (intersectionPoint.x - spheres[index].position[0]);
      normal.y = (intersectionPoint.y - spheres[index].position[1]);
      normal.z = (intersectionPoint.z - spheres[index].position[2]);
      normalize(normal);

      //get the light vector
      l.x = light.position[0] - intersectionPoint.x;
      l.y = light.position[1] - intersectionPoint.y;
      l.z = light.position[2] - intersectionPoint.z;
      normalize(l);

      //calculate L.N
      computeDotProduct(l, normal, lDOTn);
      clampValueToZero(lDOTn);
      clampValueToOne(lDOTn);

      //complute the reflected ray
      r.x = 2*(lDOTn)*normal.x - l.x;
      r.y = 2*(lDOTn)*normal.y - l.y;
      r.z = 2*(lDOTn)*normal.z - l.z;
      normalize(r);

      //compute v
      Point camera;
      camera.x = 0.0;
      camera.y = 0.0;
      camera.z = 0.0;
      computeVector(camera, intersectionPoint, v);
      normalize(v);

      //calculate R.V
      computeDotProduct(r, v, rDOTv);
      clampValueToZero(rDOTv);
      clampValueToOne(rDOTv);

      color.r += light.color[0] * ((spheres[index].color_diffuse[0] * lDOTn) + (spheres[index].color_specular[0] * pow(rDOTv, spheres[index].shininess)));
      color.g += light.color[1] * ((spheres[index].color_diffuse[1] * lDOTn) + (spheres[index].color_specular[1] * pow(rDOTv, spheres[index].shininess)));
      color.b += light.color[2] * ((spheres[index].color_diffuse[2] * lDOTn) + (spheres[index].color_specular[2] * pow(rDOTv, spheres[index].shininess)));
    }
  }
}

void getColorFromSpheresInteraction(Ray ray, Color &finalColor, double &closestInteraction)
{
  Point intersectionPoint;
  Color color;
  bool interaction = false;

  for (int i =0; i<num_spheres; i++)
  {
    if (getIntersectionWithSphere(spheres[i], ray, intersectionPoint))
    {
      if (intersectionPoint.z > closestInteraction)
      {
        interaction = true;
        color.r = 0;
        color.g = 0;
        color.b = 0;
        computeColorForSphere(i, ray, intersectionPoint, color);
        closestInteraction = intersectionPoint.z;
      }
    }
  }

  if (interaction)
  //if there is interaction, write the new color
  {
    finalColor.r = color.r;
    finalColor.g = color.g;
    finalColor.b = color.b;
  }
}

void computeColorForTriangle(int index, Ray ray, Point intersectionPoint, Color &color)
{
  Ray shadowRay;
  bool objectLit = true;
  Light light;
  for (int l=0; l<num_lights; l++)
  {
    light.position[0] = lights[l].position[0];
    light.position[1] = lights[l].position[1];
    light.position[2] = lights[l].position[2];
    light.color[0] = lights[l].color[0];
    light.color[1] = lights[l].color[1];
    light.color[2] = lights[l].color[2];

    //get shadow ray
    shadowRay.origin.x = intersectionPoint.x;
    shadowRay.origin.y = intersectionPoint.y;
    shadowRay.origin.z = intersectionPoint.z;

    shadowRay.direction.x = light.position[0] - intersectionPoint.x;
    shadowRay.direction.y = light.position[1] - intersectionPoint.y;
    shadowRay.direction.z = light.position[2] - intersectionPoint.z;
    normalize(shadowRay.direction);

    Point shadowInteractionPoint;
    objectLit = true;
    //compute collision between other objects
    for (int i=0; i<num_spheres; i++)
    {
      if (getIntersectionWithSphere(spheres[i], shadowRay, shadowInteractionPoint))
      {
        //make sure the intersection point is on the incident ray
        Vector x,y;
        Point l;
        l.x = light.position[0];
        l.y = light.position[1];
        l.z = light.position[2];
        computeVector(shadowInteractionPoint, intersectionPoint, x);
        computeVector(l, intersectionPoint, y);
        if (computeMagnitudeOfVector(x) < computeMagnitudeOfVector(y))
        {
          objectLit = false;
          break;
        }
      }
    }

    if (objectLit)
    {
      for (int j=0; j<num_triangles; j++)
      {
        if (j == index)
        {
          continue;
        }
        if (getIntersectionWithTriangle(triangles[j], shadowRay, shadowInteractionPoint))
        {
          //make sure the intersection point is on the incident ray
          Vector x,y;
          Point l;
          l.x = light.position[0];
          l.y = light.position[1];
          l.z = light.position[2];
          computeVector(shadowInteractionPoint, intersectionPoint, x);
          computeVector(l, intersectionPoint, y);
          if (computeMagnitudeOfVector(x) < computeMagnitudeOfVector(y))
          {
            objectLit = false;
            break;
          }
        }
      }
    }

    //if the object is still lit, compute color
    if (objectLit)
    {
      Point c1, c2, c0;
      c0.x = triangles[index].v[0].position[0];
      c0.y = triangles[index].v[0].position[1];
      c0.z = triangles[index].v[0].position[2];

      c1.x = triangles[index].v[1].position[0];
      c1.y = triangles[index].v[1].position[1];
      c1.z = triangles[index].v[1].position[2];

      c2.x = triangles[index].v[2].position[0];
      c2.y = triangles[index].v[2].position[1];
      c2.z = triangles[index].v[2].position[2];

      //calculate the barycentric co efficients
      double areaTri = computeAreaOfTriangle(c0,c1,c2);
      double alpha = computeAreaOfTriangle(intersectionPoint,c1,c2)/areaTri;
      double beta = computeAreaOfTriangle(c0,intersectionPoint,c2)/areaTri;
      double gamma = 1.0 - alpha - beta;

      Vector normal,l,r,v;
      double diffuseCoeff[3];
      double specularCoeff[3];
      double shininess;
      double lDOTn;
      double rDOTv;

      //interpolate n
      normal.x = (alpha * triangles[index].v[0].normal[0]) + (beta * triangles[index].v[1].normal[0]) + (gamma * triangles[index].v[2].normal[0]);
      normal.y = (alpha * triangles[index].v[0].normal[1]) + (beta * triangles[index].v[1].normal[1]) + (gamma * triangles[index].v[2].normal[1]);
      normal.z = (alpha * triangles[index].v[0].normal[2]) + (beta * triangles[index].v[1].normal[2]) + (gamma * triangles[index].v[2].normal[2]);
      normalize(normal);

      //interpolate diffuse coeff
      diffuseCoeff[0] = (alpha * triangles[index].v[0].color_diffuse[0]) + (beta * triangles[index].v[1].color_diffuse[0]) + (gamma * triangles[index].v[2].color_diffuse[0]);
      diffuseCoeff[1] = (alpha * triangles[index].v[0].color_diffuse[1]) + (beta * triangles[index].v[1].color_diffuse[1]) + (gamma * triangles[index].v[2].color_diffuse[1]);
      diffuseCoeff[2] = (alpha * triangles[index].v[0].color_diffuse[2]) + (beta * triangles[index].v[1].color_diffuse[2]) + (gamma * triangles[index].v[2].color_diffuse[2]);

      //interpolate specular coeff
      specularCoeff[0] = (alpha * triangles[index].v[0].color_specular[0]) + (beta * triangles[index].v[1].color_specular[0]) + (gamma * triangles[index].v[2].color_specular[0]);
      specularCoeff[1] = (alpha * triangles[index].v[0].color_specular[1]) + (beta * triangles[index].v[1].color_specular[1]) + (gamma * triangles[index].v[2].color_specular[1]);
      specularCoeff[2] = (alpha * triangles[index].v[0].color_specular[2]) + (beta * triangles[index].v[1].color_specular[2]) + (gamma * triangles[index].v[2].color_specular[2]);

      //interpolate shininess
      shininess = (alpha * triangles[index].v[0].shininess) + (beta * triangles[index].v[1].shininess) + (gamma * triangles[index].v[2].shininess);

      //get the light vector
      l.x = light.position[0] - intersectionPoint.x;
      l.y = light.position[1] - intersectionPoint.y;
      l.z = light.position[2] - intersectionPoint.z;
      normalize(l);

      //calculate L.N
      computeDotProduct(l, normal, lDOTn);
      clampValueToZero(lDOTn);
      clampValueToOne(lDOTn);

      //complute the reflected ray
      r.x = 2*(lDOTn)*normal.x - l.x;
      r.y = 2*(lDOTn)*normal.y - l.y;
      r.z = 2*(lDOTn)*normal.z - l.z;
      normalize(r);

      //compute v
      Point camera;
      camera.x = 0.0;
      camera.y = 0.0;
      camera.z = 0.0;
      computeVector(camera, intersectionPoint, v);
      normalize(v);

      //calculate R.V
      computeDotProduct(r, v, rDOTv);
      clampValueToZero(rDOTv);
      clampValueToOne(rDOTv);

      color.r += light.color[0] * ((diffuseCoeff[0] * lDOTn) + (specularCoeff[0] * pow(rDOTv, shininess)));
      color.g += light.color[1] * ((diffuseCoeff[1] * lDOTn) + (specularCoeff[1] * pow(rDOTv, shininess)));
      color.b += light.color[2] * ((diffuseCoeff[2] * lDOTn) + (specularCoeff[2] * pow(rDOTv, shininess)));

    }
  }
}

void getColorFromTrianglesInteraction(Ray ray,Color &finalColor,double &closestInteraction)
{
  Point intersectionPoint;
  Color color;
  bool interaction = false;

  for (int i =0; i<num_triangles; i++)
  {
    if (getIntersectionWithTriangle(triangles[i], ray, intersectionPoint))
    {
      if (intersectionPoint.z > closestInteraction)
      {
        interaction = true;
        color.r = 0;
        color.g = 0;
        color.b = 0;

        computeColorForTriangle(i, ray, intersectionPoint, color);
        closestInteraction = intersectionPoint.z;
      }
    }
  }

  if (interaction)
  {
    finalColor.r = color.r;
    finalColor.g = color.g;
    finalColor.b = color.b;
  }
}

void draw_scene()
{
  Color finalColor;
  Ray ray;
  double closestInteraction = -pow(10,10);

  for(x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);
    glBegin(GL_POINTS);
    for(y=0; y<HEIGHT; y++)
    {
      finalColor.r = 1.0;
      finalColor.g = 1.0;
      finalColor.b = 1.0;
      closestInteraction = -pow(10,10);

      if (antialiasing)
      {
        Color color1, color2, color3, color4;

        color1.r = 1.0;
        color1.g = 1.0;
        color1.b = 1.0;

        color2.r = 1.0;
        color2.g = 1.0;
        color2.b = 1.0;

        color3.r = 1.0;
        color3.g = 1.0;
        color3.b = 1.0;

        color4.r = 1.0;
        color4.g = 1.0;
        color4.b = 1.0;

        closestInteraction = -pow(10,10);
        //calculate the ray from camera to this point
        computeRayFromCamera(x, y, ray, -0.25, -0.25);

        //check collisions with spheres
        getColorFromSpheresInteraction(ray, color1, closestInteraction);

        //check collisions with triangles
        getColorFromTrianglesInteraction(ray,color1,closestInteraction);

        //add ambient color
        color1.r += ambient_light[0];
        color1.g += ambient_light[1];
        color1.b += ambient_light[2];

        clampValueToOne(color1.r);
        clampValueToOne(color1.g);
        clampValueToOne(color1.b);

        clampValueToZero(color1.r);
        clampValueToZero(color1.g);
        clampValueToZero(color1.b);

        closestInteraction = -pow(10,10);
        //calculate the ray from camera to this point
        computeRayFromCamera(x, y, ray, -0.25, 0.25);

        //check collisions with spheres
        getColorFromSpheresInteraction(ray, color2, closestInteraction);

        //check collisions with triangles
        getColorFromTrianglesInteraction(ray,color2,closestInteraction);

        //add ambient color
        color2.r += ambient_light[0];
        color2.g += ambient_light[1];
        color2.b += ambient_light[2];

        clampValueToOne(color2.r);
        clampValueToOne(color2.g);
        clampValueToOne(color2.b);

        clampValueToZero(color2.r);
        clampValueToZero(color2.g);
        clampValueToZero(color2.b);

        closestInteraction = -pow(10,10);
        //calculate the ray from camera to this point
        computeRayFromCamera(x, y, ray, 0.25, -0.25);

        //check collisions with spheres
        getColorFromSpheresInteraction(ray, color3, closestInteraction);

        //check collisions with triangles
        getColorFromTrianglesInteraction(ray,color3,closestInteraction);

        //add ambient color
        color3.r += ambient_light[0];
        color3.g += ambient_light[1];
        color3.b += ambient_light[2];

        clampValueToOne(color3.r);
        clampValueToOne(color3.g);
        clampValueToOne(color3.b);

        clampValueToZero(color3.r);
        clampValueToZero(color3.g);
        clampValueToZero(color3.b);

        closestInteraction = -pow(10,10);
        //calculate the ray from camera to this point
        computeRayFromCamera(x, y, ray, 0.25, 0.25);

        //check collisions with spheres
        getColorFromSpheresInteraction(ray, color4, closestInteraction);

        //check collisions with triangles
        getColorFromTrianglesInteraction(ray,color4,closestInteraction);

        //add ambient color
        color4.r += ambient_light[0];
        color4.g += ambient_light[1];
        color4.b += ambient_light[2];

        clampValueToOne(color4.r);
        clampValueToOne(color4.g);
        clampValueToOne(color4.b);

        clampValueToZero(color4.r);
        clampValueToZero(color4.g);
        clampValueToZero(color4.b);

        finalColor.r = (color1.r + color2.r + color3.r + color4.r)/4;
        finalColor.g = (color1.g + color2.g + color3.g + color4.g)/4;
        finalColor.b = (color1.b + color2.b + color3.b + color4.b)/4;

        clampValueToOne(finalColor.r);
        clampValueToOne(finalColor.g);
        clampValueToOne(finalColor.b);

        clampValueToZero(finalColor.r);
        clampValueToZero(finalColor.g);
        clampValueToZero(finalColor.b);
      }
      else
      {
        //calculate the ray from camera to this point
        computeRayFromCamera(x, y, ray, 0, 0);

        //check collisions with spheres
        getColorFromSpheresInteraction(ray, finalColor, closestInteraction);

        //check collisions with triangles
        getColorFromTrianglesInteraction(ray,finalColor,closestInteraction);

        //add ambient color
        finalColor.r += ambient_light[0];
        finalColor.g += ambient_light[1];
        finalColor.b += ambient_light[2];

        clampValueToOne(finalColor.r);
        clampValueToOne(finalColor.g);
        clampValueToOne(finalColor.b);

        clampValueToZero(finalColor.r);
        clampValueToZero(finalColor.g);
        clampValueToZero(finalColor.b);
      }

      plot_pixel(x, y, finalColor.r*255, finalColor.g*255, finalColor.b*255);
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
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
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 4))
  {
    printf ("Usage: %s <input scenefile> useAntialiasing [output jpegname] \n", argv[0]);
    exit(0);
  }
  if(argc == 4)
  {
    mode = MODE_JPEG;
    filename = argv[3];
    if (strcmp(argv[2],"true")==0 || strcmp(argv[2],"True")==0 || strcmp(argv[2],"TRUE")==0 )
    {
      antialiasing = true;
    }
  }
  if(argc == 3)
  {
    mode = MODE_DISPLAY;
    if (strcmp(argv[2],"true")==0 || strcmp(argv[2],"True")==0 || strcmp(argv[2],"TRUE")==0 )
    {
      antialiasing = true;
    }
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
