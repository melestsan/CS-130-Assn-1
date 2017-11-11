/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>

using namespace std;

/**
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h

/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}

// Primitive structs
struct vertex {
	vec3 color;
	vec4 position;
};

struct triangle {
	vertex a;
	vertex b;
	vertex c;
};

// Global variables
bool drawmode = true; // true = triangle, false = quad
bool matmode = true; // true = proj, false = modelview
mat4 currentModelMatrix;
mat4 currentProjMatrix;
vec3 currentColor;
vector<vertex> listOfVertices;
vector<triangle> listOfTriangles;
vector<mat4> matrixStackModel;
vector<mat4> matrixStackProj;

/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
 
/**
 * A function to help determine the bounding box of the passed
 * in triangle
 */
void determineBoundingBox(const triangle& tri,
			  unsigned& start, 
			  unsigned& end, 
			  const float pixWidth, 
			  const float pixHeight, 
			  const MGLsize width, 
			  const MGLsize height)
{	
	MGLfloat lowestX = min(tri.a.position[0],min(tri.b.position[0], tri.c.position[0]));
	MGLfloat lowestY = min(tri.a.position[1],min(tri.b.position[1], tri.c.position[1]));
	MGLfloat highestX = max(tri.a.position[0],max(tri.b.position[0], tri.c.position[0]));
	MGLfloat highestY = max(tri.a.position[1],max(tri.b.position[1], tri.c.position[1]));
	
	int starti = max((lowestX + 1) / pixWidth, 0.0f);
	int startj = max((lowestY + 1) / pixHeight, 0.0f);
	int endi = min((highestX + 1) / pixWidth, static_cast<float>(width-1));
	int endj = min((highestY + 1) / pixHeight, static_cast<float>(height-1));


	start = floor(width*startj + starti);
	end = ceil(width*endj + endi);
	
	return;
}

/**
 * Determines if a pixel in within the bounding box
 */
bool pointNotInBoundingBox(const unsigned pixel, 
			   const unsigned start, 
			   const unsigned end, 
			   const MGLsize width) 
{
	unsigned screenI = pixel % width;
	unsigned screenJ = pixel / width;
	
	return (screenI < (start % width) || screenI > (end % width) ||
			screenJ < (start / width) || screenJ > (end / width));
}

/**
 * Determines the area of a triangle using only vertices
 * Note: Multiplying by 1/2 is dropped due to this
 * function's use of determining barycentric coordinates via
 * ratios between triangle areas.
 */
float area(const vertex a, 
	   const vertex b, 
	   const vertex c)
{
	return (a.position[0] * (b.position[1] - c.position[1])) + 
	       (a.position[1] * (c.position[0] - b.position[0])) + 
	       ((b.position[0] * c.position[1]) - (b.position[1]*c.position[0]));
}

/**
 * A function that determine if a screen pixel coordinate is
 * inside the triangle in world space
 */
bool pointInTriangle(const triangle& tri,
		     const int i, 
		     const int j, 
		     const float pixWidth,
		     const float pixHeight) 
{
	float worldX = (i + 0.5)*pixWidth - 1;
	float worldY = (j + 0.5)*pixHeight - 1;
	
	vec4 pointPos = {worldX, worldY, 0, 1};
	vertex point;
	point.position = pointPos;
	point.color = currentColor;
	
	float areaOfTriangle = area(tri.a,tri.b,tri.c);
	
	float alpha = area(point, tri.b, tri.c) / areaOfTriangle;
	float beta = area(tri.a, point,tri.c) / areaOfTriangle;
	float gamma = area(tri.a, tri.b, point) / areaOfTriangle;
	
	if(alpha < 0.0f || beta < 0.0f || gamma < 0.0f) {
		return false;
	}

	return true;
}
 
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{	
	float pixWidth = 2.0 / width;
	float pixHeight = 2.0 / height;
	
	unsigned startpix = 0;
	unsigned endpix = 0;
	
	for(vector<triangle>::iterator t = listOfTriangles.begin(); t != listOfTriangles.end(); t++) {
		determineBoundingBox(*t, startpix, endpix, pixWidth, pixHeight, width, height);		
		for(unsigned pixel = startpix; pixel < endpix; pixel++) {
			if(pointNotInBoundingBox(pixel, startpix, endpix, width)) {
				//data[pixel] = Make_Pixel(0,0,0); // make black
				continue;
			}
			currentColor = t->a.color;
			
			int i = pixel % width;
			int j = pixel / width;
			
			if(pointInTriangle(*t, i, j, pixWidth, pixHeight)) {
				data[pixel] = Make_Pixel(currentColor[0], currentColor[1], currentColor[2]);
			}
		}
	}
}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
	if(mode == MGL_TRIANGLES) {
		drawmode = true;
	} else {
		drawmode = false;
	}
}


/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{	
	if(drawmode) { // handles triangles
		vertex coords[3];
		for(unsigned int i = 0; i < listOfVertices.size(); i++) {
			triangle newTri;
			for(unsigned j = 0; j < 3; j++, i++) {
				if(i >= listOfVertices.size()) {
					goto skip;
				}
				coords[j] = listOfVertices.at(i);
			}
			newTri.a = coords[0];
			newTri.b = coords[1];
			newTri.c = coords[2];
			listOfTriangles.push_back(newTri);
		}
	} else { // handles quads
		vertex coords[4];
		for(unsigned i = 0; i < listOfVertices.size(); i++) {
			triangle newTri1;
			triangle newTri2;
			for(unsigned j = 0; j < 4; j++, i++) {
				if(i >= listOfVertices.size()) {
					goto skip;
				}
				coords[j] = listOfVertices.at(i);
			}
			newTri1.a = coords[0];
			newTri1.b = coords[1];
			newTri1.c = coords[2];
			newTri2.a = coords[0];
			newTri2.b = coords[2];
			newTri2.c = coords[3];
			listOfTriangles.push_back(newTri1);
			listOfTriangles.push_back(newTri2);
		}
	}		

	skip:
	listOfVertices.clear();
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
	mglVertex3(x,y,0.0f);
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
	vec4 position = {x,y,z,1};

	position = currentProjMatrix * position;
	position = currentModelMatrix * position;
	
	vertex newVertex;
	
	newVertex.position = position;
	newVertex.color = currentColor;

	listOfVertices.push_back(newVertex);
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
	if(mode == MGL_PROJECTION) {
		matmode = true;
	} else {
		matmode = false;
	}
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
	if(matmode) {
		matrixStackProj.push_back(currentProjMatrix);
	} else {
		matrixStackModel.push_back(currentModelMatrix);
	}
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
	if(matmode) {
		if(!matrixStackProj.empty()) {
			matrixStackProj.erase(matrixStackProj.end());
		}
	} else {
		if(!matrixStackProj.empty()) {
			matrixStackModel.erase(matrixStackModel.end());
		}
	}
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
	if(matmode) {
		currentProjMatrix = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
	} else {
		currentModelMatrix = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
	}
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
	if(matmode) {
		currentProjMatrix = {matrix[0],matrix[4],matrix[8],matrix[12],
						 matrix[1],matrix[5],matrix[9],matrix[13],
						 matrix[2],matrix[6],matrix[10],matrix[14],
						 matrix[3],matrix[7],matrix[11],matrix[15]};
	} else {
		currentModelMatrix = {matrix[0],matrix[4],matrix[8],matrix[12],
						 matrix[1],matrix[5],matrix[9],matrix[13],
						 matrix[2],matrix[6],matrix[10],matrix[14],
						 matrix[3],matrix[7],matrix[11],matrix[15]};
	}
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
	mat4 multMatrix = {matrix[0],matrix[1],matrix[2],matrix[3],
					   matrix[4],matrix[5],matrix[6],matrix[7],
					   matrix[8],matrix[9],matrix[10],matrix[11],
					   matrix[12],matrix[13],matrix[14],matrix[15]};
    if(matmode) {					   
		currentProjMatrix = currentProjMatrix * multMatrix;
	} else {
		currentModelMatrix = currentModelMatrix * multMatrix;
	}
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
	//scale
	MGLfloat scaleX = 2 / (right - left);
	MGLfloat scaleY = 2 / (top - bottom);
	MGLfloat scaleZ = -2 / (far - near);
	
	//translate
	MGLfloat transX = -(right + left) / (right - left);
	MGLfloat transY = -(top + bottom) / (top - bottom);
	MGLfloat transZ = -(far + near) / (far - near);
	
	if(matmode) {
		currentProjMatrix = {scaleX, currentProjMatrix(0,1), currentProjMatrix(0,2), currentProjMatrix(0,3),
						     currentProjMatrix(1,0), scaleY, currentProjMatrix(1,2), currentProjMatrix(1,3),
						     currentProjMatrix(2,0), currentProjMatrix(2, 1), scaleZ, currentProjMatrix(2,3),
							 transX, transY, transZ, currentProjMatrix(3,3)};
	} else {
		currentModelMatrix = {scaleX, currentModelMatrix(0,1), currentModelMatrix(0,2), currentModelMatrix(0,3),
							  currentModelMatrix(1,0), scaleY, currentModelMatrix(1,2), currentModelMatrix(1,3),
							  currentModelMatrix(2,0), currentModelMatrix(2, 1), scaleZ, currentModelMatrix(2,3),
							  transX, transY, transZ, currentModelMatrix(3,3)};
	}
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
	currentColor = {red*255, green*255, blue*255};
}
