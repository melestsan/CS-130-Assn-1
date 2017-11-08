# CS130 Fall 2017 Assignment 1i
======
Implementation of a simplified 3D rendering pipeline with flat shading.
1. Vertex and viewing transformations
2. Rasterization and interpolation
3. Clipping
4. Using a z-buffer for hidden surfaces
## Sofware Framebuffer
------
A software framebuffer that stores both a 2D array of pixel colors and a 2D array of z(depth) values. Rasterization routine will write into this software buffer.
## Rasterization and z-buffer depth test
------
Implement a routine that draws a filled triangled specified by three vertices into the software framebuffer. The routine should write its result to the software framebuffer, doing the appropriate z-buffer check.
## Vertex and viewing transformations
------
Similar to OpenGL, both a projection matrix stack and a modelview stack will be maintained. When the user draws a vertex, modelview and projection will be applied to obtain the transformed geometry. Followed by a divide by w, and a viewport transform. The transformed geometry will be stored for rasterization. The transformations are part of the current state and will not necessarily persist.
