# CSCI420-Computer-Graphics-Ray-Tracing

Problem Statement: http://run.usc.edu/cs420-s18/assignments/assign3/assign3.html

In this assignment, you will be building a ray tracer. Your ray tracer will be able to handle opaque surfaces with lighting and shadows. Provided for you will be starter code that will load scene data from a file.

Step 1: Uniformly send out rays from the camera location. Since the camera does not have to move, you can assume that its location is (0,0,0). You should use backwards ray tracing where rays are sent from the camera, one ray per pixel. The final images should be 640x480, but for debugging you should use smaller resolutions with faster rendering times. For example, if you halve each dimension, you would send out 1/4th of the number of rays. You can use the field of view of 60 degrees.

Step 2: Write the intersection code. The mathematical solutions for the intersection code are provided in the lecture notes.

Step 3: Implement the illumination equations (provided below).

Step 4: Create still images showing off your ray tracer.

The source code gives you two methods of plotting a pixel. You can either only plot to the screen, or both to the screen and a JPEG file. You control this behavior by providing the second command line argument to your program. When the second argument is provided, output is generated both to screen and the filename given by the second argument. If you do not provide a second argument, plotting will be to screen only.

Functionality Requirements
This is the list of requirements for this assignment:
Triangle intersection (20 points)
Sphere intersection (20 points)
Triangle Phong shading (15 points)
Sphere Phong shading (15 points)
Shadows rays (15 points)
Still images (15 points)

Extra credit implemented:
Good antialiasing (10 points)
- add a cmd arg (after scene file arg) to enable or disable antialiasing
./hw3 test.scene true test.jpg
- the color is calculated after taking an average of the colors calculated after shooting 4 rays to the point
