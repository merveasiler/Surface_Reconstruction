Implemented the paper "[Surface Reconstruction from Unorganized Points.](https://dl.acm.org/doi/pdf/10.1145/133994.134011)" [1]

I have tried many things in this task, and I think I obtained acceptable results:

- I calculated vertex normals by getting the average of neighbor triangle normals. I also tried the weighted average of those normals by taking the inverse of area of the triangle as weight but did not change the accuracy much..
  
- I constructed my grid and for each vertex of it I checked whether it is inside the object or outside the object. To do this, I used sdf. For each vertex of the grid I found the nearest vertex of the mesh and checked the cosine of the angle between the normal of mesh vertex and the vector from the mesh vertex to the grid vertex. If this values is positive, then it means the grid vertex is outside the object, if it is negative, then the grid vertex is inside the object. 

- To find the nearest vertex, initially I used kd-tree. It was working actually once upon a time, but later on (I don't know why) it was corrupted somehow and I could not repair it on its time. Therefore, I changed my detecting nearest point implementation into private one-by-one control.

- As a result, I obtained valid results if we neglect the noise. To avoid from the noise, in the paper it shares some test results and as much as I understand there does not exist a generic rule to avoid from noise (okey it says "compare the noise+density value with the distance from the projected grid vertex on the tangent plane of the nearest mesh vertex to the mesh vertex itself.", but I tried, it does not work). Previously, I had constructed my grid such that each cell edge has a length of the 2% of the whole grid edge (the grid starts from the minimum x,y,z values of vertices and ends at maximum x,y,z values). Thus, I have decided to use some multiple of this cell edge length as a threshold to avoid from noise. I compared this threshold with the distance of nearest vertex to grid vertex and if this distance is bigger than my threshold, I assumed that grid vertex is outside the object even if its sdf is negative.

- Hence there remains only one question, what will we choose that multiple (I mean the multiple of cell edge length)? I tried some coefficients, and I obtained the mostly optimized results with "1.7" for most of the inputs.

- You can find my screenshots in the folder, and you can simply build my implementation by typing \<make\> and run by typing \<./Main\>

- It will ask the input file name, and then the coefficient that I mentioned above. In this way, you can also try different thresholds by giving different coefficients other than 1.7.

- Also you can see the measured time in the execution.


[1] Hoppe, H., DeRose, T., Duchamp, T., McDonald, J., & Stuetzle, W. (1992, July). Surface reconstruction from unorganized points. In Proceedings of the 19th annual conference on computer graphics and interactive techniques (pp. 71-78).


