This directory contains the polynomial-based closed-form codes.

There are three approaches: 
1. Direct sum
2. When roots are clustered, apply a contour integration-based approach.
3. Taylor expansion about the center of the clustered roots.

Note: by setting a very small distance definition for clustering, Method 2 becomes Method 1, as the contour approach is never activated.

Taylor expansion-based approach is the most accurate. 

The Mathematica file computes both the derivative and the potential. The results are plotted using different MATLAB codes.
