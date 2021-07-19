# LBM_2D_Flow_Over_Obstacle
A Python code to simulate the flow over a user-defined obstacle using the Lattice Boltzmann Method

The code allows the user to simulate the flow over different obstacles - Circle, Square, Rotated Square, Ellipse and a symmetrical NACA Airfoil at a specified angle of attack. The user can change the obstacle deefinition in line 91 of the code.

Numpy and Matplotlib modules need to be installed. The images and files are created in the same folder as the code. Currently, the code produces the post-processing files in .dat format, which can be post-processed using Tecplot. However, files can also be generated in other formats such as CSV by changing the definition in lines 148 to 163 of the code.
