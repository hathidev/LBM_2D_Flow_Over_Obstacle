# LBM_2D_Flow_Over_Obstacle
A Python code to simulate the flow over a user-defined obstacle using the Lattice Boltzmann Method

The collision step uses the Bhatnagar-Gross-Krook(BGK) Model and a single time-relaxation parameter, that is calculated using the Reynolds Number. The drag and lift forces are calculated using the Momentum Exchange method applied at Bounce-back boundaries.

The code allows the user to simulate the flow over different obstacles - Circle, Square, Rotated Square, Ellipse and a symmetrical NACA Airfoil at a specified angle of attack and at a given Reynolds Number. The user can change the Reynolds Number and the obstacle definition in lines 12 and 158 of the code respectively.

Numpy and Matplotlib modules need to be installed. The images and files are created in the same folder as the code. Currently, the code produces the post-processing files in .dat format, which can be post-processed using Tecplot. However, files can also be generated in other formats such as CSV by changing the definition in lines 134 to 142 of the code.
