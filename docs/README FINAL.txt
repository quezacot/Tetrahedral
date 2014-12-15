Relevant Files:
Simulations are ran from simulation.cpp.

One ball simulation:
Simulation Overview:
In this simulation, we simulate a single ball bouncing inside an impenetrable box. Users can interact with the ball via increasing or decreasing wind force and by dragging or flinging the ball inside the box.
Make sure #define MESH2 is commented out in simulation.cpp. Otherwise, this will simulate 2 balls. 

Command Line Code Example:
make simulation
./simulation data/sphere87.*
./simulation data/sphere676.*

Functionalities:
1.	Rotate camera using keyboard arrow keys.
a.	http://youtu.be/gSsR8oq4tdA

2.	Reset z-axis as the vertical axis using keyboard button ¡§L.¡¨
a.	http://youtu.be/gSsR8oq4tdA at 0:25-0:28

3.	Wind force manipulations
a.	¡§D¡¨ to increase wind speed
b.	¡§A¡¨ to decrease wind speed
c.	¡§W¡¨ to resume wind force
d.	¡§S¡¨ to stop wind force
e.	http://youtu.be/O4mwjM8yCU8 (Look at the command line output to see commands)

4.	Use mouse left-click to hold an object (make sure the camera is close to the object)
a.	http://youtu.be/rJStBwuFs_A (Look at the command line output when the ball is caught)

5.	While holding an object, make a swipe action to drag or fling the object.
a.	http://youtu.be/rJStBwuFs_A 
?
Two balls simulation:
Simulation Overview:
In this simulation, we simulate two balls bouncing inside an impenetrable box.
Make sure #define MESH2 is not commented out in simulation.cpp. Otherwise, this will simulate 1 ball. 

make simulation
./simulation data/sphere87.*

Functionalities:
Same as the above simulation of one ball.
http://youtu.be/ja5GfRd1OKI

Note: When the forces are too strong, there are some chance that balls become unstable and will get meshed up together. Also, we can only run the two ball simulation with sphere87 because we haven¡¦t found the correct values for spring constant, damping constant, and K in volume penalty to properly simulate bigger spheres. Please just restart the simulation when this occurs.
