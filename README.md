# gravity-discs
A simple rigid body physics engine for entertaining eccentrically inclined youths. It is fancier than the usual physics engine in two ways:
* The engine parallelizes the motion of discs and their interactions using rayon
* The integration time-step is dynamically updated so that the end of each time-step occurs at or prior to the next elastic (or inelastic) collision if particles are assumed to move at constant velocity over each time-step.
  This is done by solving the quadratic equation that describes when two circles traveling at constant velocity overlap.
  - This assumption is exactly true when Euler's Method is used, so collision handling does not proliferate error beyond the error due to floating point arithmetic when Euler's Method is used. 
  - This assumption is only approximately true when Velocity Verlet due to the presence of acceleration. In principle, this could still be solved exactly. However, the equation that describes when two circles traveling at constant acceleration overlap is quartic and Rust's current packages for quartic equation solving are limited.

## Features for Users
Below is a discussion of the features that gravity-discs offers to users. 
Interaction with gravity-discs is mediated by mouse-clicks and keyboard buttons. 
Much of the information below is shown in gravity-discs itself whenever the number of entities/discs being simulated is zero.

### Creating and Removing Discs
Discs are created at the location of the cursor when the left mouse button is clicked. 
Discs overlapping with the cursor are removed when the right mouse button is clicked. 
All discs are removed when <kbd>backspace</kbd> is pressed.
The number of discs being simulated is shown in the top left corner of the gravity-discs window after the word *entities*.

The size of the disc to be created (in diameter) is represented by the bar on the left side of the gravity-discs window, which is exactly the same length as the diameter, visually. 
It is also shown numerically in the top left corner of the gravity-discs window after the phrase *disc r*.

Disc diameter can be increased by pressing <kbd>w</kbd> and decreased by pressing <kbd>s</kbd> on the keyboard.
In gravity-discs, the diameter can be "negative", conferring negative mass to the disc. 
This does not impact collision mechanics but does impact gravity: a disc with positive mass and a disc with negative mass will repel each other.

### Coefficient of Restitution
The [Coefficient of Restitution](https://en.wikipedia.org/wiki/Coefficient_of_restitution) is represented by the bar on the bottom of the gravity-discs window as well as by the percentage in the bottom-left corner. 
It can be increased by pressing <kbd>d</kbd> and decreased by pressing <kbd>a</kbd>. 
The smaller the coefficient of restitution, the faster discs will lose energy upon collision and the quicker they will become trapped in each other's gravitational wells -- 
forming planets.

Note that even though two discs may be trapped in contact with each other due to gravity, coefficient of restitution does not impact tangential velocity (this purpose is served by [friction](https://en.wikipedia.org/wiki/Friction), which is unimplemented).
As a result, discs in contact with each other may spin around each other in perpetuity with relatively preserved tangential motion -- giving the appearance that one is rolling around the other.
In actuality, this is not the case: one is *sliding* over the other. Angular velocity and momentum transfer are not implemented.

### Integration Strategy
The integration method and associated time-step are displayed in the top-center of the gravity-discs window.

As of now, [Euler's Method](https://en.wikipedia.org/wiki/Euler_method) and [Velocity Verlet](https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet) are implemented as possible integration methods. 
They can be toggled between by pressing <kbd>e</kbd> on the keyboard. The time-step (*dt*) can be increased by pressing <kbd>Shift</kbd> + <kbd>d</kbd> and decreased by pressing <kbd>Shift</kbd> + <kbd>a</kbd>. 

Note that, while Velocity Verlet solves the equations of motion more accurately the Euler's Method, the current collision handling implementation is only perfect when Euler's Method is used.

### Simulation Speed and Pausing the Simulation
The frames per second of the simulation can be seen in the to-left corner of the gravity-discs window after the phrase *FPS*.

The simulation can be paused by pressing <kbd>p</kbd> on the keyboard and unpaused by pressing it again. 
The paused/unpaused status is shown in the top-right corner of the gravity-discs window.
