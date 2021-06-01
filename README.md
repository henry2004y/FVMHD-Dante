# Dante
3D finite volume ideal MHD solver on structured grid.

**NOTE**
This code has been completely rewritten in [Julia](https://github.com/henry2004y/Dante). Future versions will be maintained there.

Practice is more important than theory. Dante is a toy 3D finite-volume code for learning. This first version is implemented in MATLAB.
Some main features of the code include:
* Use primitive variables for visualization and conservative variables for computation;
* Assimilate the idea of object and class for defining variables and function checks;
* Apply to 1D/2D/3D Cartesian grid with the same kernel code;
* 1st/2nd order Rusanov/HLLE scheme with forwar-Euler/modified-Euler explicit timestepping;
* User-friendly input parameter control;
* Optimized matrix operation, avoid nested loops as much as possible;
* One logical CPU/GPU switch.

This code solves ideal MHD equations on a Cartesian mesh in 1/2/3D. It has the option of solving the energy equation (conservative form) or the pressure equation (nonconservative form). Currently only 1st and 2nd order face values are implemented. Rusanov and HLLE fluxes are implemented for face flux calculation.

One highlight feature is that there are no diverse branches for different dimensions: no additional "if" statement for 3D compared with 1D cases. The code is also organized in OOP, allowing potentially flexible extensions.

## Getting Started

To get started, type

```
Main
```

All the parameters are set in [`Parameters.m`](Parameters.m). The initial condition setup needs to be improved, so that we don't need to modify inside [`State.m`](State.m) to switch between Riemann problems and custom tests.

Note that moments instead of velocities are stored in state_GV. I may need to change the index (Ux_,Uy_,Uz_) to (RhoUx_,RhoUy_,RhoUz_) for clarification. Also, there are no DivB control besides the eight-wave scheme. General Lagrange multiplier is for sure something to try.

GPU-enabled part is easily implemented in MATLAB. For 100,000 cells in 1D, this gives roughly 10 times speedup using GTX1080 compared with 4 core intel i7 6700k.

### Prerequisites

MATLAB, Parallel Computing Toolbox, GPU-enabled driver (if using gpuArray).

## Authors

* **Hongyang Zhou** - *Initial work* - [henry2004y](https://github.com/henry2004y)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* All the nice guys who share their codes


