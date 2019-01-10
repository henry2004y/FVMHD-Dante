# Dante
3D finite volume ideal MHD solver on structured grid.

## Getting Started
This code solves ideal MHD equations on a Cartesian mesh in 1/2/3D. It has the option of solving the energy equation (conservative form) or the pressure equation (nonconservative form). Currently only 1st and 2nd order face values are implemented. Rusanov and HLLE fluxes are implemented for face flux calculation.

One highlight feature is that there are no diverse branches for different dimensions: no additional "if" statement for 3D compared with 1D cases. The code is also organized in OOP, allowing potentially flexible extensions.

To get started, type

```
Main

```

All the parameters are set in 'Parameters.m'. The initial condition setup needs to be improved, so that we don't need to modify inside State.m to switch between Riemann problems and custom tests.

Note that moments instead of velocities are stored in state_GV. I may need to change the index (Ux_,Uy_,Uz_) to (RhoUx_,RhoUy_,RhoUz_) for clarification. Also, there are no DivB control besides the eight-wave scheme.

GPU-enabled part is easily implemented in Matlab. For 100000 cells in 1D, this gives roughly 10 times speedup.

### Prerequisites

MATLAB, Parallel Computing Toolbox, GPU-enabled driver (if using gpuArray).

## Authors

* **Hongyang Zhou** - *Initial work* - [henry2004y](https://github.com/henry2004y)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* All the nice guys who share their codes


