# Dante
3D finite volume ideal MHD solver on structured grid.

## Getting Started
This code solves ideal MHD equations on a Cartesian mesh in 1/2/3D. It has the option of solving the energy equation (conservative form) or the pressure equation (nonconservative form). Currently only 1st and 2nd order face values and Rusanov flux is implemented.

One highlight feature is that there are no diverse branches for different dimensions: no additional "if" statement for 3D compared with 1D cases. The code is also organized in OOP, allowing potentially flexible extensions.

To get started, type

```
Main

```

All the parameters are set in 'Parameters.m'. The initial condition setup needs to be improved, so that we don't need to modify inside State.m to switch between Riemann problems and custom tests.

### Prerequisites

MATLAB

## Authors

* **Hongyang Zhou** - *Initial work* - [henry2004y](https://github.com/henry2004y)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* All the nice guys who share their codes


