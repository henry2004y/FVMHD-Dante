classdef Parameters
   %Main The definition of general information in the simulation.
   %
   % GridType  Type of simulation grid. Options:
   %    (1) Cartesian  (2) Spherical
   % GridSize  Size of simulation grid
   % nG        Number of ghost cells
   
   properties (Constant)
      GridType   char = 'Cartesian'
      GridSize   double {mustBeInteger} = [100,1,1]
      nI         double {mustBeInteger} = Parameters.GridSize(1)
      nJ         double {mustBeInteger} = Parameters.GridSize(2)
      nK         double {mustBeInteger} = Parameters.GridSize(3)
      nG         double {mustBeInteger} = 1
      xyzMinMax(3,2) double {mustBeReal}= [0 1;0 1; 0 1]
      % Size including ghost cells
      FullSize   double {mustBeInteger} = ...
         Parameters.GridSize + 2*Parameters.nG
      % Indexes for physical cells
      iMin       double {mustBeInteger} = 1 + Parameters.nG
      iMax       double {mustBeInteger} = Parameters.nI + Parameters.nG
      jMin       double {mustBeInteger} = 1 + Parameters.nG
      jMax       double {mustBeInteger} = Parameters.nJ + Parameters.nG
      kMin       double {mustBeInteger} = 1 + Parameters.nG
      kMax       double {mustBeInteger} = Parameters.nK + Parameters.nG
      % Indexes for all cells including ghost cells
      iMinAll    double {mustBeInteger} = 1
      iMaxAll    double {mustBeInteger} = Parameters.nI + 2*Parameters.nG
      jMinAll    double {mustBeInteger} = 1
      jMaxAll    double {mustBeInteger} = Parameters.nJ + 2*Parameters.nG
      kMinAll    double {mustBeInteger} = 1
      kMaxAll    double {mustBeInteger} = Parameters.nK + 2*Parameters.nG
      
      
      Scheme     char {mustBeMember(Scheme,{'Rusanov','HLLE'})}= 'Rusanov'
      Order      double {mustBeMember(Order,[1,2])} = 1
      CFL        double = 0.9
      limiter    char = 'MM'
      TimeAccurate logical = true
      UseConservative logical = false
      nStage     double {mustBeInteger} = Parameters.Order
      nVar       double {mustBeInteger} = 8
      Rho_       = 1
      Ux_        = 2
      Uy_        = 3
      Uz_        = 4
      Bx_        = 5
      By_        = 6
      Bz_        = 7
      P_         = 8
      E_         = Parameters.P_
      U_         = [Parameters.Ux_ Parameters.Uy_ Parameters.Uz_]
      B_         = [Parameters.Bx_ Parameters.By_ Parameters.Bz_]
      
      BC         char = 'periodic'
      IC         char {mustBeMember(IC,{'density wave','square wave', ...
         'contact discontinuity','shocktube','Riemann'})}= 'density wave'
      RiemannProblemType double {mustBeInteger} = 1;
      
      DoAdvanceTime logical = false
      nStep      double {mustBeInteger} = 20
      tEnd       double {mustBeGreaterThan(tEnd, 0)} = 0.1
      
      PlotVar    char = 'rho'
   end
   
end