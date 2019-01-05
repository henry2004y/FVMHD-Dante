% Finite Volume MHD solver: DANTE
% Designed with OOP, readability and efficiency.
%
% Hongyang Zhou, hyzhou@umich.edu 05/31/2018
% First version finished in 06/06/2018
% First series of tests passed in 07/29/2018
% Modified 12/30/2018, for fully OOP and better data structures.

clear; clc;% close all
%% Initialization

% Constants set in class NamedConst
% Input parameters set in class Parameters

disp('3D Finite Volume MHD Simulation')
disp('Simulation starts...')

% Display parameters
%printParameters
details(Parameters)

% Init grid
grid = Grid;

% Init variables and arrays
state = State(false,grid);
state = state.SetState;
state.plot('rho',grid,0)

boundary  = Boundary;
faceValue = FaceValue;
faceFlux  = FaceFlux;
source    = Source;
time      = Time;

disp('Initialization finished...')

%% Advance

% while 1
%    DoStop = stop_condition_true;
%    if DoStop
%       break
%    end
% end

% Set T if not already set by set_init
if ~exist('tEnd','var'); tEnd = Parameters.tEnd; end

t = 0; it = 0;
tic
switch Parameters.Order
   case 1
      while t < tEnd % 1st order method
         % Set boundary conditions
         state = boundary.set_cell_boundary(state);
         
         % Calculate face value
         faceValue.calc_face_value(state);
         
         % Calculate face flux
         faceFlux.calc_face_flux(faceValue);
         
         % Calculate source
         source.calc_source(grid,state);
         
         % Calculate time step
         time.calc_timestep(grid,state);
         
         if t + time.dt > tEnd; time.dt = tEnd - t; end
         
         % Update state
         state = state.update_state(grid,faceFlux,source,time);
         
         t  = t + time.dt;
         it = it + 1;
         
         state = state.GetState;
         state.plot('rho',grid,it)
         %state.plot('ux',grid,it)
         %pause(.05)
         
         fprintf('it,t=%d,%f\n',it,t)
      end
      
   case 2
      while t < tEnd % 2nd order method
         % 1st stage of modified timestepping
         
         % Set boundary conditions
         state = boundary.set_cell_boundary(state);
         
         % Calculate face value
         faceValue.calc_face_value(state);
         
         % Calculate face flux
         faceFlux.calc_face_flux(faceValue);
         
         % Calculate source
         source.calc_source(grid,state);
         
         % Calculate time step
         time.calc_timestep(grid,state);
         
         % Update state in the 1st stage
         time.dt = 0.5*time.dt;
         state1 = state.update_state(grid,faceFlux,source,time);
         
         
         % 2nd stage of modified timestepping
         
         % Set boundary conditions
         state1 = boundary.set_cell_boundary(state1);
         
         % Calculate face value
         faceValue.calc_face_value(state1);
         
         % Calculate face flux
         faceFlux.calc_face_flux(faceValue);
         
         % Calculate source
         source.calc_source(grid,state1);
         
         % Update state
         time.dt = 2*time.dt;
         state = state.update_state(grid,faceFlux,source,time);
         
         if t + time.dt > tEnd; time.dt = tEnd - t; end
         
         t = t + time.dt;
         it = it + 1;
         
         state = state.GetState;
         state.plot('rho',grid,it)
         pause(.05)
         
         fprintf('it,t=%d,%f\n',it,t)
      end
      clearvars state1
   otherwise
      error('Higher Order schemes not yet implemented!')
end
cputime = toc;

% for iStep = 1:Parameters.nStep
%
%    for iStage = 1:Parameters.nStage
%       % Set boundary conditions
%       state_VG = set_cell_boundary(state_VG);
%
%       % Calculate face value
%       %faceValue = calc_face_value(faceValue,state);
%       faceValue = faceValue.calc_face_value(state_VG);
%
%       % Calculate face flux
%       faceFlux = faceFlux.calc_face_flux(faceValue);
%
%       % Calculate source
%       source_VG = calc_source(grid,state_VG);
%
%       % Calculate time step
%       dt = calc_timestep(grid,state_VG);
%       %timestep = calc_timestep(grid,faceValue);
%
%       % Update state
%       state_VG = update_state(grid,state_VG,faceFlux,source_VG,dt);
%
%    end
%
%    %state_VG(:,3:3,2,2)
%
%    state = state.GetState(state_VG);
%    state.plot('rho',iStep)
%    pause(.1)
%
%    fprintf('iStep=%d\n',iStep)
%
% end

disp('advance finished...')

%% Visualization


if strcmp(Parameters.IC,'Riemann')
   nI = Parameters.nI;
   
   % Exact solution
   [xe,re,ue,pe,ee,te,Me,se] = ...
      EulerExact(1,0,1, 0.125,0,0.1,0.1, 3);
   Ee = pe./((Const.gamma-1)*re)+0.5*ue.^2;
   
   figure(2); hold on
   x = grid.getX;
   plot(xe,re,'--')
end
