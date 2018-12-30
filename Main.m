% Finite Volume MHD solver: DANTE
% Designed with OOP, readability and efficiency.
%
% Hongyang Zhou, 05/31/2018
% First version finished in 06/06/2018
% First series of tests passed in 07/29/2018
% Modified 12/30/2018, for fully OOP and better data structures.

clear;clc; close all
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
state_GV = state.SetState;
state.plot('rho',0)

boundary  = Boundary;
faceValue = FaceValue;
faceFlux  = FaceFlux(faceValue);
source    = Source;

clearvars density velocity Bfield pressure

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
   state_GV = boundary.set_cell_boundary(state_GV);
   
   % Calculate face value
   faceValue = faceValue.calc_face_value(state_GV);
   
   % Calculate face flux
   faceFlux = faceFlux.calc_face_flux;
   
   % Calculate source
   source.calc_source(grid,state_GV);
   
   % Calculate time step
   dt = calc_timestep(grid,state_GV);
   if t+dt>tEnd; dt=tEnd-t; end
      
   % Update state
   state_GV = update_state(grid,state_GV,faceFlux,source,dt); 
   
   t = t + dt;
   it = it + 1;
   
   state = state.GetState(state_GV);
   state.plot('rho',it)
   %state.plot('ux',it)
   pause(.1)
   
   fprintf('it,t=%d,%f\n',it,t)
end

   case 2
while t < tEnd % 2nd order method
   % 1st stage of modified timestepping
   
   % Set boundary conditions
   state_GV = set_cell_boundary(state_GV);
   
   % Calculate face value
   faceValue = faceValue.calc_face_value(state_GV);
   
   % Calculate face flux
   faceFlux = faceFlux.calc_face_flux(faceValue);
   
   % Calculate source
   source.calc_source(grid,state_GV);
   
   % Calculate time step
   dt = calc_timestep(grid,state_GV);
   if t+dt>tEnd; dt=tEnd-t; end
   
   state1_VG = update_state(grid,state_GV,faceFlux,source,0.5*dt);
   
   % 2nd stage of modified timestepping
   
   % Set boundary conditions
   state1_VG = set_cell_boundary(state1_VG);
   
   % Calculate face value
   faceValue = faceValue.calc_face_value(state1_VG);
   
   % Calculate face flux
   faceFlux = faceFlux.calc_face_flux(faceValue);
   
   % Calculate source
   source_VG = calc_source(grid,state1_VG);
   
 
      
   % Update state
   state_GV = update_state(grid,state_GV,faceFlux,source_VG,dt);
   
   
   t = t + dt;
   it = it + 1;
   
   state = state.GetState(state_GV);
   state.plot('rho',it)
   pause(.1)
   
   fprintf('it,t=%d,%f\n',it,t)
end

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

% state = state.GetState(state_VG);
% 
% state.plot('rho',iStep)

