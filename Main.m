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

if Parameters.UseGPU
   numDevices = gpuDeviceCount
   origDevice = gpuDevice
end

% Init grid
grid = Grid;

% Init variables and arrays
state = State(false,grid);
state = state.SetState;
state.plot(Parameters.PlotVar,grid,0)

boundary  = Boundary;
faceValue = FaceValue;
faceFlux  = FaceFlux;
source    = Source;
time      = Time;

t = 0; it = 0;
% Set T if not already set by set_init
if ~exist('tEnd','var'); tEnd = Parameters.tEnd; end

disp('Initialization finished...')

%% Advance

tic
if Parameters.DoAdvanceTime % Advance with time
   switch Parameters.Order
      case 1 % 1st order method
         while t < tEnd
            % Set boundary conditions
            state = boundary.set_cell_boundary(state);
            
            % Calculate face value
            faceValue.calc_face_value(state);
            
            % Calculate face flux
            faceFlux.calc_face_flux(faceValue);
            
            % Calculate source
            source.calc_source(grid,state);
            
            % Calculate time step
            time.calc_timestep(grid,faceFlux);
            
            if t + time.dt > tEnd; time.dt = tEnd - t; end
            
            % Update state
            state = state.update_state(grid,faceFlux,source,time);
            
            t  = t + time.dt;
            it = it + 1;
            
            if mod(it,Parameters.PlotInterval) == 0
               state.plot(Parameters.PlotVar,grid,it)
            end
            
            fprintf('it,t=%d,%f\n',it,t)
         end
         
      case 2 % 2nd order method
         while t < tEnd
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
            time.calc_timestep(grid,faceFlux);
            
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
            
            if mod(it,Parameters.PlotInterval) == 0
               state.plot(Parameters.PlotVar,grid,it)
            end
            
            fprintf('it,t=%d,%f\n',it,t)
         end
         clearvars state1
      otherwise
         error('Higher Order schemes not yet implemented!')
   end
   state.plot(Parameters.PlotVar,grid,it)
   
else % Advance with steps
   switch Parameters.Order
      case 1 % 1st order method
         for iStep = 1:Parameters.nStep
            % Set boundary conditions
            state = boundary.set_cell_boundary(state);
            
            % Calculate face value
            faceValue.calc_face_value(state);
            
            % Calculate face flux
            faceFlux.calc_face_flux(faceValue);
            
            % Calculate source
            source.calc_source(grid,state);
            
            % Calculate time step
            time.calc_timestep(grid,faceFlux);
            
            % Update state
            state = state.update_state(grid,faceFlux,source,time);
            
            t  = t + time.dt;
            
            if mod(iStep,Parameters.PlotInterval) == 0
               state.plot(Parameters.PlotVar,grid,iStep)
            end
            
            fprintf('it,t=%d,%f\n',iStep,t)
         end
         
      case 2 % 2nd order method
         for iStep = 1:Parameters.nStep
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
            time.calc_timestep(grid,faceFlux);
            
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
            
            t = t + time.dt;
            
            if mod(iStep,Parameters.PlotInterval) == 0
               state.plot(Parameters.PlotVar,grid,iStep)
            end
            
            fprintf('it,t=%d,%f\n',iStep,t)
         end
         clearvars state1
      otherwise
         error('Higher Order schemes not yet implemented!')
   end
   
   state.plot(Parameters.PlotVar,grid,iStep)
end
cputime = toc;
disp('advance finished...')

%% Visualization

if Parameters.UseGPU, t = gather(t); end

if strcmp(Parameters.IC,'Riemann')
   nI = Parameters.nI;
   [Rho,U,~,P,tEnd,~] = state.set_init_Riemann(grid);
   % Exact solution
   [xe,re,ue,pe,ee,te,Me,se] = ...
      EulerExact(Rho(1),U(1),P(1), ...
      Rho(end),U(end,1,1,1),P(end),t, 3);
   
   Ee = pe./((Const.gamma-1)*re)+0.5*ue.^2;
   
   clearvars Rho U B P CFL
   figure(2); hold on
   x = grid.getX;
   switch Parameters.PlotVar
      case 'rho'
         plot(xe,re,'--')
      case 'p'
         plot(xe,pe,'--')
      case 'ux'
         plot(xe,ue,'--')
   end
end
