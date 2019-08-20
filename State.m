classdef State %< handle
   %STATE The class of physical states in MHD.
   %   Primitive variables include:
   % rho: density
   % u:   velocity (3D)
   % B:   magnetic field (3D)
   % p:   pressure
   % e:   energy
   %
   % For higher order schemes, this should be a normal class instead of
   % handle class.
   
   %======================== MEMBERS =================================
   properties
      Rho(:,:,:)  %double {mustBeReal, mustBeGreaterThan(Rho, 0)}
      U(:,:,:,3)  %double {mustBeReal}
      B(:,:,:,3)  %double {mustBeReal}
      P(:,:,:)    %double {mustBeReal, mustBeGreaterThan(P, 0)}
      
      %state_GV(:,:,:,:) {typeCheck}? %double
      state_GV(:,:,:,:)
   end
   
   properties (Dependent)
      % Energy
      Energy double {mustBeReal, mustBeGreaterThan(Energy, 0)}
   end
   
   %======================== CONSTRUCTORS ============================
   methods
      function obj = State(UseDefault,grid)
         % Construct an instance of this class
         %
         
         disp('Initializing state variables...')
         
         % Set initial conditions
         switch Parameters.IC
            case {'density wave','contact discontinuity','shocktube',...
                  'square wave'}
               [density,velocity,Bfield,pressure] = obj.set_init;
            case 'Riemann'
               [density,velocity,Bfield,pressure,tEnd] = ...
                  obj.set_init_Riemann(grid);
            otherwise
               error('unknown initial condition!')
         end
         
         if nargin == 0
            disp('Using default grid')
            obj.Rho = ones(GridGen.GridSize);
            % May be I don't need initialization this way.
         elseif ~UseDefault
            % Call asset constructor
            obj.Rho = density;
            obj.U   = velocity;
            obj.B   = Bfield;
            obj.P   = pressure;
            obj = obj.SetState;
         end
      end
   end
   
   methods (Static)
      function [density,velocity,Bfield,pressure] = set_init
         %SET_INIT Initialization of state variables.
         %
         
         if Parameters.UseGPU
            density  = ones(Parameters.FullSize,'gpuArray');
            velocity = zeros([Parameters.FullSize,3],'gpuArray');
            Bfield   = zeros([Parameters.FullSize,3],'gpuArray');
            pressure = ones(Parameters.FullSize,'gpuArray');
         else
            density  = ones(Parameters.FullSize);
            velocity = zeros([Parameters.FullSize,3]);
            Bfield   = zeros([Parameters.FullSize,3]);
            pressure = ones(Parameters.FullSize);
         end
         
         nI = Parameters.nI;
         
         switch Parameters.IC
            case 'contact discontinuity'
               velocity(:,:,:,:) = 0.0;
               density(1:nI/2,:,:) = 2.0;
               density(nI/2+1:end,:,:) = 1.0;
            case 'density wave'
               density(nI/2:nI/2,:,:) = 2.0;
               velocity(:,:,:,1) = 1.0;
               pressure(:) = 0.01;
            case 'square wave'
               density(nI/2-10:nI/2,:,:) = 2.0;
               velocity(:,:,:,1) = 1.0;
               pressure(:) = 0.01;
            case 'shocktube'
               
            otherwise
               error('unknown initial condition type!')
         end
      end
      
      function [Rho,U,B,P,tEnd,CFL] = set_init_Riemann(grid)
         %SET_INIT_RIEMANN Load the IC of a classical 1D Riemann Problems.
         %   Detailed explanation goes here
         
         
         % Riemann Problems
         switch Parameters.RiemannProblemType
            case{1} % Configuration 1, Sod's Problem
               fprintf('Case 1: Sods problem \n');
               p   = [1    0.1  ];
               u   = [0    0    ];
               rho = [1    0.125];
               tEnd = 0.1; CFL = 0.90;
               
            case{2} % Configuration 2, Left Expansion and right strong shock
               fprintf('Case 2: Left Expansion and right strong shock \n');
               p   = [1000 0.1  ];
               u   = [0    0    ];
               rho = [3    2    ];
               tEnd = 0.02; CFL = 0.90;
               
            case{3} % Configuration 3, Right Expansion and left strong shock
               fprintf('Case 3: Right Expansion and left strong shock \n');
               p   = [7    10   ];
               u   = [0    0    ];
               rho = [1    1    ];
               tEnd = 0.1; CFL = 0.90;
               
            case{4} % Configuration 4, Double Shock
               fprintf('Case 4: Double Shock \n');
               p   = [450  45   ];
               u   = [20   -6   ];
               rho = [6    6    ];
               tEnd = 0.02; CFL = 0.90;
               
            case{5} % Configuration 5, Double Expansion
               fprintf('Case 5: Double Expansion \n');
               p   = [40   40   ];
               u   = [-2   2    ];
               rho = [1    2.5  ];
               tEnd = 0.03; CFL = 0.90;
               
            case{6} % Configuration 6, Cavitation
               fprintf('Case 6: Cavitation \n');
               p   = [0.4  0.4  ];
               u   = [-2    2   ];
               rho = [ 1    1   ];
               tEnd = 0.1; CFL = 0.90;
               
            case{7} % Shocktube problem of G.A. Sod, JCP 27:1, 1978
               fprintf('Shocktube problem of G.A. Sod, JCP 27:1, 1978');
               p   = [1.0  0.1  ];
               u   = [0.75 0    ];
               rho = [1    0.125];
               tEnd = 0.17; CFL = 0.90;
               
            case{8} % Lax test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997
               fprintf('Lax test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997');
               p   = [3.528 0.571];
               u   = [0.698 0    ];
               rho = [0.445 0.5  ];
               tEnd = 0.15; CFL = 0.90;
               
            case{9} % Mach = 3 test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997
               fprintf('Mach = 3 test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997');
               p   = [10.333  1  ];
               u   = [ 0.92  3.55];
               rho = [ 3.857  1  ];
               tEnd = 0.09; CFL = 0.90;
               
            case{10} % Shocktube problem with supersonic zone
               fprintf('Shocktube problem with supersonic zone');
               p   = [1  0.02];
               u   = [0  0.00];
               rho = [1  0.02];
               tEnd = 0.162; CFL = 0.90;
               
            case{11} % Contact discontinuity
               fprintf('Contact discontinuity');
               p   = [0.5 0.5];
               u   = [0.0 0.0];
               rho = [1.0 0.6];
               tEnd = 1; CFL = 0.90;
               
            case{12} % Stationary shock
               fprintf('Stationary shock');
               p   = [ 1.0  0.1 ];
               u   = [-2.0 -2.0 ];
               rho = [ 1.0 0.125];
               tEnd = 0.1; CFL = 0.28;
               
            case{13} % left side of 2-d Riemman case 12
               fprintf('left side of 2-d Riemman case 13');
               p   = [ 1.0  1.0 ];
               u   = [ 0.0  0.0 ];
               rho = [ 0.8  1.0 ];
               tEnd = 0.1; CFL = 0.25;
               
            case{14} % right side of 2-d Riemman case 12
               fprintf('right side of 2-d Riemman case 14');
               p   = [ 1.0  0.4 ];
               u   = [ 0.7276  0.0 ];
               rho = [ 1.0  0.5313 ];
               tEnd = 0.1; CFL = 0.25;
               
            case{15} % Stationary Shock
               fprintf('right side of 2-d Riemman case 15');
               p   = [ 0.1  0.676 ];
               u   = [ 1.2  0.723966942148760 ];
               rho = [ 1.0  1.657534246575342 ];
               tEnd = 0.1; CFL = 0.25;
               
            otherwise
               error('Case not available');
               
         end
         % Print for Riemann Problems
         fprintf('\n');
         fprintf('density (L): %1.3f\n',rho(1));
         fprintf('velocity(L): %1.3f\n',u(1));
         fprintf('Pressure(L): %1.3f\n',p(1));
         fprintf('\n');
         fprintf('density (R): %1.3f\n',rho(2));
         fprintf('velocity(R): %1.3f\n',u(2));
         fprintf('Pressure(R): %1.3f\n',p(2));
         fprintf('\n');
         
         % Load Selected case Initial condition:
         % Pre-Allocate variables
         Rho = ones(Parameters.FullSize);
         U   = zeros([Parameters.FullSize,3]);
         B   = zeros([Parameters.FullSize,3]);
         P   = ones(Parameters.FullSize);
         
         % Parameters of regions dimensions
         x_middle = 0.5*...
            (Parameters.xyzMinMax(1,2) - Parameters.xyzMinMax(1,1));
         L = find(grid.X<=x_middle);
         R = find(grid.X>x_middle);
         
         % Initial Condition for our 1D domain
         % Density
         Rho(L) = rho(1); % region 1
         Rho(R) = rho(2); % region 2
         % Velocity in x
         U(L) = u(1); % region 1
         U(R) = u(2); % region 2
         % Pressure
         P(L) = p(1); % region 1
         P(R) = p(2); % region 2
         
      end
      
   end
   
   %======================== METHODS =================================
   methods
      function obj = SetState(obj)
         % Reorganize data structure.
         
         if Parameters.UseGPU
            obj.state_GV = Inf([Parameters.FullSize,Parameters.nVar],...
               'gpuArray');
         else
            obj.state_GV = Inf([Parameters.FullSize,Parameters.nVar]);
         end
         
         obj.state_GV(:,:,:,Parameters.Rho_) = obj.Rho;
         obj.state_GV(:,:,:,Parameters.U_)   = obj.U .*...
            obj.state_GV(:,:,:,Parameters.Rho_);
         obj.state_GV(:,:,:,Parameters.B_) = obj.B;
         obj.state_GV(:,:,:,Parameters.P_) = obj.P;
      end
      
      function obj = GetState(obj)
         % Reorganize data structure for post-processing.
         
         obj.Rho = obj.state_GV(:,:,:,Parameters.Rho_);
         obj.U   = obj.state_GV(:,:,:,Parameters.U_) ./...
            obj.state_GV(:,:,:,Parameters.Rho_);
         obj.B   = obj.state_GV(:,:,:,Parameters.B_);
         obj.P   = obj.state_GV(:,:,:,Parameters.P_);
         
      end
      
      function obj = get.Energy(obj)
         %GET.ENERGY Calculate the energy in each cell
         
         obj.Energy = obj.P / (Const.gamma-1) + ...
            0.5*obj.Rho.*sum(obj.U.^2,4) + 0.5*sum(obj.B.^2,4);
      end
      
      function obj = update_state(obj,grid,faceFlux,source,time)
         %UPDATE_STATE Update the state variables in one timestep.
         %
         %INPUTS:
         % grid:     class of grid
         % state_VG: state variables
         % faceFlux: class of fluxes
         % source:   class of sources
         % dt:       timestep
         %OUTPUTS:
         % stateNew_VG: updated states
         %
         % Hongyang Zhou, hyzhou@umich.edu
         
         stateNew_GV = obj.state_GV;
         state_GV    = obj.state_GV;
         source_GV   = source.source_GV;
         dt          = time.dt;
         
         CellSize_D = grid.CellSize_D;
         
         iMin = Parameters.iMin;
         iMax = Parameters.iMax;
         jMin = Parameters.jMin;
         jMax = Parameters.jMax;
         kMin = Parameters.kMin;
         kMax = Parameters.kMax;
         Rho_ = Parameters.Rho_;
         P_   = Parameters.P_;
         E_   = Parameters.E_;
         U_   = Parameters.U_;
         B_   = Parameters.B_;
         Bz_  = Parameters.Bz_;
         
         Flux_XV = faceFlux.Flux_XV;
         Flux_YV = faceFlux.Flux_YV;
         Flux_ZV = faceFlux.Flux_ZV;
         
         if strcmp(Parameters.GridType,'Cartesian')
            % No need for volume and face if the grid is uniform Cartesian
            
            if ~Parameters.UseConservative
               stateNew_GV(iMin:iMax,jMin:jMax,kMin:kMax,:) = ...
                  state_GV(iMin:iMax,jMin:jMax,kMin:kMax,:) - dt.*(...
                  (Flux_XV(2:end,:,:,:) - Flux_XV(1:end-1,:,:,:))/...
                  CellSize_D(1) +...
                  (Flux_YV(:,2:end,:,:) - Flux_YV(:,1:end-1,:,:))/...
                  CellSize_D(2) +...
                  (Flux_ZV(:,:,2:end,:) - Flux_ZV(:,:,1:end-1,:))/...
                  CellSize_D(3) -...
                  source_GV);            
            else
               stateNew_GV(iMin:iMax,jMin:jMax,kMin:kMax,Rho_:Bz_) = ...
                  state_GV(iMin:iMax,jMin:jMax,kMin:kMax,Rho_:Bz_) - ...
                  dt.*(...
                  (Flux_XV(2:end  ,:,:,Rho_:Bz_) -...
                   Flux_XV(1:end-1,:,:,Rho_:Bz_)) / CellSize_D(1) + ...
                  (Flux_YV(:,2:end  ,:,Rho_:Bz_) - ...
                   Flux_YV(:,1:end-1,:,Rho_:Bz_)) / CellSize_D(2) + ...
                  (Flux_ZV(:,:,2:end  ,Rho_:Bz_) - ...
                   Flux_ZV(:,:,1:end-1,Rho_:Bz_)) / CellSize_D(3) + ...
                  source_GV(:,:,:,Rho_:Bz_));
               
               gamma = Const.gamma;
               
               state_GV(iMin:iMax,jMin:jMax,kMin:kMax,E_) = ...
                  state_GV(iMin:iMax,jMin:jMax,kMin:kMax,P_) / ...
                  (gamma-1) + ...
                  0.5./state_GV(iMin:iMax,jMin:jMax,kMin:kMax,Rho_).*...
                  sum(state_GV(iMin:iMax,jMin:jMax,kMin:kMax,U_).^2,4) +...
                  0.5*sum(state_GV(iMin:iMax,jMin:jMax,kMin:kMax,B_).^2,4);
               
               stateNew_GV(iMin:iMax,jMin:jMax,kMin:kMax,E_) = ...
                  state_GV(iMin:iMax,jMin:jMax,kMin:kMax,E_) - dt.*(...
                  (Flux_XV(2:end,:,:,E_) - Flux_XV(1:end-1,:,:,E_))/...
                  CellSize_D(1) + ...
                  (Flux_YV(:,2:end,:,E_) - Flux_YV(:,1:end-1,:,E_))/...
                  CellSize_D(2) + ...
                  (Flux_ZV(:,:,2:end,E_) - Flux_ZV(:,:,1:end-1,E_))/...
                  CellSize_D(3) + ...
                  source_GV(:,:,:,E_));
               
               stateNew_GV(iMin:iMax,jMin:jMax,kMin:kMax,P_) = ...
                  (gamma-1)* ...
                  (stateNew_GV(iMin:iMax,jMin:jMax,kMin:kMax,E_) - ...
                  0.5./stateNew_GV(iMin:iMax,jMin:jMax,kMin:kMax,Rho_).*...
                  sum(...
                  stateNew_GV(iMin:iMax,jMin:jMax,kMin:kMax,U_).^2,4) -...
                  0.5*sum(...
                  stateNew_GV(iMin:iMax,jMin:jMax,kMin:kMax,B_).^2,4) );
            end
         else
            % Need volume and face
            stateNew_GV = 0;
            
         end
         
         obj.state_GV = stateNew_GV;
      end
      
      function plot(obj,varargin)
         
         obj = obj.GetState;
         
         iMin = Parameters.iMin;
         iMax = Parameters.iMax;
         jMin = Parameters.jMin;
         jMax = Parameters.jMax;
         kMin = Parameters.kMin;
         kMax = Parameters.kMax;
         
         varname = varargin{1};
         grid    = varargin{2};
         iStep   = varargin{3};
         
         switch varname
            case 'rho'
               var = obj.Rho;
            case 'ux'
               var = obj.U(:,:,:,1);
            case 'p'
               var = obj.P;
            otherwise
               error('unknown plotting varname!')
         end
         
         x = grid.getX;
         
         var = squeeze(var); % Unify dimensions
         var = var(iMin:iMax,jMin:jMax,kMin:kMax); % Remove ghost cells
         
         figure(2)
         plot(x,var,'-*','LineWidth',1.5)
         
         title(sprintf('iStep=%d',iStep))
         legend(varname)
         %axis([1 201 0 2.0])
         
         %M(k) = getframe;
      end
      
      
   end
end
