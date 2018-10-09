classdef StateClass < handle
   %State The class of physical states in MHD
   %   Primitive variables include:
   % rho: density
   % u:   velocity (3D)
   % B:   magnetic field (3D)
   % p:   pressure
   % e:   energy
   
   properties
      Rho(:,:,:)  double {mustBeReal, mustBeGreaterThan(Rho, 0)}
      U(3,:,:,:)  double {mustBeReal}
      B(3,:,:,:)  double {mustBeReal}
      P(:,:,:)    double {mustBeReal, mustBeGreaterThan(P, 0)}
   end
   
   properties (Dependent)
      % Energy
      Energy double {mustBeReal, mustBeGreaterThan(Energy, 0)}
   end
      
   methods
      function obj = StateClass(density,velocity,Bfield,pressure)
         %State Construct an instance of this class
         %   Detailed explanation goes here
         
         disp('Creating state variables...')
         
         if nargin == 0
            disp('Using default grid')
            obj.Rho = ones(GridGen.GridSize);
            % May be I don't need initialization this way.
         end
         % Call asset constructor
         obj.Rho = density;
         obj.U   = velocity;
         obj.B   = Bfield;
         obj.P   = pressure;
         
      end
      
      function state_VG = SetState(obj)
         % Reorganize data structure for efficiency
         state_VG = Inf([Parameters.nVar,Parameters.FullSize]);
         
         state_VG(Parameters.Rho_,:,:,:) = obj.Rho;
         state_VG(Parameters.U_,:,:,:) = obj.U .*...
            state_VG(Parameters.Rho_,:,:,:);
         state_VG(Parameters.B_,:,:,:) = obj.B;
         state_VG(Parameters.P_,:,:,:) = obj.P;
         
      end
  
      function obj = GetState(obj,state_VG)
         % Reorganize data structure for post-processing
         
         obj.Rho = state_VG(Parameters.Rho_,:,:,:);
         obj.U   = state_VG(Parameters.U_,:,:,:) ./...
            state_VG(Parameters.Rho_,:,:,:);
         obj.B   = state_VG(Parameters.B_,:,:,:);
         obj.P   = state_VG(Parameters.P_,:,:,:);
         
      end      
           
      function obj = get.Energy(obj)
         %GET.ENERGY Calculate the energy in each cell

         obj.Energy = obj.P / (Const.gamma-1) + ...
            0.5*obj.Rho.*sum(obj.U.^2,1) + 0.5*sum(obj.B.^2,1);
      end
            
      function plot(obj,varargin)
         
         iMin = Parameters.iMin;
         iMax = Parameters.iMax;
         jMin = Parameters.jMin;
         jMax = Parameters.jMax;
         kMin = Parameters.kMin;
         kMax = Parameters.kMax;         
         
         varname = varargin{1};
         iStep = varargin{2};
         
         switch varname
            case 'rho'
               y = obj.Rho;
            case 'ux'
               y = obj.U(1,:,:,:);
            otherwise
               error('unknown plotting varname!')
         end
         
         y = squeeze(y);
         y = y(iMin:iMax,jMin:jMax,kMin:kMax);
         
         %figure;
         plot(y,'LineWidth',1.5)
         
         title(sprintf('iStep=%d',iStep))
         axis([1 101 0 2.0])
         
         %M(k) = getframe;
      end
      
      
   end
end

