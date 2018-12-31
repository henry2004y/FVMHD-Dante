classdef Time < handle
   %TIME Class of timestep control.
   % This class contains data and operations pertinent to timesteps.
   
   %======================== MEMBERS =================================
   properties
      dt
   end
   
   %======================== CONSTRUCTORS ============================
   methods
      function obj = Time()
         %TIME Construct an instance of this class
         %
         obj.dt = 0;
      end
   end
   
   %======================== METHODS =================================
   methods (Access = public)
      function calc_timestep(obj,grid,state_GV)
         %CALC_TIMESTEP Calculate the timestep.
         %  Return the timestep under CFL condition.
         
         CFL = Parameters.CFL;
         CellSize_D = grid.CellSize_D;
         
         [Cmax_GX,Cmax_GY,Cmax_GZ] = obj.get_speed_max(state_GV);
         
         % [Cmax_XF,Cmax_YF,Cmax_ZF] = faceValue.get_speed_max;
         %
         % Cmax_XG = max(Cmax_XF(1,1:end-1,:,:),Cmax_XF(1,2:end,:,:));
         % Cmax_YG = max(Cmax_YF(1,:,1:end-1,:),Cmax_YF(1,:,2:end,:));
         % Cmax_ZG = max(Cmax_ZF(1,:,:,1:end-1),Cmax_ZF(1,:,:,2:end));
         
         if Parameters.TimeAccurate
            timestep = CFL ./ ...
               (Cmax_GX/CellSize_D(1) + ...
               Cmax_GY/CellSize_D(2) + ...
               Cmax_GZ/CellSize_D(3));
            obj.dt = min(timestep(:));
         else
            obj.dt = CFL / ...
               (Cmax_GX/CellSize_D(1) + ...
               Cmax_GY/CellSize_D(2) + ...
               Cmax_GZ/CellSize_D(3));
         end
      end
   end
   
   methods (Static)
      function [Cmax_GX,Cmax_GY,Cmax_GZ] = get_speed_max(state_GV)
         %GET_SPEED_MAX Get the cell-centered maximum speed.
          
         gamma = Const.gamma;
         Rho_ = Parameters.Rho_;
         Ux_  = Parameters.Ux_;
         Uy_  = Parameters.Uy_;
         Uz_  = Parameters.Uz_;
         Bx_  = Parameters.Bx_;
         By_  = Parameters.By_;
         Bz_  = Parameters.Bz_;
         P_   = Parameters.P_;
         U_   = Parameters.U_;
         B_   = Parameters.B_;
         
         Cs2_G = gamma* state_GV(:,:,:,P_) ./ state_GV(:,:,:,Rho_);
         Ca2_G = (state_GV(:,:,:,Bx_).^2 + state_GV(:,:,:,By_).^2 + ...
            state_GV(:,:,:,Bz_).^2) ./ state_GV(:,:,:,Rho_);
         Can2_XG = (state_GV(:,:,:,Bx_).^2 ) ./ state_GV(:,:,:,Rho_);
         Can2_YG = (state_GV(:,:,:,By_).^2 ) ./ state_GV(:,:,:,Rho_);
         Can2_ZG = (state_GV(:,:,:,Bz_).^2 ) ./ state_GV(:,:,:,Rho_);
         
         Cmax_GX = abs(state_GV(:,:,:,Ux_)./state_GV(:,:,:,Rho_)) + ...
            sqrt(0.5*(Cs2_G + Ca2_G + sqrt((Cs2_G + Ca2_G).^2 - ...
            4*Cs2_G.*Can2_XG)));
         Cmax_GY = abs(state_GV(:,:,:,Uy_)./state_GV(:,:,:,Rho_)) + ...
            sqrt(0.5*(Cs2_G + Ca2_G + sqrt((Cs2_G + Ca2_G).^2 - ...
            4*Cs2_G.*Can2_YG)));
         Cmax_GZ = abs(state_GV(:,:,:,Uz_)./state_GV(:,:,:,Rho_)) + ...
            sqrt(0.5*(Cs2_G + Ca2_G + sqrt((Cs2_G + Ca2_G).^2 - ...
            4*Cs2_G.*Can2_ZG)));
      end
      
   end
end

