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
      function calc_timestep(obj,grid,faceFlux)
         %CALC_TIMESTEP Calculate the timestep.
         %  Return the timestep under CFL condition.
         
         CFL = Parameters.CFL;
         CellSize_D = grid.CellSize_D;
         
         Cmax_XF = faceFlux.Cmax_XF;
         Cmax_YF = faceFlux.Cmax_YF;
         Cmax_ZF = faceFlux.Cmax_ZF;
         
         nI = Parameters.nI; nJ = Parameters.nJ; nK = Parameters.nK;
         if Parameters.TimeAccurate      
            Cmax_XG = max(Cmax_XF(2:end,:,:),Cmax_XF(1:end-1,:,:));
            Cmax_YG = max(Cmax_YF(:,2:end,:),Cmax_YF(:,1:end-1,:));
            Cmax_ZG = max(Cmax_ZF(:,:,2:end),Cmax_ZF(:,:,1:end-1));
            
            time_G = CFL ./ (...
               (nI > 1)*Cmax_XG/CellSize_D(1) + ...
               (nJ > 1)*Cmax_YG/CellSize_D(2) + ...
               (nK > 1)*Cmax_ZG/CellSize_D(3));
            
            obj.dt = min(time_G(:));
         else
            obj.dt = CFL ./ ( ...
               (nI > 1)*Cmax_XF/CellSize_D(1) + ...
               (nJ > 1)*Cmax_YF/CellSize_D(2) + ...
               (nK > 1)*Cmax_ZF/CellSize_D(3));
         end
      end
   end
   
end