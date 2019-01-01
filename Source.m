classdef Source < handle
   %SOURCE Class of source terms.
   % This class holds all data and operations pertinent to the sources.
   
   %======================== MEMBERS =================================
   properties
      source_GV(:,:,:,:) double
   end
   
   %======================== METHODS =================================
   methods
      function calc_source(obj,grid,state)
         %CALC_SOURCE Calculate all the source terms.
         %
         % There seems to be no change for 2nd order schemes?
         
         state_GV = state.state_GV;
         
         gamma = Const.gamma;
         
         nVar = Parameters.nVar;
         GridSize = Parameters.GridSize;
         Rho_ = Parameters.Rho_;
         P_   = Parameters.P_;
         E_   = Parameters.E_;
         U_   = Parameters.U_;
         B_   = Parameters.B_;
         
         iMin = Parameters.iMin;
         iMax = Parameters.iMax;
         jMin = Parameters.jMin;
         jMax = Parameters.jMax;
         kMin = Parameters.kMin;
         kMax = Parameters.kMax;
         
         % Calculate divergence of B using central difference
         DivB = divergence_ndgrid(grid.X,grid.Y,grid.Z,state_GV(:,:,:,B_));
         
         DivU = divergence_ndgrid(grid.X,grid.Y,grid.Z,...
            state_GV(:,:,:,U_)./state_GV(:,:,:,Rho_));
         
         obj.source_GV = Inf([GridSize nVar]);
         
         obj.source_GV(:,:,:,Rho_) = 0;
         obj.source_GV(:,:,:,U_) = ...
            -state_GV(iMin:iMax,jMin:jMax,kMin:kMax,B_).*...
            DivB(iMin:iMax,jMin:jMax,kMin:kMax);
         
         obj.source_GV(:,:,:,B_) = ...
            state_GV(iMin:iMax,jMin:jMax,kMin:kMax,U_).*...
            DivB(iMin:iMax,jMin:jMax,kMin:kMax);
         
         if ~Parameters.UseConservative
            obj.source_GV(:,:,:,P_) = -(gamma-1) * ...
               state_GV(iMin:iMax,jMin:jMax,kMin:kMax,P_).*...
               DivU(iMin:iMax,jMin:jMax,kMin:kMax);
         else
            obj.source_GV(:,:,:,E_) = ...
               -sum(state_GV(iMin:iMax,jMin:jMax,kMin:kMax,U_) .*...
               state_GV(iMin:iMax,jMin:jMax,kMin:kMax,B_),4) .* ...
               DivB(iMin:iMax,jMin:jMax,kMin:kMax);
         end
      end
      
   end
end

