classdef Source < handle
   %SOURCE Class of source terms.
   % This class holds all data and operations pertinent to the sources.
   
   %======================== MEMBERS =================================
   properties
      source_VG(:,:,:,:) double {mustBeReal}
   end
   
   %======================== METHODS =================================
   methods
      function calc_source(obj,grid,state_VG)
         %CALC_SOURCE Calculate all the source terms.
         %
         
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
         DivB = divergence_ndgrid(grid.X,grid.Y,grid.Z,state_VG(B_,:,:,:));
         % Add a singleton dimension for the sake of matrix operation later
         DivB = reshape(DivB,[1 size(DivB)]);
         
         DivU = divergence_ndgrid(grid.X,grid.Y,grid.Z,...
            state_VG(U_,:,:,:)./state_VG(Rho_,:,:,:));
         % Add a singleton dimension for the sake of matrix operation later
         DivU = reshape(DivU,[1 size(DivU)]);
         
         obj.source_VG = Inf([nVar GridSize]);
         
         if Parameters.Order == 1
            obj.source_VG(Rho_,:,:,:) = 0;
            obj.source_VG(U_,:,:,:) = ...
               -state_VG(B_,iMin:iMax,jMin:jMax,kMin:kMax).*...
               DivB(1,iMin:iMax,jMin:jMax,kMin:kMax);
            
            obj.source_VG(B_,:,:,:) = ...
               state_VG(U_,iMin:iMax,jMin:jMax,kMin:kMax).*...
               DivB(1,iMin:iMax,jMin:jMax,kMin:kMax);
            
            if ~Parameters.UseConservative
               obj.source_VG(P_,:,:,:) = -(gamma-1) * ...
                  state_VG(P_,iMin:iMax,jMin:jMax,kMin:kMax).*...
                  DivU(1,iMin:iMax,jMin:jMax,kMin:kMax);
            else
               obj.source_VG(E_,:,:,:) = ...
                  -sum(state_VG(U_,iMin:iMax,jMin:jMax,kMin:kMax) .*...
                  state_VG(B_,iMin:iMax,jMin:jMax,kMin:kMax),1) .* ...
                  DivB(1,iMin:iMax,jMin:jMax,kMin:kMax);
            end
            
         elseif Parameters.Order == 2
            
         else
            
         end
      end
   end
end

