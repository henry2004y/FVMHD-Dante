function [source_VG] = calc_source(grid,state_VG)
%calc_source Summary of this function goes here
%   Detailed explanation goes here

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

source_VG = Inf([nVar GridSize]);

if Parameters.Order == 1
   source_VG(Rho_,:,:,:) = 0;
    source_VG(U_,:,:,:) = ...
       -state_VG(B_,iMin:iMax,jMin:jMax,kMin:kMax).*...
       DivB(1,iMin:iMax,jMin:jMax,kMin:kMax);
   
   source_VG(B_,:,:,:) = ...
      state_VG(U_,iMin:iMax,jMin:jMax,kMin:kMax).*...
      DivB(1,iMin:iMax,jMin:jMax,kMin:kMax);
 
   if ~Parameters.UseConservative
      source_VG(P_,:,:,:) = ...
         -(gamma-1) * state_VG(P_,iMin:iMax,jMin:jMax,kMin:kMax).*...
         DivU(1,iMin:iMax,jMin:jMax,kMin:kMax);
   else
      source_VG(E_,:,:,:) = -sum(state_VG(U_,iMin:iMax,jMin:jMax,kMin:kMax)...
         .*state_VG(B_,iMin:iMax,jMin:jMax,kMin:kMax),1) .* ...
         DivB(1,iMin:iMax,jMin:jMax,kMin:kMax);
   end
   
elseif Parameters.Order == 2
   
else
   
end

end

