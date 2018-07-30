function [source_VG] = calc_source(grid,state_VG)
%calc_source Summary of this function goes here
%   Detailed explanation goes here

gamma = Const.gamma;

nVar = Parameters.nVar;
GridSize = Parameters.GridSize;
%FullSize = Parameters.FullSize;
Rho_ = Parameters.Rho_;
Ux_  = Parameters.Ux_;
Uy_  = Parameters.Uy_;
Uz_  = Parameters.Uz_;
Bx_  = Parameters.Bx_;
By_  = Parameters.By_;
Bz_  = Parameters.Bz_;
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
% ndgrid to meshgrid
% X = permute(grid.X,[2 1 3]);
% Y = permute(grid.Y,[2 1 3]);
% Z = permute(grid.Z,[2 1 3]);
% This may not work well and should be check or rewritten!!!
% The problem for this is that it requires meshgrid format instead of 
% ndgrid
% DivB = divergence(X,Y,Z,...
%    squeeze(state_VG(Bx_,:,:,:)),...
%    squeeze(state_VG(By_,:,:,:)),...
%    squeeze(state_VG(Bz_,:,:,:)));
DivB = divergence_ndgrid(grid.X,grid.Y,grid.Z,state_VG(B_,:,:,:));

% Add a singleton dimension for the sake of matrix operation later
DivB = reshape(DivB,[1 size(DivB)]);

% DivU = divergence(X,Y,Z,...
%    squeeze(state_VG(Ux_,:,:,:)),...
%    squeeze(state_VG(Uy_,:,:,:)),...
%    squeeze(state_VG(Uz_,:,:,:)));
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

%    source_VG(Ux_,:,:,:) = ...
%       -squeeze(state_VG(Bx_,iMin:iMax,jMin:jMax,kMin:kMax)).*...
%       DivB(iMin:iMax,jMin:jMax,kMin:kMax);
%    source_VG(Uy_,:,:,:) = ...
%       -squeeze(state_VG(By_,iMin:iMax,jMin:jMax,kMin:kMax)).*...
%       DivB(iMin:iMax,jMin:jMax,kMin:kMax);
%    source_VG(Uz_,:,:,:) = ...
%       -squeeze(state_VG(Bz_,iMin:iMax,jMin:jMax,kMin:kMax)).*...
%       DivB(iMin:iMax,jMin:jMax,kMin:kMax);   
   
   
   source_VG(B_,:,:,:) = ...
      state_VG(U_,iMin:iMax,jMin:jMax,kMin:kMax).*...
      DivB(1,iMin:iMax,jMin:jMax,kMin:kMax);

%    source_VG(Bx_,:,:,:) = ...
%       -squeeze(state_VG(Ux_,iMin:iMax,jMin:jMax,kMin:kMax)).*...
%       DivB(iMin:iMax,jMin:jMax,kMin:kMax);
%    source_VG(By_,:,:,:) = ...
%       -squeeze(state_VG(Uy_,iMin:iMax,jMin:jMax,kMin:kMax)).*...
%       DivB(iMin:iMax,jMin:jMax,kMin:kMax);
%    source_VG(Bz_,:,:,:) = ...
%       -squeeze(state_VG(Uz_,iMin:iMax,jMin:jMax,kMin:kMax)).*...
%       DivB(iMin:iMax,jMin:jMax,kMin:kMax); 


   if ~Parameters.UseConservative
      source_VG(P_,:,:,:) = ...
         -(gamma-1) * state_VG(P_,iMin:iMax,jMin:jMax,kMin:kMax).*...
         DivU(1,iMin:iMax,jMin:jMax,kMin:kMax);
   
      % I searched the bug in my code for a long time, and it finally turns
      % out that the source term in pressure due to divU is the key here. 
      % This term represents work done by compressible flow. If I set this
      % to 0, at least I can run many shock and discontinuity tests. This 
      % also explains why I don't have issue with stationary shock --- 
      % because velocity is always zero!
      %source_VG(P_,:,:,:) = 0.0;
   else
      source_VG(E_,:,:,:) = -sum(state_VG(U_,iMin:iMax,jMin:jMax,kMin:kMax)...
         .*state_VG(B_,iMin:iMax,jMin:jMax,kMin:kMax),1) .* ...
         DivB(1,iMin:iMax,jMin:jMax,kMin:kMax);
   end
   
elseif Parameters.Order == 2
   
else
   
end

end

