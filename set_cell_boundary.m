function [state_VG] = set_cell_boundary(state_VG)
%set_cell_boundary Summary of this function goes here
%   Detailed explanation goes here

nG = Parameters.nG;
iMin = Parameters.iMin;
iMax = Parameters.iMax;
jMin = Parameters.jMin;
jMax = Parameters.jMax;
kMin = Parameters.kMin;
kMax = Parameters.kMax;

iMaxAll = Parameters.iMaxAll;
jMaxAll = Parameters.jMaxAll;
kMaxAll = Parameters.kMaxAll;

switch Parameters.BC
   case 'periodic'
      state_VG(:,1:nG,jMin:jMax,kMin:kMax) = ...
         state_VG(:,iMax-nG+1:iMax,jMin:jMax,kMin:kMax);
      state_VG(:,iMin:iMax,1:nG,kMin:kMax) = ...
         state_VG(:,iMin:iMax,jMax-nG+1:jMax,kMin:kMax);
      state_VG(:,iMin:iMax,jMin:jMax,1:nG) = ...
         state_VG(:,iMin:iMax,jMin:jMax,kMax-nG+1:kMax);
      state_VG(:,iMax+1:iMaxAll,jMin:jMax,kMin:kMax) = ...
         state_VG(:,iMin:iMin+nG-1,jMin:jMax,kMin:kMax);  
      state_VG(:,iMin:iMax,jMax+1:jMaxAll,kMin:kMax) = ...
         state_VG(:,iMin:iMax,jMin:jMin+nG-1,kMin:kMax);
      state_VG(:,iMin:iMax,jMin:jMax,kMax+1:kMaxAll) = ...
         state_VG(:,iMin:iMax,jMin:jMax,kMin:kMin+nG-1);
   case 'fixed'
      CellState_VG = [2 0 0 0 0 0 0 1]';
      state_VG(:,1:nG,jMin:jMax,kMin:kMax) = CellState_VG;
      CellState_VG = [1 0 0 0 0 0 0 1]';
      state_VG(:,iMax+1:iMaxAll,jMin:jMax,kMin:kMax) = CellState_VG;
      state_VG(:,iMin:iMax,1:nG,kMin:kMax) = ...
         state_VG(:,iMin:iMax,jMin-nG+1:jMax,kMin:kMax);
      state_VG(:,iMin:iMax,jMin:jMax,1:nG) = ...
         state_VG(:,iMin:iMax,jMin:jMax,kMin-nG+1:kMax);
      state_VG(:,iMin:iMax,jMax+1:jMaxAll,kMin:kMax) = ...
         state_VG(:,iMin:iMax,jMin:jMin+nG-1,kMin:kMax);
      state_VG(:,iMin:iMax,jMin:jMax,kMax+1:kMaxAll) = ...
         state_VG(:,iMin:iMax,jMin:jMax,kMin:kMin+nG-1);
   case 'float'
      state_VG(:,1:nG,jMin:jMax,kMin:kMax) = ...
         repmat(state_VG(:,iMin,jMin:jMax,kMin:kMax),[1 nG 1 1]);
      state_VG(:,iMin:iMax,1:nG,kMin:kMax) = ...
         repmat(state_VG(:,iMin:iMax,jMin,kMin:kMax),[1 1 nG 1]);
      state_VG(:,iMin:iMax,jMin:jMax,1:nG) = ...
         repmat(state_VG(:,iMin:iMax,jMin:jMax,kMin),[1 1 1 nG]);
      state_VG(:,iMax+1:iMaxAll,jMin:jMax,kMin:kMax) = ...
         repmat(state_VG(:,iMax,jMin:jMax,kMin:kMax),[1 nG 1 1]);  
      state_VG(:,iMin:iMax,jMax+1:jMaxAll,kMin:kMax) = ...
         repmat(state_VG(:,iMin:iMax,jMax,kMin:kMax),[1 1 nG 1]);
      state_VG(:,iMin:iMax,jMin:jMax,kMax+1:kMaxAll) = ...
         repmat(state_VG(:,iMin:iMax,jMin:jMax,kMax),[1 1 1 nG]);      
   otherwise
      error('unknown boundary type!')
end
      
      
end

