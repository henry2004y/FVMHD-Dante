classdef Boundary < handle
   %BOUNDARY Class of boundaries.
   % This class hold all data and operations pertinent to boundary
   % settings.
   
   %======================== MEMBERS =================================
   properties
      TypeBC char
   end
   
   %======================== CONSTRUCTORS ============================
   methods
      function obj = Boundary()
         %BOUNDARY Construct an instance of this class
         obj.TypeBC = Parameters.BC;
      end
   end
   
   %======================== METHODS =================================
   methods (Access = public)
      function state_GV = set_cell_boundary(obj,state_GV)
         %SET_CELL_BOUNDARY applies boundary conditions to the chosen grid.
         
         
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
         
         switch obj.TypeBC
            case 'periodic'
               state_GV(1:nG,jMin:jMax,kMin:kMax,:) = ...
                  state_GV(iMax-nG+1:iMax,jMin:jMax,kMin:kMax,:);
               state_GV(iMin:iMax,1:nG,kMin:kMax,:) = ...
                  state_GV(iMin:iMax,jMax-nG+1:jMax,kMin:kMax,:);
               state_GV(iMin:iMax,jMin:jMax,1:nG,:) = ...
                  state_GV(iMin:iMax,jMin:jMax,kMax-nG+1:kMax,:);
               state_GV(iMax+1:iMaxAll,jMin:jMax,kMin:kMax,:) = ...
                  state_GV(iMin:iMin+nG-1,jMin:jMax,kMin:kMax,:);
               state_GV(iMin:iMax,jMax+1:jMaxAll,kMin:kMax,:) = ...
                  state_GV(iMin:iMax,jMin:jMin+nG-1,kMin:kMax,:);
               state_GV(iMin:iMax,jMin:jMax,kMax+1:kMaxAll,:) = ...
                  state_GV(iMin:iMax,jMin:jMax,kMin:kMin+nG-1,:);               
            case 'fixed'
               % Note: this needs to be improved! No magic numbers inside source
               % codes!
               CellState_VG = [2 0 0 0 0 0 0 1]';
               state_GV(:,1:nG,jMin:jMax,kMin:kMax) = CellState_VG;
               CellState_VG = [1 0 0 0 0 0 0 1]';
               state_GV(:,iMax+1:iMaxAll,jMin:jMax,kMin:kMax) = CellState_VG;
               state_GV(:,iMin:iMax,1:nG,kMin:kMax) = ...
                  state_GV(:,iMin:iMax,jMin-nG+1:jMax,kMin:kMax);
               state_GV(:,iMin:iMax,jMin:jMax,1:nG) = ...
                  state_GV(:,iMin:iMax,jMin:jMax,kMin-nG+1:kMax);
               state_GV(:,iMin:iMax,jMax+1:jMaxAll,kMin:kMax) = ...
                  state_GV(:,iMin:iMax,jMin:jMin+nG-1,kMin:kMax);
               state_GV(:,iMin:iMax,jMin:jMax,kMax+1:kMaxAll) = ...
                  state_GV(:,iMin:iMax,jMin:jMax,kMin:kMin+nG-1);
            case 'float'
               % There might be some minor problem in implementing high order
               % float boundary condition.
               state_GV(:,1:nG,jMin:jMax,kMin:kMax) = ...
                  repmat(state_GV(:,iMin,jMin:jMax,kMin:kMax),[1 nG 1 1]);
               state_GV(:,iMin:iMax,1:nG,kMin:kMax) = ...
                  repmat(state_GV(:,iMin:iMax,jMin,kMin:kMax),[1 1 nG 1]);
               state_GV(:,iMin:iMax,jMin:jMax,1:nG) = ...
                  repmat(state_GV(:,iMin:iMax,jMin:jMax,kMin),[1 1 1 nG]);
               state_GV(:,iMax+1:iMaxAll,jMin:jMax,kMin:kMax) = ...
                  repmat(state_GV(:,iMax,jMin:jMax,kMin:kMax),[1 nG 1 1]);
               state_GV(:,iMin:iMax,jMax+1:jMaxAll,kMin:kMax) = ...
                  repmat(state_GV(:,iMin:iMax,jMax,kMin:kMax),[1 1 nG 1]);
               state_GV(:,iMin:iMax,jMin:jMax,kMax+1:kMaxAll) = ...
                  repmat(state_GV(:,iMin:iMax,jMin:jMax,kMax),[1 1 1 nG]);
               
               
            otherwise
               error('unknown boundary type!')
         end
         
      end
   end
end

