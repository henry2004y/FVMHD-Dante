classdef FaceValue < handle
   %FaceValue Class of left and right state on the faces.
   % This class holds all data and operations pertinent to the face values.
   
   %======================== MEMBERS =================================
   properties
      LState_XV(:,:,:,:)  %double {mustBeReal}
      RState_XV(:,:,:,:)  %double {mustBeReal}
      LState_YV(:,:,:,:)  %double {mustBeReal}
      RState_YV(:,:,:,:)  %double {mustBeReal}
      LState_ZV(:,:,:,:)  %double {mustBeReal}
      RState_ZV(:,:,:,:)  %double {mustBeReal}
   end
   
   %======================== CONSTRUCTORS ============================
   methods
      function obj = FaceValue
         %FaceValue Construct an instance of this class
         %
      end
   end
   
   %======================== METHODS =================================
   methods (Access = public)
      function calc_face_value(obj,state)
         %CALC_FACE_VALUE Calculate the face values.
         % For uniform grid, there is no need to take grid sizes into
         % account.
         
         state_GV = state.state_GV;
         
         iMin = Parameters.iMin;
         iMax = Parameters.iMax;
         jMin = Parameters.jMin;
         jMax = Parameters.jMax;
         kMin = Parameters.kMin;
         kMax = Parameters.kMax;
         
         if Parameters.Order == 1 % 1st order
            obj.LState_XV = state_GV(iMin-1:iMax,jMin:jMax,kMin:kMax,:);
            obj.RState_XV = state_GV(iMin:iMax+1,jMin:jMax,kMin:kMax,:);
            obj.LState_YV = state_GV(iMin:iMax,jMin-1:jMax,kMin:kMax,:);
            obj.RState_YV = state_GV(iMin:iMax,jMin:jMax+1,kMin:kMax,:);
            obj.LState_ZV = state_GV(iMin:iMax,jMin:jMax,kMin-1:kMax,:);
            obj.RState_ZV = state_GV(iMin:iMax,jMin:jMax,kMin:kMax+1,:);
         else % 2nd order
            % Compute and limit slopes
            
            % Get slope with limiters
            switch Parameters.limiter
               case 'MC'
                  % Find dq_j = minmod{fwd diff, bwd diff, cntrl diff}
                  dqR_X = ...
                     state_GV(iMin  :iMax+2,jMin:jMax,kMin:kMax,:) -...
                     state_GV(iMin-1:iMax+1,jMin:jMax,kMin:kMax,:);
                  dqL_X = ...
                     state_GV(iMin-1:iMax+1,jMin:jMax,kMin:kMax,:) - ...
                     state_GV(iMin-2:iMax  ,jMin:jMax,kMin:kMax,:);
                  dqC_X = ...
                     state_GV(iMin  :iMax+2,jMin:jMax,kMin:kMax,:) - ...
                     state_GV(iMin-2:iMax  ,jMin:jMax,kMin:kMax,:);
                  
                  dq_X = obj.minmod(dqR_X,dqL_X,dqC_X);
                  
                  dqR_Y = ...
                     state_GV(iMin:iMax,jMin  :jMax+2,kMin:kMax,:) -...
                     state_GV(iMin:iMax,jMin-1:jMax+1,kMin:kMax,:);
                  dqL_Y = ...
                     state_GV(iMin:iMax,jMin-1:jMax+1,kMin:kMax,:) - ...
                     state_GV(iMin:iMax,jMin-2:jMax  ,kMin:kMax,:);
                  
                  dq_Y = obj.minmod(dqR_Y,dqL_Y,dqC_Y);
                  
                  dqR_Z = ...
                     state_GV(iMin:iMax,jMin:jMax,kMin  :kMax+2,:) - ...
                     state_GV(iMin:iMax,jMin:jMax,kMin-1:kMax+1,:);
                  dqL_Z = ...
                     state_GV(iMin:iMax,jMin:jMax,kMin-1:kMax+1,:) - ...
                     state_GV(iMin:iMax,jMin:jMax,kMin-2:kMax  ,:);
                  
                  dq_Z = obj.minmod(dqR_Z,dqL_Z,dqC_Z);
               case 'MM' % Minmod limiter
                  % Find dq_j = minmod{fwd diff, bwd diff}
                  
                  dqR_X = ...
                     state_GV(iMin  :iMax+2,jMin:jMax,kMin:kMax,:) -...
                     state_GV(iMin-1:iMax+1,jMin:jMax,kMin:kMax,:);
                  dqL_X = ...
                     state_GV(iMin-1:iMax+1,jMin:jMax,kMin:kMax,:) - ...
                     state_GV(iMin-2:iMax  ,jMin:jMax,kMin:kMax,:);
                  
                  dq_X = obj.minmod(dqR_X,dqL_X);
                  
                  dqR_Y = ...
                     state_GV(iMin:iMax,jMin  :jMax+2,kMin:kMax,:) -...
                     state_GV(iMin:iMax,jMin-1:jMax+1,kMin:kMax,:);
                  dqL_Y = ...
                     state_GV(iMin:iMax,jMin-1:jMax+1,kMin:kMax,:) - ...
                     state_GV(iMin:iMax,jMin-2:jMax  ,kMin:kMax,:);
                  
                  dq_Y = obj.minmod(dqR_Y,dqL_Y);
                  
                  dqR_Z = ...
                     state_GV(iMin:iMax,jMin:jMax,kMin  :kMax+2,:) - ...
                     state_GV(iMin:iMax,jMin:jMax,kMin-1:kMax+1,:);
                  dqL_Z = ...
                     state_GV(iMin:iMax,jMin:jMax,kMin-1:kMax+1,:) - ...
                     state_GV(iMin:iMax,jMin:jMax,kMin-2:kMax  ,:);
                  
                  dq_Z = obj.minmod(dqR_Z,dqL_Z);
            end
            
            % Linear interpolation onto edge centers
            obj.LState_XV = state_GV(iMin-1:iMax,jMin:jMax,kMin:kMax,:)+...
               0.5*dq_X(1:end-1,:,:,:);
            obj.RState_XV = state_GV(iMin:iMax+1,jMin:jMax,kMin:kMax,:)-...
               0.5*dq_X(2:end  ,:,:,:);
            obj.LState_YV = state_GV(iMin:iMax,jMin-1:jMax,kMin:kMax,:)+...
               0.5*dq_Y(:,1:end-1,:,:);
            obj.RState_YV = state_GV(iMin:iMax,jMin:jMax+1,kMin:kMax,:)-...
               0.5*dq_Y(:,2:end  ,:,:);
            obj.LState_ZV = state_GV(iMin:iMax,jMin:jMax,kMin-1:kMax,:)+...
               0.5*dq_Z(:,:,1:end-1,:);
            obj.RState_ZV = state_GV(iMin:iMax,jMin:jMax,kMin:kMax+1,:)-...
               0.5*dq_Z(:,:,2:end  ,:);
            
         end
      end
      
   end
   
   methods (Static)
      function m = minmod(a,b,c)
         %MINMOD
         % For three inputs, use Harten's generalized definition.
         %OUTPUT:
         % m: zero if opposite sign, otherwise the one of smaller magnitude.
         
         if nargin == 2  % Two input arguments
            m = (sign(a) + sign(b))/2.*min(abs(a),abs(b));
         elseif nargin == 3  % Three input arguments
            s = (sign(a) + sign(b) + sign(c))/3;
            if abs(s)==1
               m = s*min(abs(a),abs(b),abs(c));
            else
               m = zeros(size(a));
            end
         end
      end
      
   end
   
end
