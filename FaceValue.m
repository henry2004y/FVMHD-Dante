classdef FaceValue < handle
   %FaceValue Class of left and right state on the faces.
   % This class holds all data and operations pertinent to the face values.
   
   %======================== MEMBERS =================================
   properties
      LState_VX(:,:,:,:)  double {mustBeReal}
      RState_VX(:,:,:,:)  double {mustBeReal}
      LState_VY(:,:,:,:)  double {mustBeReal}
      RState_VY(:,:,:,:)  double {mustBeReal}
      LState_VZ(:,:,:,:)  double {mustBeReal}
      RState_VZ(:,:,:,:)  double {mustBeReal}
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
      function obj = calc_face_value(obj,state_VG)
         
         iMin = Parameters.iMin;
         iMax = Parameters.iMax;
         jMin = Parameters.jMin;
         jMax = Parameters.jMax;
         kMin = Parameters.kMin;
         kMax = Parameters.kMax;
                
         if Parameters.Order == 1 % 1st order
            obj.LState_VX = state_VG(:,iMin-1:iMax,jMin:jMax,kMin:kMax);
            obj.RState_VX = state_VG(:,iMin:iMax+1,jMin:jMax,kMin:kMax);
            obj.LState_VY = state_VG(:,iMin:iMax,jMin-1:jMax,kMin:kMax);
            obj.RState_VY = state_VG(:,iMin:iMax,jMin:jMax+1,kMin:kMax);
            obj.LState_VZ = state_VG(:,iMin:iMax,jMin:jMax,kMin-1:kMax);
            obj.RState_VZ = state_VG(:,iMin:iMax,jMin:jMax,kMin:kMax+1);
         else % 2nd order
            % Compute and limit slopes
            %dq_X=zeros(Parameters.nVar,Parameters.nI+1);

            switch Parameters.limiter
               case 'MC'
                  
               case 'MM' % Minmod limiter
                  % Find dq_j = minmod{fwd diff, bwd diff}                  
                  dqR_X = squeeze(...
                     state_VG(:,iMin+1:iMax+2,jMin:jMax,kMin:kMax) -...
                     state_VG(:,iMin:iMax+1,jMin:jMax,kMin:kMax));
                  dqL_X = squeeze(...
                     state_VG(:,iMin:iMax+1,jMin:jMax,kMin:kMax) - ...
                     state_VG(:,iMin-1:iMax,jMin:jMax,kMin:kMax));
                  
                  dq_X = minmod(dqR_X,dqL_X);
            end
            
            % Get the gradient of states
            
            % Linear interpolation onto edge centers (with limiters)
            
            
         end       
      end
     
   end
   
end

