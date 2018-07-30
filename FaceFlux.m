classdef FaceFlux < handle
   %FaceClass The face flux variable class
   %   Detailed explanation goes here
   
   properties
      Flux_VX(:,:,:,:)  double {mustBeReal}
      Flux_VY(:,:,:,:)  double {mustBeReal}
      Flux_VZ(:,:,:,:)  double {mustBeReal}
   end
   
   methods
      function obj = FaceFlux
         %FaceClass Construct an instance of this class
         %   Detailed explanation goes here
         
         nVar = Parameters.nVar;
         GridSize = Parameters.GridSize;
         
%          obj.Flux_VX = Inf([nVar,GridSize+1]);
%          obj.Flux_VY = Inf([nVar,GridSize+1]);
%          obj.Flux_VZ = Inf([nVar,GridSize+1]);
      end
      
      function obj = calc_face_flux(obj,faceValue)
                  
         if Parameters.Order == 1
         
            [obj.Flux_VX, obj.Flux_VY, obj.Flux_VZ] = ...
               faceValue.get_physical_flux;
            
            %obj.Flux_VX(8,1:5)
            
            [obj.Flux_VX, obj.Flux_VY, obj.Flux_VZ] = ...
               faceValue.add_numerical_flux(...
               obj.Flux_VX,obj.Flux_VY,obj.Flux_VZ);
         
            %obj.Flux_VX(8,1:5)
         
         end
      end
      
   end
   
end

