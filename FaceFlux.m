classdef FaceFlux < handle
   %FaceFlux Class of fluxes
   %
   
   properties
      Flux_VX(:,:,:,:)  double {mustBeReal}
      Flux_VY(:,:,:,:)  double {mustBeReal}
      Flux_VZ(:,:,:,:)  double {mustBeReal}
   end
   
   methods
      function obj = FaceFlux
         %FACEFLUX Construct an instance of this class
         %
         
         nVar = Parameters.nVar;
         GridSize = Parameters.GridSize;
         
%          obj.Flux_VX = Inf([nVar,GridSize+1]);
%          obj.Flux_VY = Inf([nVar,GridSize+1]);
%          obj.Flux_VZ = Inf([nVar,GridSize+1]);
      end
      
      function obj = calc_face_flux(obj,faceValue)
         % Calculate the face fluxes in generalized coordinates from face 
         % values on the left and right of the face     
         
         [obj.Flux_VX, obj.Flux_VY, obj.Flux_VZ] = ...
            faceValue.get_physical_flux;
         
         [obj.Flux_VX, obj.Flux_VY, obj.Flux_VZ] = ...
            faceValue.add_numerical_flux(...
            obj.Flux_VX,obj.Flux_VY,obj.Flux_VZ);
      end
      
   end
   
end

