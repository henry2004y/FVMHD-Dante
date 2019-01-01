classdef Grid < handle
   %grid The class of simulation grid
   %  This class holds grid information.
   
   %======================== MEMBERS =================================
   properties (SetAccess = private)
      CoordMinMax_D(3,2) double {mustBeReal}
      %xyz_DG(3,:,:,:) double {mustBeReal}
      X(:,:,:)        double {mustBeReal}
      Y(:,:,:)        double {mustBeReal}
      Z(:,:,:)        double {mustBeReal}
      %xyz_DN(3,:,:,:) double {mustBeReal}
      %FaceNormal_DDF(3,3,:,:,:) double {mustBeReal}
      CellSize_D(3,1)  double {mustBeGreaterThan(CellSize_D,0)} = 1
      CellVolume_G(:,:,:) double {mustBeGreaterThan(CellVolume_G,0)}
   end
   
   %======================== CONSTRUCTORS ============================
   methods
      function obj = Grid
         %grid Construct an instance of this class
         %   Detailed explanation goes here
         
         switch Parameters.GridType
            case 'Cartesian'
               disp('Using Cartesian grid')
               
               nG = Parameters.nG;
               GridSize = Parameters.GridSize;
               FullSize = Parameters.FullSize;
               iMinAll = Parameters.iMinAll;
               iMaxAll = Parameters.iMaxAll;
               jMinAll = Parameters.jMinAll;
               jMaxAll = Parameters.jMaxAll;               
               kMinAll = Parameters.kMinAll;
               kMaxAll = Parameters.kMaxAll;
               
               obj.CoordMinMax_D = Parameters.xyzMinMax;
               obj.CellSize_D = (obj.CoordMinMax_D(:,2) - ...
                  obj.CoordMinMax_D(:,1)) ./ Parameters.GridSize';
               
%                obj.xyz_DG = Inf([3,FullSize]);
%                obj.xyz_DN = Inf([3,FullSize]);
%                for k = kMinAll:kMaxAll
%                for j = jMinAll:jMaxAll
%                for i = iMinAll:iMaxAll
%                   obj.xyz_DG(:,i,j,k) = obj.CoordMinMax_D(:,1) + ...
%                      ([i-nG,j-nG,k-nG]'- 0.5) .* obj.CellSize_D;
%                end;end;end
         
               x = linspace(...
                  obj.CoordMinMax_D(1,1) - (0.5 - nG)*obj.CellSize_D(1),...
                  obj.CoordMinMax_D(1,2) + (nG - 0.5)*obj.CellSize_D(1),...
                  FullSize(1));
               y = linspace(...
                  obj.CoordMinMax_D(1,1) - (0.5 - nG)*obj.CellSize_D(2),...
                  obj.CoordMinMax_D(1,2) + (nG - 0.5)*obj.CellSize_D(2),...
                  FullSize(2));               
               z = linspace(...
                  obj.CoordMinMax_D(1,1) - (0.5 - nG)*obj.CellSize_D(3),...
                  obj.CoordMinMax_D(1,2) + (nG - 0.5)*obj.CellSize_D(3),...
                  FullSize(3));
               
               [obj.X,obj.Y,obj.Z] = ndgrid(x,y,z);
         
%                for k = kMinAll:kMaxAll
%                for j = jMinAll:jMaxAll
%                for i = iMinAll:iMaxAll
%                   obj.xyz_DN(:,i,j,k) = obj.CoordMinMax_D(:,1) + ...
%                      ([i-nG,j-nG,k-nG]'- 1.0) .* obj.CellSize_D;
                  
                  % FaceNormal is probably not needed in Cartesian
%                end;end;end
          
            case 'Spherical'
               error('Not yet implemented!')
            otherwise
               error('Unknown grid type=%s',Parameters.Gridtype)
         end

      end
      
      function x = getX(obj)
      %GETX Get the x coordinates.
         nG = Parameters.nG;
         x = obj.X(nG+1:end-nG,nG+1,nG+1);
      end
   end
end

