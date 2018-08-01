classdef FaceValue < handle
   %FaceValue Class of left and right state on the faces
   %   Detailed explanation goes here
   
   properties
      LState_VX(:,:,:,:)  double {mustBeReal}
      RState_VX(:,:,:,:)  double {mustBeReal}
      LState_VY(:,:,:,:)  double {mustBeReal}
      RState_VY(:,:,:,:)  double {mustBeReal}
      LState_VZ(:,:,:,:)  double {mustBeReal}
      RState_VZ(:,:,:,:)  double {mustBeReal}
   end
   
   methods
      function obj = FaceValue
         %FaceValue Construct an instance of this class
         %
         
%          nVar = Parameters.nVar;
%          GridSize = Parameters.GridSize;
%          
%          obj.LState_VX = Inf([nVar,GridSize+1]);
%          obj.RState_VX = Inf([nVar,GridSize+1]);
%          obj.LState_VY = Inf([nVar,GridSize+1]);
%          obj.RState_VY = Inf([nVar,GridSize+1]);
%          obj.LState_VZ = Inf([nVar,GridSize+1]);
%          obj.RState_VZ = Inf([nVar,GridSize+1]);
      end
      
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
                  dqR_X = squeeze(state_VG(:,iMin+1:iMax+2,jMin:jMax,kMin:kMax) - state_VG(:,iMin:iMax+1,jMin:jMax,kMin:kMax));
                  dqL_X = squeeze(state_VG(:,iMin:iMax+1,jMin:jMax,kMin:kMax) - state_VG(:,iMin-1:iMax,jMin:jMax,kMin:kMax));
                  
                  dq_X = minmod(dqR_X,dqL_X);
            end
            
            % Get the gradient of states
            
            % Linear interpolation onto edge centers (with limiters)
            
            
         end
         
      end
      
      
      function [Flux_VX,Flux_VY,Flux_VZ] = get_physical_flux(obj)
         % Calculate the physical fluxes
         %INPUT:
         %  FaceValue Class object
         %OUTPUT:
         %  Fluxes in 3 directions
         
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
         
         LFlux_VX = Inf([nVar,GridSize+[1 0 0]]); 
         RFlux_VX = Inf([nVar,GridSize+[1 0 0]]);
         LFlux_VY = Inf([nVar,GridSize+[0 1 0]]); 
         RFlux_VY = Inf([nVar,GridSize+[0 1 0]]);
         LFlux_VZ = Inf([nVar,GridSize+[0 0 1]]); 
         RFlux_VZ = Inf([nVar,GridSize+[0 0 1]]);
         
         
         % Density flux
         LFlux_VX(Rho_,:,:,:) = obj.LState_VX(Ux_,:,:,:);
         RFlux_VX(Rho_,:,:,:) = obj.RState_VX(Ux_,:,:,:);           
         LFlux_VY(Rho_,:,:,:) = obj.LState_VY(Uy_,:,:,:);
         RFlux_VY(Rho_,:,:,:) = obj.RState_VY(Uy_,:,:,:);         
         LFlux_VZ(Rho_,:,:,:) = obj.LState_VZ(Uz_,:,:,:);
         RFlux_VZ(Rho_,:,:,:) = obj.RState_VZ(Uz_,:,:,:);
         
         
         % Momentum flux
         LFlux_VX(Ux_,:,:,:) = ...
            obj.LState_VX(Ux_,:,:,:).^2 ./ obj.LState_VX(Rho_,:,:,:) + ...
            obj.LState_VX(P_,:,:,:) + ...
            0.5*sum(obj.LState_VX(B_,:,:,:).^2,1) - ...
            obj.LState_VX(Bx_,:,:,:).^2;
         RFlux_VX(Ux_,:,:,:) = ...
            obj.RState_VX(Ux_,:,:,:).^2 ./ obj.RState_VX(Rho_,:,:,:) + ...
            obj.RState_VX(P_,:,:,:) + ...
            0.5*sum(obj.RState_VX(B_,:,:,:).^2,1) - ...
            obj.RState_VX(Bx_,:,:,:).^2;         
         LFlux_VX(Uy_,:,:,:) = obj.LState_VX(Ux_,:,:,:) .* ...
            obj.LState_VX(Uy_,:,:,:) ./ obj.LState_VX(Rho_,:,:,:) - ...
            obj.LState_VX(Bx_,:,:,:).*obj.LState_VX(By_,:,:,:);
         RFlux_VX(Uy_,:,:,:) = obj.RState_VX(Ux_,:,:,:) .* ...
            obj.RState_VX(Uy_,:,:,:) ./ obj.RState_VX(Rho_,:,:,:) - ...
            obj.RState_VX(Bx_,:,:,:).*obj.RState_VX(By_,:,:,:);
         LFlux_VX(Uz_,:,:,:) = obj.LState_VX(Ux_,:,:,:) .* ...
            obj.LState_VX(Uz_,:,:,:) ./ obj.LState_VX(Rho_,:,:,:)- ...
            obj.LState_VX(Bx_,:,:,:).*obj.LState_VX(Bz_,:,:,:);
         RFlux_VX(Uz_,:,:,:) = obj.RState_VX(Ux_,:,:,:) .* ...
            obj.RState_VX(Uz_,:,:,:) ./ obj.RState_VX(Rho_,:,:,:) - ...
            obj.RState_VX(Bx_,:,:,:).*obj.RState_VX(Bz_,:,:,:);
         
        
         LFlux_VY(Ux_,:,:,:) = obj.LState_VY(Uy_,:,:,:) .* ...
            obj.LState_VY(Ux_,:,:,:) ./ obj.LState_VY(Rho_,:,:,:) - ...
            obj.LState_VY(By_,:,:,:).*obj.LState_VY(Bx_,:,:,:);
         RFlux_VY(Ux_,:,:,:) = obj.RState_VY(Uy_,:,:,:) .* ...
            obj.RState_VY(Ux_,:,:,:) ./ obj.RState_VY(Rho_,:,:,:) - ...
            obj.RState_VY(By_,:,:,:).*obj.RState_VY(Bx_,:,:,:);
         LFlux_VY(Uy_,:,:,:) = ...
            obj.LState_VY(Uy_,:,:,:).^2 ./ obj.LState_VY(Rho_,:,:,:) + ...
            obj.LState_VY(P_,:,:,:) + ...
            0.5*sum(obj.LState_VY(B_,:,:,:).^2,1) - ...
            obj.LState_VY(By_,:,:,:).^2;
         RFlux_VY(Uy_,:,:,:) = ...
            obj.RState_VY(Uy_,:,:,:).^2 ./ obj.RState_VY(Rho_,:,:,:) + ...
            obj.RState_VY(P_,:,:,:) + ...
            0.5*sum(obj.RState_VY(B_,:,:,:).^2,1) - ...
            obj.RState_VY(By_,:,:,:).^2;
         LFlux_VY(Uz_,:,:,:) = obj.LState_VY(Uy_,:,:,:) .* ...
            obj.LState_VY(Uz_,:,:,:) ./ obj.LState_VY(Rho_,:,:,:)- ...
            obj.LState_VY(Bx_,:,:,:).*obj.LState_VY(Bz_,:,:,:);
         RFlux_VY(Uz_,:,:,:) = obj.RState_VY(Uy_,:,:,:) .* ...
            obj.RState_VY(Uz_,:,:,:) ./ obj.RState_VY(Rho_,:,:,:)- ...
            obj.RState_VY(Bx_,:,:,:).*obj.RState_VY(Bz_,:,:,:);         
         
         
         LFlux_VZ(Ux_,:,:,:) = obj.LState_VZ(Uz_,:,:,:) .* ...
            obj.LState_VZ(Ux_,:,:,:) ./ obj.LState_VZ(Rho_,:,:,:)- ...
            obj.LState_VZ(Bz_,:,:,:).*obj.LState_VZ(Bx_,:,:,:);
         RFlux_VZ(Ux_,:,:,:) = obj.RState_VZ(Uz_,:,:,:) .* ...
            obj.RState_VZ(Ux_,:,:,:) ./ obj.RState_VZ(Rho_,:,:,:)- ...
            obj.RState_VZ(Bz_,:,:,:).*obj.RState_VZ(Bx_,:,:,:);
         LFlux_VZ(Uy_,:,:,:) = obj.LState_VZ(Uz_,:,:,:) .* ...
            obj.LState_VZ(Uy_,:,:,:) ./ obj.LState_VZ(Rho_,:,:,:)- ...
            obj.LState_VZ(Bz_,:,:,:).*obj.LState_VZ(By_,:,:,:);
         RFlux_VZ(Uy_,:,:,:) = obj.RState_VZ(Uz_,:,:,:) .* ...
            obj.RState_VZ(Uy_,:,:,:) ./ obj.RState_VZ(Rho_,:,:,:)- ...
            obj.RState_VZ(Bz_,:,:,:).*obj.RState_VZ(By_,:,:,:);         
         LFlux_VZ(Uz_,:,:,:) = ...
            obj.LState_VZ(Uz_,:,:,:).^2 ./ obj.LState_VZ(Rho_,:,:,:) + ...
            obj.LState_VZ(P_,:,:,:) + ...
            0.5*sum(obj.LState_VZ(B_,:,:,:).^2,1) - ...
            obj.LState_VZ(Bz_,:,:,:).^2;
         RFlux_VZ(Uz_,:,:,:) = ...
            obj.RState_VZ(Uz_,:,:,:).^2 ./ obj.RState_VZ(Rho_,:,:,:) + ...
            obj.RState_VZ(P_,:,:,:) + ...
            0.5*sum(obj.RState_VZ(B_,:,:,:).^2,1) - ...
            obj.RState_VZ(Bz_,:,:,:).^2;          
         

         % Magnetic flux
         LFlux_VX(Bx_,:,:,:) = 0;
         RFlux_VX(Bx_,:,:,:) = 0;
         LFlux_VX(By_,:,:,:) = obj.LState_VX(Ux_,:,:,:).*...
            obj.LState_VX(By_,:,:,:) - obj.LState_VX(Bx_,:,:,:).*...
            obj.LState_VX(Uy_,:,:,:);
         RFlux_VX(By_,:,:,:) = obj.RState_VX(Ux_,:,:,:).*...
            obj.RState_VX(By_,:,:,:) - obj.RState_VX(Bx_,:,:,:).*...
            obj.RState_VX(Uy_,:,:,:); 
         LFlux_VX(Bz_,:,:,:) = obj.LState_VX(Ux_,:,:,:).*...
            obj.LState_VX(Bz_,:,:,:) - obj.LState_VX(Bx_,:,:,:).*...
            obj.LState_VX(Uz_,:,:,:);
         RFlux_VX(Bz_,:,:,:) = obj.RState_VX(Ux_,:,:,:).*...
            obj.RState_VX(Bz_,:,:,:) - obj.RState_VX(Bx_,:,:,:).*...
            obj.RState_VX(Uz_,:,:,:);  
         
         
         LFlux_VY(Bx_,:,:,:) = obj.LState_VY(Uy_,:,:,:).*...
            obj.LState_VY(Bx_,:,:,:) - obj.LState_VY(By_,:,:,:).*...
            obj.LState_VY(Ux_,:,:,:);
         RFlux_VY(Bx_,:,:,:) = obj.RState_VY(Uy_,:,:,:).*...
            obj.RState_VY(Bx_,:,:,:) - obj.RState_VY(By_,:,:,:).*...
            obj.RState_VY(Ux_,:,:,:);
         LFlux_VY(By_,:,:,:) = 0;
         RFlux_VY(By_,:,:,:) = 0;
         LFlux_VY(Bz_,:,:,:) = obj.LState_VY(Uy_,:,:,:).*...
            obj.LState_VY(Bz_,:,:,:) - obj.LState_VY(By_,:,:,:).*...
            obj.LState_VY(Uz_,:,:,:);
         RFlux_VY(Bz_,:,:,:) = obj.RState_VY(Uy_,:,:,:).*...
            obj.RState_VY(Bz_,:,:,:) - obj.RState_VY(By_,:,:,:).*...
            obj.RState_VY(Uz_,:,:,:);
         
         
         LFlux_VZ(Bx_,:,:,:) = obj.LState_VZ(Uz_,:,:,:).*...
            obj.LState_VZ(Bx_,:,:,:) - obj.LState_VZ(Bz_,:,:,:).*...
            obj.LState_VZ(Ux_,:,:,:);
         RFlux_VZ(Bx_,:,:,:) = obj.RState_VZ(Uz_,:,:,:).*...
            obj.RState_VZ(Bx_,:,:,:) - obj.RState_VZ(Bz_,:,:,:).*...
            obj.RState_VZ(Ux_,:,:,:);
         LFlux_VZ(By_,:,:,:) = obj.LState_VZ(Uz_,:,:,:).*...
            obj.LState_VZ(By_,:,:,:) - obj.LState_VZ(Bz_,:,:,:).*...
            obj.LState_VZ(Uy_,:,:,:);
         RFlux_VZ(By_,:,:,:) = obj.RState_VZ(Uz_,:,:,:).*...
            obj.RState_VZ(By_,:,:,:) - obj.RState_VZ(Bz_,:,:,:).*...
            obj.RState_VZ(Uy_,:,:,:);
         LFlux_VZ(Bz_,:,:,:) = 0;
         RFlux_VZ(Bz_,:,:,:) = 0;          
         
         % Pressure flux / energy flux
         if ~Parameters.UseConservative
            LFlux_VX(P_,:,:,:) = obj.LState_VX(P_,:,:,:).*...
               obj.LState_VX(Ux_,:,:,:)./obj.LState_VX(Rho_,:,:,:);
            RFlux_VX(P_,:,:,:) = obj.RState_VX(P_,:,:,:).*...
               obj.RState_VX(Ux_,:,:,:)./obj.RState_VX(Rho_,:,:,:);
            LFlux_VY(P_,:,:,:) = obj.LState_VY(P_,:,:,:).*...
               obj.LState_VY(Uy_,:,:,:)./obj.LState_VY(Rho_,:,:,:);
            RFlux_VY(P_,:,:,:) = obj.RState_VY(P_,:,:,:).*...
               obj.RState_VY(Uy_,:,:,:)./obj.RState_VY(Rho_,:,:,:);
            LFlux_VZ(P_,:,:,:) = obj.LState_VZ(P_,:,:,:).*...
               obj.LState_VZ(Uz_,:,:,:)./obj.LState_VZ(Rho_,:,:,:);
            RFlux_VZ(P_,:,:,:) = obj.RState_VZ(P_,:,:,:).*...
               obj.RState_VZ(Uz_,:,:,:)./obj.RState_VZ(Rho_,:,:,:);
         else
            % Currently I am using the same index for pressure/energy
            gamma = Const.gamma;
            
            LFlux_VX(E_,:,:,:) = obj.LState_VX(Ux_,:,:,:)./...
               obj.LState_VX(Rho_,:,:,:).* ...
               ((obj.LState_VX(P_,:,:,:) / (gamma-1) + ...
               0.5 ./ obj.LState_VX(Rho_,:,:,:).*...
               sum(obj.LState_VX(U_,:,:,:).^2,1) + ...
               0.5 .* sum(obj.LState_VX(B_,:,:,:).^2,1) + ...
               obj.LState_VX(P_,:,:,:) + ...
               0.5 .* sum(obj.LState_VX(B_,:,:,:).^2,1))) - ...
               sum(obj.LState_VX(U_,:,:,:).*obj.LState_VX(B_,:,:,:),1).*...
               obj.LState_VX(Bx_,:,:,:);
  
            RFlux_VX(E_,:,:,:) = obj.RState_VX(Ux_,:,:,:)./...
               obj.RState_VX(Rho_,:,:,:).* ...
               ((obj.RState_VX(P_,:,:,:) / (gamma-1) + ...
               0.5 ./ obj.RState_VX(Rho_,:,:,:).*...
               sum(obj.RState_VX(U_,:,:,:).^2,1) + ...
               0.5 .* sum(obj.RState_VX(B_,:,:,:).^2,1) + ...
               obj.RState_VX(P_,:,:,:) + ...
               0.5 .* sum(obj.RState_VX(B_,:,:,:).^2,1))) - ...
               sum(obj.RState_VX(U_,:,:,:).*obj.RState_VX(B_,:,:,:),1).*...
               obj.RState_VX(Bx_,:,:,:);
            
            LFlux_VY(E_,:,:,:) = obj.LState_VY(Uy_,:,:,:)./...
               obj.LState_VY(Rho_,:,:,:).* ...
               ((obj.LState_VY(P_,:,:,:) / (gamma-1) + ...
               0.5 ./ obj.LState_VY(Rho_,:,:,:).*...
               sum(obj.LState_VY(U_,:,:,:).^2,1) + ...
               0.5 .* sum(obj.LState_VY(B_,:,:,:).^2,1) + ...
               obj.LState_VY(P_,:,:,:) + ...
               0.5 .* sum(obj.LState_VY(B_,:,:,:).^2,1))) - ...
               sum(obj.LState_VY(U_,:,:,:).*obj.LState_VY(B_,:,:,:),1).*...
               obj.LState_VY(By_,:,:,:);
  
            RFlux_VY(E_,:,:,:) = obj.RState_VY(Uy_,:,:,:)./...
               obj.RState_VY(Rho_,:,:,:).* ...
               ((obj.RState_VY(P_,:,:,:) / (gamma-1) + ...
               0.5 ./ obj.RState_VY(Rho_,:,:,:).*...
               sum(obj.RState_VY(U_,:,:,:).^2,1) + ...
               0.5 .* sum(obj.RState_VY(B_,:,:,:).^2,1) + ...
               obj.RState_VY(P_,:,:,:) + ...
               0.5 .* sum(obj.RState_VY(B_,:,:,:).^2,1))) - ...
               sum(obj.RState_VY(U_,:,:,:).*obj.RState_VY(B_,:,:,:),1).*...
               obj.RState_VY(By_,:,:,:);
           
           LFlux_VZ(E_,:,:,:) = obj.LState_VZ(Uz_,:,:,:)./...
               obj.LState_VZ(Rho_,:,:,:).* ...
               ((obj.LState_VZ(P_,:,:,:) / (gamma-1) + ...
               0.5 ./ obj.LState_VZ(Rho_,:,:,:).*...
               sum(obj.LState_VZ(U_,:,:,:).^2,1) + ...
               0.5 .* sum(obj.LState_VZ(B_,:,:,:).^2,1) + ...
               obj.LState_VZ(P_,:,:,:) + ...
               0.5 .* sum(obj.LState_VZ(B_,:,:,:).^2,1))) - ...
               sum(obj.LState_VZ(U_,:,:,:).*obj.LState_VZ(B_,:,:,:),1).*...
               obj.LState_VZ(Bz_,:,:,:);
  
           RFlux_VZ(E_,:,:,:) = obj.RState_VZ(Uz_,:,:,:)./...
               obj.RState_VZ(Rho_,:,:,:).* ...
               ((obj.RState_VZ(P_,:,:,:) / (gamma-1) + ...
               0.5 ./ obj.RState_VZ(Rho_,:,:,:).*...
               sum(obj.RState_VZ(U_,:,:,:).^2,1) + ...
               0.5 .* sum(obj.RState_VZ(B_,:,:,:).^2,1) + ...
               obj.RState_VZ(P_,:,:,:) + ...
               0.5 .* sum(obj.RState_VZ(B_,:,:,:).^2,1))) - ...
               sum(obj.RState_VZ(U_,:,:,:).*obj.RState_VZ(B_,:,:,:),1).*...
               obj.RState_VZ(Bz_,:,:,:);                       
         end
         
         % Collect all the physical fluxes
         Flux_VX = 0.5 * (LFlux_VX + RFlux_VX);
         Flux_VY = 0.5 * (LFlux_VY + RFlux_VY);
         Flux_VZ = 0.5 * (LFlux_VZ + RFlux_VZ);
      end
      
      
      function [Flux_VX,Flux_VY,Flux_VZ] = add_numerical_flux(obj,...
            Flux_VX,Flux_VY,Flux_VZ)
         % Calculate numerical fluxes and add to the physical fluxes.
         %INPUTS:
         % FaceValue class object
         % physical face fluxes in 3 directions
         %OUTPUTS:
         % total face fluxes in 3 directions
         
         if Parameters.Scheme == 'Rusanov'
            
            [cmax_XF,cmax_YF,cmax_ZF] = obj.get_speed_max;
            
            if ~Parameters.UseConservative
               Flux_VX = Flux_VX - 0.5*cmax_XF.*(obj.RState_VX - obj.LState_VX);
               Flux_VY = Flux_VY - 0.5*cmax_YF.*(obj.RState_VY - obj.LState_VY);
               Flux_VZ = Flux_VZ - 0.5*cmax_ZF.*(obj.RState_VZ - obj.LState_VZ);       
            else
            % If I solve energy equation instead of pressure, there's
            % duplicate calculation above, even though the expression looks
            % compact. That's why I use an if-else statement.
               Rho_ = Parameters.Rho_;
               Bz_  = Parameters.Bz_;
               P_   = Parameters.P_;
               E_   = Parameters.E_;
               U_   = Parameters.U_;
               B_   = Parameters.B_;
               gamma= Const.gamma;
            
               Flux_VX(Rho_:Bz_,:,:,:) = Flux_VX(Rho_:Bz_,:,:,:) - ...
                  0.5*cmax_XF.*(obj.RState_VX(Rho_:Bz_,:,:,:) - ...
                  obj.LState_VX(Rho_:Bz_,:,:,:));
               Flux_VX(E_,:,:,:) = Flux_VX(E_,:,:,:) - ...
                  0.5*cmax_XF.* (...
                  (obj.RState_VX(P_,:,:,:) / (gamma-1) + ...
                  0.5./obj.RState_VX(Rho_,:,:,:).*...
                  sum(obj.RState_VX(U_,:,:,:).^2,1) + ...
                  0.5*sum(obj.RState_VX(B_,:,:,:).^2,1)) - ...
                  (obj.LState_VX(P_,:,:,:) / (gamma-1) + ...
                  0.5./obj.LState_VX(Rho_,:,:,:).*...
                  sum(obj.LState_VX(U_,:,:,:).^2,1) + ...
                  0.5*sum(obj.LState_VX(B_,:,:,:).^2,1)));               
               
               Flux_VY(Rho_:Bz_,:,:,:) = Flux_VY(Rho_:Bz_,:,:,:) - ...
                  0.5*cmax_YF.*(obj.RState_VY(Rho_:Bz_,:,:,:) - ...
                  obj.LState_VY(Rho_:Bz_,:,:,:));
               Flux_VY(E_,:,:,:) = Flux_VY(E_,:,:,:) - ...
                  0.5*cmax_YF.* (...
                  (obj.RState_VY(P_,:,:,:) / (gamma-1) + ...
                  0.5./obj.RState_VY(Rho_,:,:,:).*...
                  sum(obj.RState_VY(U_,:,:,:).^2,1) + ...
                  0.5*sum(obj.RState_VY(B_,:,:,:).^2,1)) - ...
                  (obj.LState_VY(P_,:,:,:) / (gamma-1) + ...
                  0.5./obj.LState_VY(Rho_,:,:,:).*...
                  sum(obj.LState_VY(U_,:,:,:).^2,1) + ...
                  0.5*sum(obj.LState_VY(B_,:,:,:).^2,1)));                             
               
               Flux_VZ(Rho_:Bz_,:,:,:) = Flux_VZ(Rho_:Bz_,:,:,:) - ...
                  0.5*cmax_ZF.*(obj.RState_VZ(Rho_:Bz_,:,:,:) - ...
                  obj.LState_VZ(Rho_:Bz_,:,:,:));   
               Flux_VZ(E_,:,:,:) = Flux_VZ(E_,:,:,:) - ...
                  0.5*cmax_ZF.* (...
                  (obj.RState_VZ(P_,:,:,:) / (gamma-1) + ...
                  0.5./obj.RState_VZ(Rho_,:,:,:).*...
                  sum(obj.RState_VZ(U_,:,:,:).^2,1) + ...
                  0.5*sum(obj.RState_VZ(B_,:,:,:).^2,1)) - ...
                  (obj.LState_VZ(P_,:,:,:) / (gamma-1) + ...
                  0.5./obj.LState_VZ(Rho_,:,:,:).*...
                  sum(obj.LState_VZ(U_,:,:,:).^2,1) + ...
                  0.5*sum(obj.LState_VZ(B_,:,:,:).^2,1)));               
            end
            
         end
      end
        
   end
   
   methods (Access = private)
      function mm = minmod(v)
         % Using Harten's generalized definition
         % minmod: zero if opposite sign, otherwise the one of smaller magnitude.
         %m=size(v,1); mm=zeros(size(v,2),1); s=sum(sign(v),2)/m; ids=find(abs(s)==1);
         %if(~isempty(ids)); mm(ids)=s(ids).*min(abs(v(ids,:)),[],2); end
         s = sum(sign(v))/numel(v);
         if abs(s)==1; mm = s*min(abs(v(:))); else mm=0; end
      end
      
      function [cmax_XF,cmax_YF,cmax_ZF] = get_speed_max(obj)
         gamma = Const.gamma;
         
         LS_VX = obj.LState_VX; RS_VX = obj.RState_VX;
         LS_VY = obj.LState_VY; RS_VY = obj.RState_VY;
         LS_VZ = obj.LState_VZ; RS_VZ = obj.RState_VZ;
         
         Rho_ = Parameters.Rho_;
         Ux_  = Parameters.Ux_;
         Uy_  = Parameters.Uy_;
         Uz_  = Parameters.Uz_;
         Bx_  = Parameters.Bx_;
         By_  = Parameters.By_;
         Bz_  = Parameters.Bz_;
         P_   = Parameters.P_;
         %U_   = Parameters.U_;
         %B_   = Parameters.B_;
         
         % Maybe use this for speed?
         %rho_VX = LState_VX(Rho_,:,:,:) + RState_VX(Rho_,:,:,:);
         
         Cs2_XF = gamma*(LS_VX(P_,:,:,:) + RS_VX(P_,:,:,:)) ./...
            (LS_VX(Rho_,:,:,:) + RS_VX(Rho_,:,:,:));
         Cs2_YF = gamma*(LS_VY(P_,:,:,:) + RS_VY(P_,:,:,:)) ./...
            (LS_VY(Rho_,:,:,:) + RS_VY(Rho_,:,:,:));
         Cs2_ZF = gamma*(LS_VZ(P_,:,:,:) + RS_VZ(P_,:,:,:)) ./...
            (LS_VZ(Rho_,:,:,:) + RS_VZ(Rho_,:,:,:));
         
         Ca2_XF = ( (LS_VX(Bx_,:,:,:) + RS_VX(Bx_,:,:,:)).^2 + ...
            (LS_VX(By_,:,:,:) + RS_VX(By_,:,:,:)).^2 + ...
            (LS_VX(Bz_,:,:,:) + RS_VX(Bz_,:,:,:)).^2 ) ./ ...
            (2 * LS_VX(Rho_,:,:,:) + RS_VX(Rho_,:,:,:));
         Ca2_YF = ( (LS_VY(Bx_,:,:,:) + RS_VY(Bx_,:,:,:)).^2 + ...
            (LS_VY(By_,:,:,:) + RS_VY(By_,:,:,:)).^2 + ...
            (LS_VY(Bz_,:,:,:) + RS_VY(Bz_,:,:,:)).^2 ) ./ ...
            (2 * LS_VY(Rho_,:,:,:) + RS_VY(Rho_,:,:,:));
         Ca2_ZF = ( (LS_VZ(Bx_,:,:,:) + RS_VZ(Bx_,:,:,:)).^2 + ...
            (LS_VZ(By_,:,:,:) + RS_VZ(By_,:,:,:)).^2 + ...
            (LS_VZ(Bz_,:,:,:) + RS_VZ(Bz_,:,:,:)).^2 ) ./ ...
            (2 * LS_VZ(Rho_,:,:,:) + RS_VZ(Rho_,:,:,:));
         
         Can2_XF = ( (LS_VX(Bx_,:,:,:) + RS_VX(Bx_,:,:,:)).^2 ) ./ ...
            (2 * LS_VX(Rho_,:,:,:) + RS_VX(Rho_,:,:,:));
         Can2_YF = ( (LS_VY(By_,:,:,:) + RS_VY(By_,:,:,:)).^2 ) ./ ...
            (2 * LS_VY(Rho_,:,:,:) + RS_VY(Rho_,:,:,:));
         Can2_ZF = ( (LS_VZ(Bz_,:,:,:) + RS_VZ(Bz_,:,:,:)).^2 ) ./ ...
            (2 * LS_VZ(Rho_,:,:,:) + RS_VZ(Rho_,:,:,:));
         
         cmax_XF = ...
            0.5 * abs(LS_VX(Ux_,:,:,:)./LS_VX(Rho_,:,:,:) + ...
            RS_VX(Ux_,:,:,:)./RS_VX(Rho_,:,:,:)) + ...
            sqrt( 0.5*(Cs2_XF + Ca2_XF + ...
            sqrt((Cs2_XF + Ca2_XF).^2-4*Cs2_XF.*Can2_XF)) );
         
         cmax_YF = ...
            0.5 * abs(LS_VY(Uy_,:,:,:)./LS_VY(Rho_,:,:,:) + ...
            RS_VY(Uy_,:,:,:)./RS_VY(Rho_,:,:,:)) + ...
            sqrt( 0.5*(Cs2_YF + Ca2_YF + ...
            sqrt((Cs2_YF + Ca2_YF).^2-4*Cs2_YF.*Can2_YF)) );
         
         cmax_ZF = ...
            0.5 * abs(LS_VZ(Uz_,:,:,:)./LS_VZ(Rho_,:,:,:) + ...
            RS_VZ(Uz_,:,:,:)./RS_VZ(Rho_,:,:,:)) + ...
            sqrt( 0.5*(Cs2_ZF + Ca2_ZF + ...
            sqrt((Cs2_ZF + Ca2_ZF).^2-4*Cs2_ZF.*Can2_ZF)) );
      end   
      
   end
   
end

