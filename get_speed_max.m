function [Cmax_GX,Cmax_GY,Cmax_GZ] = get_speed_max(state_GV)
%GET_SPEED_MAX Get the cell-centered maximum speed.


gamma = Const.gamma;
Rho_ = Parameters.Rho_;
Ux_  = Parameters.Ux_;
Uy_  = Parameters.Uy_;
Uz_  = Parameters.Uz_;
Bx_  = Parameters.Bx_;
By_  = Parameters.By_;
Bz_  = Parameters.Bz_;
P_   = Parameters.P_;
U_   = Parameters.U_;
B_   = Parameters.B_;

Cs2_G = gamma* state_GV(:,:,:,P_) ./ state_GV(:,:,:,Rho_);
Ca2_G = (state_GV(:,:,:,Bx_).^2 + state_GV(:,:,:,By_).^2 + ...
   state_GV(:,:,:,Bz_).^2) ./ state_GV(:,:,:,Rho_);
Can2_XG = (state_GV(:,:,:,Bx_).^2 ) ./ state_GV(:,:,:,Rho_);
Can2_YG = (state_GV(:,:,:,By_).^2 ) ./ state_GV(:,:,:,Rho_);
Can2_ZG = (state_GV(:,:,:,Bz_).^2 ) ./ state_GV(:,:,:,Rho_);

Cmax_GX = abs(state_GV(:,:,:,Ux_)./state_GV(:,:,:,Rho_)) + ...
   sqrt(0.5*(Cs2_G + Ca2_G + sqrt((Cs2_G + Ca2_G).^2 - 4*Cs2_G.*Can2_XG)));
Cmax_GY = abs(state_GV(:,:,:,Uy_)./state_GV(:,:,:,Rho_)) + ...
   sqrt(0.5*(Cs2_G + Ca2_G + sqrt((Cs2_G + Ca2_G).^2 - 4*Cs2_G.*Can2_YG)));
Cmax_GZ = abs(state_GV(:,:,:,Uz_)./state_GV(:,:,:,Rho_)) + ...
   sqrt(0.5*(Cs2_G + Ca2_G + sqrt((Cs2_G + Ca2_G).^2 - 4*Cs2_G.*Can2_ZG)));


end

