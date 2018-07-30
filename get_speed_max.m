function [Cmax_XG,Cmax_YG,Cmax_ZG] = get_speed_max(state_VG)
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

Cs2_G = gamma* state_VG(P_,:,:,:) ./ state_VG(Rho_,:,:,:);
Ca2_G = (state_VG(Bx_,:,:,:).^2 + state_VG(By_,:,:,:).^2 + ...
   state_VG(Bz_,:,:,:).^2) ./ state_VG(Rho_,:,:,:);
Can2_XG = (state_VG(Bx_,:,:,:).^2 ) ./ state_VG(Rho_,:,:,:);
Can2_YG = (state_VG(By_,:,:,:).^2 ) ./ state_VG(Rho_,:,:,:);
Can2_ZG = (state_VG(Bz_,:,:,:).^2 ) ./ state_VG(Rho_,:,:,:);

Cmax_XG = abs(state_VG(Ux_,:,:,:)./state_VG(Rho_,:,:,:)) + ...
   sqrt(0.5*(Cs2_G + Ca2_G + sqrt((Cs2_G + Ca2_G).^2 - 4*Cs2_G.*Can2_XG)));
Cmax_YG = abs(state_VG(Uy_,:,:,:)./state_VG(Rho_,:,:,:)) + ...
   sqrt(0.5*(Cs2_G + Ca2_G + sqrt((Cs2_G + Ca2_G).^2 - 4*Cs2_G.*Can2_YG)));
Cmax_ZG = abs(state_VG(Uz_,:,:,:)./state_VG(Rho_,:,:,:)) + ...
   sqrt(0.5*(Cs2_G + Ca2_G + sqrt((Cs2_G + Ca2_G).^2 - 4*Cs2_G.*Can2_ZG)));


end

