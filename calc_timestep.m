function [timestep] = calc_timestep(grid,state_GV)
%calc_timestep Summary of this function goes here
%   Detailed explanation goes here

CFL = Parameters.CFL;
CellSize_D = grid.CellSize_D;

[Cmax_GX,Cmax_GY,Cmax_GZ] = get_speed_max(state_GV);

% [Cmax_XF,Cmax_YF,Cmax_ZF] = faceValue.get_speed_max;
% 
% Cmax_XG = max(Cmax_XF(1,1:end-1,:,:),Cmax_XF(1,2:end,:,:));
% Cmax_YG = max(Cmax_YF(1,:,1:end-1,:),Cmax_YF(1,:,2:end,:));
% Cmax_ZG = max(Cmax_ZF(1,:,:,1:end-1),Cmax_ZF(1,:,:,2:end));

if Parameters.TimeAccurate
   timestep = CFL ./ ...
      (Cmax_GX/CellSize_D(1)+Cmax_GY/CellSize_D(2)+Cmax_GZ/CellSize_D(3));
   timestep = min(timestep(:)); 
else
   timestep = CFL / ...
      (Cmax_GX/CellSize_D(1)+Cmax_GY/CellSize_D(2)+Cmax_GZ/CellSize_D(3));
end

end

