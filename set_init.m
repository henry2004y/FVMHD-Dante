function [density,velocity,Bfield,pressure] = set_init
%set_init Summary of this function goes here
%   Detailed explanation goes here

density  = ones(Parameters.FullSize);
velocity = zeros([3,Parameters.FullSize]);
Bfield   = zeros([3,Parameters.FullSize]);
pressure = ones(Parameters.FullSize);

nI = Parameters.nI;

switch Parameters.IC
   case 'contact discontinuity'
      velocity(:,:,:,:) = 0.0;
      density(1:nI/2,:,:) = 2.0;
      density(nI/2+1:end,:,:) = 1.0;
   case 'density wave'
      density(nI/2:nI/2,:,:) = 2.0;
      velocity(1,:,:,:) = 1.0;
      
      
   case 'shocktube'
      
   otherwise
      error('unknown initial condition type!')
end

end

