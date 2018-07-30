function div=divergence_ndgrid(varargin)
%DIVERGENCE_NDGRID Divergence of a vector field
%   DIV = DIVERGENCE_NDGRID(X,Y,Z,VAR) computes the divergence of a 3-D
%   vector field VAR. The arrays X,Y,Z define the coordinates for
%   U,V,W and must be monotonic and 3-D plaid (as if produced by
%   NDGRID).
%
% Hongyang Zhou, 06/03/2018

narginchk(4,4);
x = varargin{1};
y = varargin{2};
z = varargin{3};
var = varargin{4};

hx = x(:,1,1);
hy = y(1,:,1);
hz = z(1,1,:);

% [px, ~, ~] = gradient(squeeze(var(1,:,:,:)), hx, hy, hz);
% [~, qy, ~] = gradient(squeeze(var(2,:,:,:)), hx, hy, hz);
% [~, ~, rz] = gradient(squeeze(var(3,:,:,:)), hx, hy, hz);

% Take forward differences on left and right edges
% if n > 1
%    g(1,:) = (f(2,:) - f(1,:))/(h(2)-h(1));
%    g(n,:) = (f(n,:) - f(n-1,:))/(h(end)-h(end-1));
% end

siz = size(var); siz = siz(2:end);

px = zeros(siz,class(var));
qy = zeros(siz,class(var));
rz = zeros(siz,class(var));

n = size(hx,1);
% Right now do nothing for the ghost cells; maybe needed later!
% Take central differences on interior points
if n > 2
   px(2:n-1,:) = (squeeze(var(1,3:n,:)-var(1,1:n-2,:))) ./ ...
      (hx(3:n) - hx(1:n-2));
end

n = size(hy,2);
if n > 2
   qy(:,2:n-1,:) = (squeeze(var(2,:,3:n,:)-var(2,:,1:n-2,:))) ./ ...
      (hy(3:n) - hy(1:n-2));
end

n = size(hz,3);
if n > 2
   rz(:,:,2:n-1) = (squeeze(var(3,:,:,3:n)-var(3,:,:,1:n-2))) ./ ...
      (hz(3:n) - hz(1:n-2));
end

div = px+qy+rz;

end

