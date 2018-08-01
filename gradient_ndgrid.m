function varargout = gradient_ndgrid(varargin)
%GRADIENT_NDGRID Gradient of vector
%   Not currently used or fully tested! Be careful!

narginchk(4,4);
f = varargin{1};
x = varargin{1};
y = varargin{2};
z = varargin{3};

% Loop over each dimension. 

varargout = cell(1,ndim);
siz = size(f);
% first dimension 
g  = zeros(size(f),class(f)); % case of singleton dimension
h = loc{1}(:); 
n = siz(1);

% Take forward differences on left and right edges
if n > 1
   g(1,:) = (f(2,:) - f(1,:))/(h(2)-h(1));
   g(n,:) = (f(n,:) - f(n-1,:))/(h(end)-h(end-1));
end

% Take centered differences on interior points
if n > 2
   g(2:n-1,:) = (f(3:n,:)-f(1:n-2,:)) ./ (h(3:n) - h(1:n-2));
end


end

