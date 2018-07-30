function [Rho,U,B,P,tEnd] = set_init_Riemann(grid)
%SET_INIT_RIEMANN Load the IC of a classical 1D Riemann Problems.
%   Detailed explanation goes here



%% Riemann Problems
switch Parameters.RiemannProblemType
    case{1} % Configuration 1, Sod's Problem
        fprintf('Case 1: Sods problem \n');
        p   = [1    0.1  ];
        u   = [0    0    ];
        rho = [1    0.125];
        tEnd = 0.1; cfl = 0.90;
        
    case{2} % Configuration 2, Left Expansion and right strong shock
        fprintf('Case 2: Left Expansion and right strong shock \n');
        p   = [1000 0.1  ];
        u   = [0    0    ];
        rho = [3    2    ];
        tEnd = 0.02; cfl = 0.90;
        
    case{3} % Configuration 3, Right Expansion and left strong shock
        fprintf('Case 3: Right Expansion and left strong shock \n');
        p   = [7    10   ];
        u   = [0    0    ];
        rho = [1    1    ];
        tEnd = 0.1; cfl = 0.90;
        
    case{4} % Configuration 4, Double Shock
        fprintf('Case 4: Double Shock \n');
        p   = [450  45   ];
        u   = [20   -6   ];
        rho = [6    6    ];
        tEnd = 0.02; cfl = 0.90;
        
    case{5} % Configuration 5, Double Expansion
        fprintf('Case 5: Double Expansion \n');
        p   = [40   40   ];
        u   = [-2   2    ];
        rho = [1    2.5  ];
        tEnd = 0.03; cfl = 0.90;

    case{6} % Configuration 6, Cavitation
        fprintf('Case 6: Cavitation \n');
        p   = [0.4  0.4  ];
        u   = [-2    2   ];
        rho = [ 1    1   ];
        tEnd = 0.1; cfl = 0.90;
        
    case{7} % Shocktube problem of G.A. Sod, JCP 27:1, 1978 
        fprintf('Shocktube problem of G.A. Sod, JCP 27:1, 1978');
        p   = [1.0  0.1  ];
        u   = [0.75 0    ];
        rho = [1    0.125];
        tEnd = 0.17; cfl = 0.90;
        
    case{8} % Lax test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997
        fprintf('Lax test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997');
        p   = [3.528 0.571];
        u   = [0.698 0    ];
        rho = [0.445 0.5  ];
        tEnd = 0.15; cfl = 0.90; 
      
    case{9} % Mach = 3 test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997
        fprintf('Mach = 3 test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997');
        p   = [10.333  1  ];
        u   = [ 0.92  3.55];
        rho = [ 3.857  1  ];
        tEnd = 0.09; cfl = 0.90;
        
    case{10} % Shocktube problem with supersonic zone
        fprintf('Shocktube problem with supersonic zone');
        p   = [1  0.02];
        u   = [0  0.00];
        rho = [1  0.02];
        tEnd = 0.162; cfl = 0.90; 

    case{11} % Contact discontinuity
        fprintf('Contact discontinuity');
        p   = [0.5 0.5];
        u   = [0.0 0.0];
        rho = [1.0 0.6];
        tEnd = 1; cfl = 0.90; 

    case{12} % Stationary shock
        fprintf('Stationary shock');
        p   = [ 1.0  0.1 ];
        u   = [-2.0 -2.0 ];
        rho = [ 1.0 0.125];
        tEnd = 0.1; cfl = 0.28; 
        
    case{13} % left side of 2-d Riemman case 12
        fprintf('left side of 2-d Riemman case 13');
        p   = [ 1.0  1.0 ];
        u   = [ 0.0  0.0 ];
        rho = [ 0.8  1.0 ];
        tEnd = 0.1; cfl = 0.25; 
        
	 case{14} % right side of 2-d Riemman case 12
        fprintf('right side of 2-d Riemman case 14');
        p   = [ 1.0  0.4 ];
        u   = [ 0.7276  0.0 ];
        rho = [ 1.0  0.5313 ];
        tEnd = 0.1; cfl = 0.25; 
        
    case{15} % Stationary Shock
        fprintf('right side of 2-d Riemman case 15');
        p   = [ 0.1  0.676 ];
        u   = [ 1.2  0.723966942148760 ];
        rho = [ 1.0  1.657534246575342 ];
        tEnd = 0.1; cfl = 0.25;
        
    otherwise 
        error('Case not available');
        
end
% Print for Riemann Problems
fprintf('\n');
fprintf('density (L): %1.3f\n',rho(1));
fprintf('velocity(L): %1.3f\n',u(1));
fprintf('Presure (L): %1.3f\n',p(1));
fprintf('\n');
fprintf('density (R): %1.3f\n',rho(2));
fprintf('velocity(R): %1.3f\n',u(2));
fprintf('Presure (R): %1.3f\n',p(2));
fprintf('\n');

%% Load Selected case Initial condition:
% Pre-Allocate variables
Rho = ones(Parameters.FullSize);
U   = zeros([3,Parameters.FullSize]);
B   = zeros([3,Parameters.FullSize]);
P   = ones(Parameters.FullSize);

% Parameters of regions dimensions
x_middle = (Parameters.xyzMinMax(1,2) - Parameters.xyzMinMax(1,1))/2;
L = find(grid.X<=x_middle);
R = find(grid.X>x_middle);

% Initial Condition for our 1D domain
% Density
Rho(L) = rho(1); % region 1
Rho(R) = rho(2); % region 2
% Velocity in x
U(L) = u(1); % region 1
U(R) = u(2); % region 2
% Pressure
P(L) = p(1); % region 1
P(R) = p(2); % region 2

end

