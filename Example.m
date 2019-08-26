%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 Example                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clc; 

% Set parameters 
type = 'BT';                 % example type (BT or UP)
n = 50;                      % image size 
savename = 'TestSol';        % savename 
maxIter = 10000;              % max. number of iterations
tol = 1e-6;                  % error tolerance 
maxIterProj = 10;            % max. number of iterations for projection 
tolProj = 1e-6;              % error tolerance for projection
alpha = 0.5;                 % alpha (for BT)
a = 5;                       % a (for UP)
epsilon = 0.4;               % epsilon (for UP)
example = '1To2Points';      % Example (see getExample.m)

%-------------------------------------------------------------------------
% Compute solution 
if ( strcmp(type,'BT') )
    [v,phi1,phi2,phi3,u] = BranchedTransportSolver(n,savename,maxIter,tol,maxIterProj,tolProj,alpha,example);
elseif ( strcmp(type,'UP') )
    [v,phi1,phi2,phi3,u] = UrbanPlanningSolver(n,savename,maxIter,tol,maxIterProj,tolProj,a,epsilon,example);
end
%-------------------------------------------------------------------------
% Show result 
figure(); imagesc(u); colormap gray; axis image; 
%-------------------------------------------------------------------------