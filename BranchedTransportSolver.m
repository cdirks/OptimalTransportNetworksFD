%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Branched transport solver                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [v,phi1,phi2,phi3,u] = BranchedTransportSolver(n,savename,maxIter,tol,maxIterProj,tolProj,alpha,example)
disp(['+++++ Branched transport with alpha = ',num2str(alpha),' +++++']);
%-------------------------------------------------------------------------
% Create initial data and grid 
[v_init,range,x,y,bdy] = getExample(n,example); % TODO
[n1,n2,n3] = size(v_init);
[~,~,Z] = meshgrid(x,y,range);
hz = Z(1,1,2)-Z(1,1,1);
%-------------------------------------------------------------------------
% Compute gradient matrix 
gradientComplete = generateForwardGradient([n1,n2,n3],[1,1,1]);
gradX = gradientComplete{1};
gradY = gradientComplete{2};
gradZ = gradientComplete{3};
gradMatrix = [abs(gradX);abs(gradY);abs(gradZ)];
%-------------------------------------------------------------------------
% Parameters 
theta = 1; 
tauFactor = 2;
sigma = sum(gradMatrix,2); sigma = (1/max(sigma))*tauFactor;
tau = sum(gradMatrix,1); tau = (1/max(tau))/tauFactor;
%-------------------------------------------------------------------------
% Initial values 
v_ = v_init(:); v = v_init(:);
phi1 = zeros(size(v)); phi2 = phi1; phi3 = phi1; 
%-------------------------------------------------------------------------
% ITERATION 
error = 1; iter = 0;
while ( error > tol && iter <= maxIter )
    % Save old iterate
    vOld = v; 
    phi1Old = phi1; phi2Old = phi2; phi3Old = phi3; 
    % Update phi
    phi1 = phi1(:) + sigma*gradX*v_; phi2 = phi2(:) + sigma*gradY*v_; phi3 = phi3(:) + sigma*gradZ*v_;
    phi1 = reshape(phi1,size(v_init)); phi2 = reshape(phi2,size(v_init)); phi3 = reshape(phi3,size(v_init));
    [phi1,phi2] = MSprojectBTmex(phi1,phi2,alpha,hz,maxIterProj,tolProj); 
    phi3(phi3<0) = 0;
    % Update v
    v = v - tau * ( gradX'*phi1(:)+gradY'*phi2(:)+gradZ'*phi3(:) );
    v = max(0,min(1,v));
    v(bdy(:)) = v_init(bdy(:));
    % Overrelaxation 
    v_ = v + theta * (v - vOld);
    % Compute primal dual residual 
    if ( mod(iter,100)==0 ) 
        P = 1/tau*(vOld-v) - (gradX'*(phi1Old(:)-phi1(:))+gradY'*(phi2Old(:)-phi2(:))+gradZ'*(phi3Old(:)-phi3(:)));
        D1 = 1/sigma*(phi1Old(:)-phi1(:)) - gradX*(vOld-v);
        D2 = 1/sigma*(phi2Old(:)-phi2(:)) - gradY*(vOld-v);
        D3 = 1/sigma*(phi3Old(:)-phi3(:)) - gradZ*(vOld-v);
        primalRes = sum(abs(P))/numel(v);
        dualRes = (sum(abs(D1))+sum(abs(D2))+sum(abs(D3)))/(3*numel(phi1));
        error = primalRes+dualRes;
        disp(['Iteration ',num2str(iter),': Error = ',num2str(error)]);
    end
    iter = iter+1;
end
%-------------------------------------------------------------------------
% Save results and delift 
v = reshape(v,size(v_init));
thresh = sum(v(:))/numel(v);
vThresh = v; vThresh(v>=thresh) = 1; vThresh(v<thresh) = 0;
u = delift(vThresh,range); 
save(savename,'v','phi1','phi2','phi3','u');
%-------------------------------------------------------------------------
end