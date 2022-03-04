function   [xnp1, yn] = oneStepSimulation(f, g, xn, un, ds, nevals, isDiscrete, varargin)
xnp1 = nan(ds.nx, nevals);
yn   = nan(ds.ny, nevals);
% f = @(x,u,nORt) sfuns.ffun(x, u, nORt, p, ds.nx, nevals);
% g = @(x,u,nORt) sfuns.gfun(x, u, nORt, p, ds.ny, nevals);

if isDiscrete
    n    = varargin{1};
    xnn  = varargin{2}; %Process noise @ step n
    ynn  = varargin{3}; %Observation noise @ step n
%     xnp1(:,:) = sfuns.ffun(xn, un, n, p, ds.nx, nevals);
%     yn(:,:)   = sfuns.gfun(xn, un, n, p, ds.ny, nevals);
    xnp1(:,:) = f(xn,un,n) + xnn; 
    yn(:,:)   = g(xn,un,n) + ynn; 
else %isContinuous TODO: noise has to be appropriately accounted for, 
     % for the continuous case too
    tn     = varargin{1};
    dt     = varargin{2};
    utpdt2 = varargin{3}; % u(t+dt/2)
    utpdt  = varargin{4}; % u(t+dt)
    % RK 34
    k1 = dt*f(xn,      un,     tn);
    k2 = dt*f(xn+k1/2, utpdt2, tn + dt/2);
    k3 = dt*f(xn+k2/2, utpdt2, tn + dt/2); 
    k4 = dt*f(xn+k3,   utpdt,  tn + dt); 
    xnp1(:,:) = xn + k1/6 + k2/3 + k3/3 + k4/6; %x(t+dt)
    yn(:,:)   = g(xn, un, tn);
    
end


end
