function varargout = mySystSim(sfuns, systp, simp, varargin)
% simp = struct('nt', 'u', 'oneStep',  ) % Always Required
%               ('xnoise', 'ynoise')     % Only if NOISE is set to 1
%               ('t', 'udt2', 'Ts')      % Only if the system is continuous
%               ('ct')                   % Only if oneStep is true
%Set default parameters
NOISE = 0; PLOT = 0;
if isfield(simp, 'NOISE')
    NOISE = simp.NOISE;
end

%Update parameters is passed as arguments
if nargin > 3
    PLOT   = varargin{1};
end

% Update simulation parameters
ds = systp.dim;  % Aggregated dimensions

    if simp.oneStep == false
        nt  = simp.nt; % The user asked for the whole simulation
        u  = simp.u;     % Input (for all the times) 
        if systp.isDiscrete
            toi = simp.t; %Was 1:simp.nt
        else
            toi = simp.t;
        end
    else                 % The user asked for a one-step simulation
        nt = 2; % x(j), x(j+1)
        indTimes = simp.ct:simp.ct+1;
        u        = simp.u(:, indTimes);     % Input for the relevant samples
        if systp.isDiscrete
            toi = simp.t(indTimes); %Was = indTimes;
        else
            toi = simp.t(indTimes);
        end
    end
    
    %         ds.nt   = simp.nt;
    x0      = systp.x0;   % Initial conditions
    p       = systp.p;    % System's parameters
    if NOISE
        xnoise = simp.xnoise;        ynoise = simp.ynoise;
    else
        xnoise = zeros(ds.nx, nt);   ynoise = zeros(ds.ny, nt);
    end
    
    if ~systp.isDiscrete
        % Will only be used for the case of a continuous system
        utpdt2   = simp.udt2;  %% u(t+dt/2) (it is needed for RK integration)
        %         t        = simp.t;
        dt       = simp.Ts;
    end
    
    
    % Variables allocation
    x = nan(ds.nx, nt); x(:,1) = x0;
    y = nan(ds.ny, nt);
    
    % Simulation
    nevals = 1; %This parameter is not used at the moment
    f = @(x,u,nORt) sfuns.ffun(x, u, nORt, p, ds.nx, nevals);
    g = @(x,u,nORt) sfuns.gfun(x, u, nORt, p, ds.ny, nevals);

    if systp.isDiscrete %Discrete time system--------------------------
        oss   = @(j, xj, uj, xnj, ynj) oneStepSimulation(f, g, xj, uj, ds, nevals, systp.isDiscrete,...
            j, xnj, ynj);
        for iss = 2:nt %Index simulation step
%             if iss>=2101
%                 disp(iss)
%             end
            it = toi(iss-1);
            [x(:,iss), y(:,iss-1)] = oss(it, x(:,iss-1), u(:,iss-1), xnoise(:, iss-1), ynoise(:,iss-1));
        end
        [~, y(:,nt)] = oss(toi(end), x(:,nt), u(:,nt), xnoise(:, nt), ynoise(:,nt));
        
    else                %Continuous time system------------------------
        oss   = @(tj, xj, uj, ujpdt2, ujpdt) oneStepSimulation(f, g, xj, uj, ds, nevals, systp.isDiscrete,...
            tj, dt, ujpdt2, ujpdt);
        for iss = 2:nt
            it = toi(iss-1);
            [x(:,iss ), y(:,iss-1)] = oss(it,x(:,iss-1), u(:,iss-1), utpdt2(:,iss-1),   u(:,iss));
        end
        [~, y(:,nt)] = oss(toi(end), x(:,nt), u(:,nt), utpdt2(:,nt), u(:,nt));
    end
    out =  struct('x', x, 'y', y);


if PLOT
    % Plot
    figure,
    subplot(3,1,1), plot(simp.t, u), legend(systp.uNames)
    subplot(3,1,2), plot(simp.t, x), legend(systp.xNames);
    subplot(3,1,3), plot(simp.t, y), legend(systp.yNames)
end
% Return
varargout{1} = out;
% varargout{2} = ds;
end

