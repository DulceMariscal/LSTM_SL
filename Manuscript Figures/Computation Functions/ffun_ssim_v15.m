function xnp1 = ffun_ssim_v15(xn, un, n, p, nx, nevals) %Returns x(n+1) from x(n), for many particles
% 5-state model------------------------
% xn(1) = x-active 1
% xn(2) = x-active 2
% xn(3) = x-active 3
% xn(4) = x-active 4
% xn(5) = w-positive-tuned motor primitive
% xn(6) = w-negative-tuned motor primitive
% xn(7) = memory of errors
%--------------------------------------
ix = struct('uf',1,'f',2,'s',3,'us',4,'zp',5,'zn',6,'ip',7,'in',8,'lp',9,'ln',10);

% Parameters extraction
% kcell = num2cell(p);
% [ar, aa, aet, bet, br, ba_fac, ba_int, d] = kcell{:};

auf   = p(1);  af     = p(2);     as  = p(3);   aus  = p(4);
b0uf  = p(5);  b0f    = p(6);     b0s = p(7);   b0us = p(8);
cuf   = p(9);  cf     = p(10);    cs  = p(11);   cus = p(12);
aet   = p(13); bet    = p(14);
mu    = p(15); sigma1 = p(16); sigma2 = p(17);
bias_uf  = p(18);
bias_f   = p(19);
bias_s   = p(20);
bias_us  = p(21);
ae = p(22); be = p(23);
th = p(24);
al = p(25); bl = p(26);
d        = p(27);

IPE   = 5;   INE = 6;
maxPA = 1; %This is a fixed parameter and won't be fit



% Initializations
xnp1           =  nan(nx, nevals);

% Computation of relevant quantities
z  = xn(1,:) + xn(2,:) + xn(3,:) + 0*xn(4,:); % Motor output for current step
e  = d*un(1,:) - z;                         % Motor error for current step (~SLA)
em = abs(e);                                % Error magnitude

% Update rules
% 1. Reactive component update
%     xnp1(1,:)  =  ar*xn(1,:) + br*e;

% 2. Active components update
xnp1(1,:)  =  auf*xn(1,:);
xnp1(2,:)  =  af*xn(2,:);
xnp1(3,:)  =  as*xn(3,:);
xnp1(4,:)  =  aus*xn(4,:);

%% fvu .0567. BIC: -40677
if e > 0
    xnp1(1,:) = xnp1(1,:)  +  (b0uf     * (1  +  20*max(0,xn(IPE,:)       + (cuf)*xn(INE,:))))*e;    % Active UltraFast, during positive errors
    xnp1(2,:) = xnp1(2,:)  +  (b0f      * (1  +  20*max(0,xn(IPE,:)       + (cf)*xn(INE,:))))*e;     % Active Fast, during positive errors
    xnp1(3,:) = xnp1(3,:)  +  0*(b0s      * (1  +  20*max(0,xn(IPE,:)     + (cs)*xn(INE,:))))*e;     % Active Slow, during positive errors
    xnp1(4,:) = xnp1(4,:)  +  0*( (b0us * (1  +  20*max(0,xn(IPE,:)       + (cus)*xn(INE,:))))*e);    % Active UltraSlow, during positive errors
    
else
    xnp1(1,:) = xnp1(1,:)  +  (b0uf    * (1   +  20*max(0,- (cuf)*xn(IPE,:)    -  xn(INE,:))))*e;   % Active UltraFast, during negative errors
    xnp1(2,:) = xnp1(2,:)  +  (b0f     * (1   +  20*max(0,- (cf)*xn(IPE,:)     -  xn(INE,:))))*e;    % Active Fast, during negative errors
    xnp1(3,:) = xnp1(3,:)  +  0*(b0s     * (1   +  20*max(0,- (cs)*xn(IPE,:)   -  xn(INE,:))))*e;    % Active Slow, during negative errors
    xnp1(4,:) = xnp1(4,:)  +  0*((b0us * (1   +  20*max(0,- (cus)*xn(IPE,:)    -  xn(INE,:))))*e);   % Active UltraSlow, during negative errors
end


% if e > 0
%     xnp1(1,:) = xnp1(1,:)  +  (b0uf     +  max(0,cuf*xn(IPE,:)        + (1-cuf)*xn(INE,:)))*e;    % Active UltraFast, during positive errors
%     xnp1(2,:) = xnp1(2,:)  +  (b0f      +  max(0,1*cf*xn(IPE,:)       + 1*(1-cf)*xn(INE,:)))*e;     % Active Fast, during positive errors
%     xnp1(3,:) = xnp1(3,:)  +  0*(b0s      +  max(0,1*cs*xn(IPE,:)       + 1*(1-cs)*xn(INE,:)))*e;     % Active Slow, during positive errors
%     xnp1(4,:) = xnp1(4,:)  +  0*( (b0us  +  max(0,cus*xn(IPE,:)     + (1-cus)*xn(INE,:)))*e);    % Active UltraSlow, during positive errors
%     
% else
%     xnp1(1,:) = xnp1(1,:)  +  (b0uf     +  max(0,- (1-cuf)*xn(IPE,:)    -  cuf*xn(INE,:)))*e;   % Active UltraFast, during negative errors
%     xnp1(2,:) = xnp1(2,:)  +  (b0f      +  max(0,- 1*(1-cf)*xn(IPE,:)   -  1*cf*xn(INE,:)))*e;    % Active Fast, during negative errors
%     xnp1(3,:) = xnp1(3,:)  +  0*(b0s      +  max(0,- 1*(1-cs)*xn(IPE,:)   -  1*cs*xn(INE,:)))*e;    % Active Slow, during negative errors
%     xnp1(4,:) = xnp1(4,:)  +  0*((b0us  +  max(0,- (1-cus)*xn(IPE,:)    -  cus*xn(INE,:)))*e);   % Active UltraSlow, during negative errors
% end


% if e > 0
%     xnp1(1,:) = xnp1(1,:)  +  (b0uf  +  cuf*xn(IPE,:)        + (1-cuf)*xn(INE,:))*e;    % Active UltraFast, during positive errors
%     xnp1(2,:) = xnp1(2,:)  +  (b0f   +  1*cf*xn(IPE,:)       + 1*(1-cf)*xn(INE,:))*e;     % Active Fast, during positive errors
%     xnp1(3,:) = xnp1(3,:)  +  (b0s   +  1*cs*xn(IPE,:)       + 1*(1-cs)*xn(INE,:))*e;     % Active Slow, during positive errors
%     xnp1(4,:) = xnp1(4,:)  + 0*( (b0us  +  cus*xn(IPE,:)     + (1-cus)*xn(INE,:))*e);    % Active UltraSlow, during positive errors
%     
% else
%     xnp1(1,:) = xnp1(1,:)  +  (b0uf  - (1-cuf)*xn(IPE,:)     - cuf*xn(INE,:))*e;   % Active UltraFast, during negative errors
%     xnp1(2,:) = xnp1(2,:)  +  (b0f   - 1*(1-cf)*xn(IPE,:)    - 1*cf*xn(INE,:))*e;    % Active Fast, during negative errors
%     xnp1(3,:) = xnp1(3,:)  +  (b0s   - 1*(1-cs)*xn(IPE,:)    - 1*cs*xn(INE,:))*e;    % Active Slow, during negative errors
%     xnp1(4,:) = xnp1(4,:)  +  0*((b0us  - (1-cus)*xn(IPE,:)  - cus*xn(INE,:))*e);   % Active UltraSlow, during negative errors
% end

% if e > 0
%     xnp1(1,:) = xnp1(1,:)  +  (b0uf  +   cuf*xn(IPE,:) + bias_uf*xn(INE,:)) *e;    % Active UltraFast, during positive errors
%     xnp1(2,:) = xnp1(2,:)  +  (b0f   +   cf*xn(IPE,:)  + bias_f*xn(INE,:))  *e;     % Active Fast, during positive errors
%     xnp1(3,:) = xnp1(3,:)  +  (b0s   +   0*cs*xn(IPE,:)  + 0*bias_s*xn(INE,:))  *e;     % Active Slow, during positive errors
%     xnp1(4,:) = xnp1(4,:)  +  0*(b0us  +   cus*xn(IPE,:) + bias_us*xn(INE,:)) *e;    % Active UltraSlow, during positive errors
%     
% else
%     xnp1(1,:) = xnp1(1,:)  +  (b0uf  -  cuf *xn(INE,:) + bias_uf*xn(IPE,:))*e;   % Active UltraFast, during negative errors
%     xnp1(2,:) = xnp1(2,:)  +  (b0f   -  cf  *xn(INE,:) + bias_f*xn(IPE,:)) *e;    % Active Fast, during negative errors
%     xnp1(3,:) = xnp1(3,:)  +  (b0s   -  0*cs  *xn(INE,:) + 0*bias_s*xn(IPE,:)) *e;    % Active Slow, during negative errors
%     xnp1(4,:) = xnp1(4,:)  +  0*(b0us  -  cus *xn(INE,:) + bias_us*xn(IPE,:))*e;   % Active UltraSlow, during negative errors
% end

%% Add bias
% xnp1(1,:) = xnp1(1,:) + bias_uf*(xn(IPE,:)  +  xn(INE,:));
% xnp1(2,:) = xnp1(2,:) + bias_f*(xn(IPE,:)  +  xn(INE,:));
% xnp1(3,:) = xnp1(3,:) + bias_s*(xn(IPE,:)  +  xn(INE,:));
% xnp1(4,:) = xnp1(4,:) + bias_us*(xn(IPE,:) +  xn(INE,:));

% xnp1(1,:) = xnp1(1,:) + bias_uf*(xn(IPE,:) +  xn(INE,:))*em;
% xnp1(2,:) = xnp1(2,:) + bias_f*(xn(IPE,:)  +  xn(INE,:))*em;
% xnp1(3,:) = xnp1(3,:) + bias_s*(xn(IPE,:)  +  xn(INE,:))*em;
% xnp1(4,:) = xnp1(4,:) + 0*bias_us*(xn(IPE,:) +  xn(INE,:))*em;
% 
% %     if e > 0
%         xnp1(2,:) = xnp1(2,:)  +  cfast*(b0f + xn(4,:) + (cint)*xn(5,:))*e;    % Active Fast, during positive errors
%         xnp1(3,:) = xnp1(3,:)  +        (b0s + xn(4,:) + (cint)*xn(5,:))*e;    % Active Slow, during positive errors
%     else
%         xnp1(2,:) = xnp1(2,:)  +  cfast*(b0f - ((cint)*xn(4,:) + xn(5,:)))*e;  % Active Fast, during negative errors
%         xnp1(3,:) = xnp1(3,:)  +        (b0s - ((cint)*xn(4,:) + xn(5,:)))*e;  % Active Slow, during negative errors
%     end


%     assert(xnp1(4,:)>=0,'x_pt<0')
%     assert(xnp1(5,:)<=0,'x_nt>0')
%     if e > 0
%         xnp1(4,:) = min([xnp1(4,:) + bet*em,  maxPA]);
%     else
%         xnp1(5,:) = max([xnp1(5,:) - bet*em, -maxPA]);
%     end
% if e > 0
%     xnp1(IPE,:) = min([xnp1(IPE,:) + bet*mySigmoid(e, c1, c2),  maxPA]);
% elseif e < 0
%     xnp1(INE,:) = max([xnp1(INE,:) + bet*mySigmoid(e, c1, c2), -maxPA]);
% else
%     % No update
% end
% eSmooth = xn(7,:);
% eSmooth = e;

%% e-learning--------------------------------------------------------------
% if eSmooth > th
%     xnp1(IPE,:) = min([xnp1(IPE,:) + bet*mySigmoid(eSmooth, sigma1, sigma2),  maxPA]);
% elseif eSmooth < -th
%     xnp1(INE,:) = max([xnp1(INE,:) + bet*mySigmoid(eSmooth, sigma1, sigma2), -maxPA]);
% else
%     % No update
% end

%% z-e learning------------------------------------------------------------
% if z*e>0 %Same sign
%     if e > th
%         xnp1(IPE,:) = min([xnp1(IPE,:) + bet*abs(z)*e,  maxPA]);
%     elseif e < -th
%         xnp1(INE,:) = max([xnp1(INE,:) + bet*abs(z)*e, -maxPA]);
%     else
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
%         % No update
%     end
% end

%% c-e learning------------------------------------------------------------
if z==0
    dz = 1;  %Consistency Factor
else
    dz = sum(xnp1(1:3,:)) - sum(xn(1:3,:));
%     dz = xnp1(1,:)-xn(1,:);
end

% if z > th
%     xnp1(IPE,:) = min([xnp1(IPE,:) + bet*z*xn(7,:),  maxPA]);
% elseif z < -th
%     xnp1(INE,:) = max([xnp1(INE,:) + bet*z*xn(7,:), -maxPA]);
% else
%     % No update
% end

% de = e - xn(7,:);
% if e*z>0
% if e > th
%     xnp1(IPE,:) = min([xnp1(IPE,:) + z*bet*(mySigmoid(abs(e), sigma1, sigma2, 0)),  maxPA]);
% elseif e < -th
%     xnp1(INE,:) = max([xnp1(INE,:) + z*bet*(mySigmoid(abs(e), sigma1, sigma2, 0)), -maxPA]);
% else
%     % No update
% end
% end

% consistency = 1./(1 + [0 abs(dz)]);
% learning    = abs((ae*e).^2./(ae*e + be*consistency.^-1));
% if e > th
%     xnp1(IPE,:) = min([xnp1(IPE,:) + learning,  maxPA]);
% elseif e < -th
%     xnp1(INE,:) = max([xnp1(INE,:) - learning, -maxPA]);
% else
%     % No update
% end
% be*z/(be + abs(dz))
% if e > th
%     xnp1(IPE,:) = min([xnp1(IPE,:) + bet*(1/abs(dz))*e/(be*e + ae*1/abs(dz)),  maxPA]);
% elseif e < -th
%     xnp1(INE,:) = max([xnp1(INE,:) + bet*(1/abs(dz))*e/(be*abs(e) + ae*1/abs(dz)), -maxPA]);
% else
%     % No update
% end

% % Long-term learning
% if e > th
%     xnp1(IPE,:) = min([xnp1(IPE,:) + (bet +  mySigmoid(e, sigma1, sigma2, 0))*e,  maxPA]);
% elseif e < -th
%     xnp1(INE,:) = max([xnp1(INE,:) + (bet +  mySigmoid(e, sigma1, sigma2, 0))*e, -maxPA]);
% else
%     % No update
% end

% if eSmooth > th
%     xnp1(IPE,:) = min([xnp1(IPE,:) + bet*myGaussSym(eSmooth, 0, mu, sigma1, sigma2, 3),  maxPA]);
% elseif eSmooth < -th
%     xnp1(INE,:) = max([xnp1(INE,:) + bet*myGaussSym(eSmooth, 0, mu, sigma1, sigma2, 3), -maxPA]);
% else
%     % No update
% end

% xnp1(7,:) = ae*xn(7,:) + be*e; %Smoothed error
% xnp1(7,:) = z;% Perturbation estimate
% xnp1(7,:) = 0*(ae*xn(7,:) + be*mySigmoid(abs(dz), sigma1, sigma2, 1));% Consistency of Perturbation Tracker



% %Update learning traces (complete spiking model)-------------------------
% % Positive lt
% xnp1(ix.lp,:)  = al*xn(ix.lp,:);
% if xn(ix.ip,:) >= th % IP-neuron spikes 
%     xnp1(ix.lp,:) = xnp1(ix.lp,:) + bl; %Transmit the spike
%     xnp1(ix.ip,:) = 0; % Reset the current to 0
% else % No spike. Keep integrating
%     xnp1(ix.ip,:) = ae*xn(ix.ip,:) + be*e*(e>0); 
% end
% 
% % Negative lt
% xnp1(ix.ln,:)  = al*xn(ix.ln,:);
% if xn(ix.in,:) <= -th %IN-neuron spikes
%     xnp1(ix.ln,:) = xnp1(ix.ln,:) - bl; %Transmit the spike
%     xnp1(ix.in,:) = 0; % Reset the current to 0
% else %No negative spike. Keep integrating
%     xnp1(ix.in,:) = ae*xn(ix.in,:) + be*e*(e<0);
% end
% End complete spiking model----------------------------------------------

% Update learning traces (reduced learning model)-------------------------
% Positive trace
xnp1(ix.lp,:) = ae*xn(ix.lp,:) + be*mySigmoid(e, sigma1, sigma2, 0)*(e>0);
% Negative trace
xnp1(ix.ln,:) = ae*xn(ix.ln,:) + be*mySigmoid(e, sigma1, sigma2, 0)*(e<0);

% Long-term perturbation learning (Motor primitives update)
xnp1(ix.zp,:) = aet*xn(ix.zp,:);
xnp1(ix.zn,:) = aet*xn(ix.zn,:);

if e > 0
    xnp1(ix.zp,:) = min([xnp1(ix.zp,:) + (bet +  xn(ix.lp))*e,  maxPA]);
elseif e < 0
    xnp1(ix.zn,:) = max([xnp1(ix.zn,:) + (bet +  abs(xn(ix.ln)))*e, -maxPA]);
else
    % No update
end
end

























