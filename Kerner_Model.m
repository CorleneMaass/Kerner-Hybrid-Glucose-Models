% Author: Corlene Maass & Zenae de Villiers 
% Course: Biomathematics 378
% Kerner's Model 
%--------------------------------------------------

% Required:
% MATLAB
% Optimization Toolbox (lsqnonlin)
% swtest.m (included in repository)

%--------------------------------------------------


%% FSIGT DATA (Kerner et al.)

% Sampling times (min)
timeData = [0 2 4 6 8 10 12 14 16 19 22 27 32 42 52 ...
62 72 82 92 102 122 142 162 182]';

% Blood glucose concentrations (mg/dL)
Gdata = [92 350 287 251 240 216 211 205 196 192 172 ...
163 142 124 105 92 84 77 82 81 82 82 85 90]';

% Plasma insulin concentrations (μU/mL)
Idata = [11 26 130 85 51 49 45 41 35 30 30 27 30 22 ...
15 15 11 10 8 11 7 8 8 7]';

%--------------------------------------------------

function dYdt = bergman_glucose(t,Y,p,timeIns,insData,Ib,Gb)

    % State variables
    G = Y(1);     % Blood glucose concentration
    X = Y(2);     % Remote insulin action
    
    % Interpolated plasma insulin input
    I = interp1(timeIns,insData,t,'linear','extrap');
    
    % Parameters
    p1 = p(1);
    p2 = p(2);
    p3 = p(3);
    
    % Bergman minimal model
    dG = p1*(Gb - G) - X*G;
    dX = p3*(I - Ib) - p2*X;
    
    % Return state derivatives
    dYdt = [dG; dX];

end

%--------------------------------------------------

%% Basal values
Gb = Gdata(end);
Ib = Idata(end);

%% Initial parameter estimates
p = [0.0399 0.02 4e-5];

%% Initial conditions
G0 = 287;
X0 = 0;
Y0 = [G0; X0];

%% Solve model
[tSol,Y] = ode45(@(t,Y) ...
    bergman_glucose(t,Y,p,timeData,Idata,Ib,Gb), ...
    [0 182],Y0);

G = Y(:,1);
X = Y(:,2);

%% Glucose profile
figure
plot(tSol,G,'b','LineWidth',2)
hold on
scatter(timeData,Gdata,40,'filled')
yline(Gb,'k--')

xlabel('Time (min)')
ylabel('Glucose (mg/dL)')
title('Kerner Model Simulation')
legend('Model','Data','Basal')
grid on

%% Remote insulin action
figure
plot(tSol,X,'LineWidth',2)

xlabel('Time (min)')
ylabel('X(t)')
title('Remote Insulin Action')
grid on

%--------------------------------------------------

function err = objective_function(p,timeData,Gdata,Idata,Gb,Ib)

    % Initial conditions
    G0 = 287;
    X0 = 0;
    Y0 = [G0; X0];
    
    % Solve model at experimental sampling times
    [~,Y] = ode45(@(t,Y) ...
    bergman_glucose(t,Y,p,timeData,Idata,Ib,Gb), ...
    timeData,Y0);
    
    % Model prediction
    Gmodel = Y(:,1);
    
    % Residual vector
    err = Gmodel - Gdata;

end

%--------------------------------------------------

%% Parameter Estimation

% Initial parameter estimates from Kerner et al.
p0 = [0.0399 0.02 4e-5];

% Parameter bounds
lb = [1e-4 1e-4 1e-8];
ub = [1 1 1e-3];

% Estimate parameters using nonlinear least squares
[p_est,resnorm,residual,exitflag,output,lambda,J] = ...
    lsqnonlin( ...
    @(p) objective_function(p,timeData,Gdata,Idata,Gb,Ib), ...
    p0,lb,ub);

%% Parameter Uncertainty
n = length(Gdata);
k = length(p_est);

sigma2 = resnorm/(n-k);

Cov = sigma2 * pinv(full(J'*J));

SE = sqrt(diag(Cov));

fprintf('\nParameter Estimates (± SE)\n')
fprintf('--------------------------\n')
fprintf('p1    = %.6g ± %.6g\n',p_est(1),SE(1))
fprintf('p2    = %.6g ± %.6g\n',p_est(2),SE(2))
fprintf('p3    = %.6g ± %.6g\n',p_est(3),SE(3))



%% Compare Initial and Final Residual Norms

initialNorm = norm(objective_function(p0,timeData,Gdata,Idata,Gb,Ib));

finalNorm = norm(objective_function(p_est,timeData,Gdata,Idata,Gb,Ib));

fprintf('\nResidual Norms\n')
fprintf('Initial residual norm = %.4f\n',initialNorm)
fprintf('Final residual norm = %.4f\n',finalNorm)

%% Simulate Model Using Estimated Parameters

[tFit,Yfit] = ode45( ...
    @(t,Y) bergman_glucose(t,Y,p_est,timeData,Idata,Ib,Gb), ...
    timeData,Y0);

Gfit = Yfit(:,1);

%% Plot Fitted Glucose Profile

figure

plot(timeData,Gfit,'b','LineWidth',2)
hold on

scatter(timeData,Gdata,50,'filled')

yline(Gb,'k--')

xlabel('Time (min)')
ylabel('Glucose (mg/dL)')
title('Fitted Kerner Model')
legend('Model Fit','Experimental Data','Basal Glucose')
grid on

%% Clinical Indices

SI = p_est(3)/p_est(2);   % Insulin sensitivity
SG = p_est(1);            % Glucose effectiveness

fprintf('\nClinical Indices\n')
fprintf('Insulin sensitivity (SI) = %g\n',SI)
fprintf('Glucose effectiveness (SG) = %g\n',SG)

%% Goodness-of-Fit Metrics

% Root Mean Squared Error
RMSE = sqrt(mean((Gfit-Gdata).^2));

% Coefficient of determination
SSres = sum((Gdata-Gfit).^2);
SStot = sum((Gdata-mean(Gdata)).^2);

R2 = 1 - SSres/SStot;

fprintf('\nGoodness-of-Fit Metrics\n')
fprintf('RMSE = %.3f\n',RMSE)
fprintf('R^2 = %.4f\n',R2)

%% Interpolated Insulin Profile

tt = linspace(0,182,1000);

II = interp1(timeData,Idata,tt,'linear');

figure

plot(tt,II,'LineWidth',2)
hold on

scatter(timeData,Idata,'filled')

xlabel('Time (min)')
ylabel('Insulin (\muU/mL)')
title('Interpolated Insulin Input')
grid on

%% Information Criteria

k = 3;                  % Number of fitted parameters
n = length(Gdata);      % Number of observations

AIC = n*log(SSres/n) + 2*k;
BIC = n*log(SSres/n) + k*log(n);

fprintf('\nInformation Criteria\n')
fprintf('AIC = %.3f\n',AIC)
fprintf('BIC = %.3f\n',BIC)


%% Distribution and Residual Analysis 

residuals = Gdata - Gfit;


[h,p,W] = swtest(residuals)


figure
plot(timeData,residuals,'o-')
yline(0,'k--')
xlabel('Time')
ylabel('Residual')
title('Residuals vs Time (Kerner Model)')
grid on