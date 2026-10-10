% Author: Corlene Maass & Zenae de Villiers 
% Course: Biomathematics 378
% Hybrid Model 
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

function dYdt = hybrid_model(t,Y,p,Ib,Gb)

    % State variables
    G = Y(1);    % Blood glucose concentration
    X = Y(2);    % Remote insulin action
    I = Y(3);    % Plasma insulin concentration
    
    % Parameters
    p1    = p(1);   % Glucose effectiveness
    p2    = p(2);   % Decay of insulin action
    p3    = p(3);   % Generation of insulin action
    gamma = p(4);   % Insulin clearance rate
    delta = p(5);   % Glucose-stimulated insulin secretion
    
    % Glucose dynamics
    dG = p1*(Gb - G) - X*G;
    
    % Remote insulin action dynamics
    dX = p3*(I - Ib) - p2*X;
    
    % Insulin dynamics
    dI = -gamma*(I - Ib) + delta*(G - Gb)*I;
    
    % Return derivatives
    dYdt = [dG; dX; dI];

end

%--------------------------------------------------


%% Basal Values

Gb = Gdata(end);
Ib = Idata(end);

%% Initial Parameter Estimates

p0 = [0.0399 0.02 4e-5 0.2 0.001];

%% Initial Conditions

G0 = 287;
% Initial glucose concentration used in the Kerner approach
% (not the maximum measured glucose value)

X0 = 0;

I0 = Ib;

Y0 = [G0; X0; I0];

%% Simulate Hybrid Model

[tSol,Y] = ode45( ...
    @(t,Y) hybrid_model(t,Y,p0,Ib,Gb), ...
    [0 182],Y0);

%% Extract State Variables

G = Y(:,1);
X = Y(:,2);
I = Y(:,3);

%% Glucose Profile

figure

plot(tSol,G,'b','LineWidth',2)
hold on

scatter(timeData,Gdata,40,'filled')

yline(Gb,'k--')

xlabel('Time (min)')
ylabel('Glucose (mg/dL)')
title('Hybrid Glucose-Insulin Model')
legend('Model','Data','Basal')
grid on

%% Remote Insulin Action

figure

plot(tSol,X,'LineWidth',2)

xlabel('Time (min)')
ylabel('X(t)')
title('Remote Insulin Action')
grid on

%% Insulin Profile

figure

plot(tSol,I,'LineWidth',2)

xlabel('Time (min)')
ylabel('Insulin (\muU/mL)')
title('Model Insulin Response')
grid on

%--------------------------------------------------

function err = objective_hybrid(p,timeData,Gdata,Gb,Ib)

    % Initial conditions
    G0 = 287;
    X0 = 0;
    I0 = Ib;
    
    Y0 = [G0; X0; I0];
    
    % Solve model at experimental sampling times
    [~,Y] = ode45( ...
        @(t,Y) hybrid_model(t,Y,p,Ib,Gb), ...
        timeData,Y0);
    
    % Model glucose prediction
    Gmodel = Y(:,1);
    
    % Residual vector
    err = Gmodel - Gdata;

end

%--------------------------------------------------

%% Basal Values

Gb = Gdata(end);
Ib = Idata(end);

%% Initial Conditions

G0 = 287;
% Initial glucose concentration used in the Kerner approach
% (not the maximum measured glucose value)

X0 = 0;

I0 = Ib;

Y0 = [G0; X0; I0];

%% Parameter Estimation

% Initial parameter estimates
p0 = [0.0399 0.02 4e-5 0.2 0.001];

% Parameter bounds
lb = [1e-4 1e-4 1e-8 1e-4 1e-6];
ub = [1 1 1e-3 10 1];

% Estimate parameters using nonlinear least squares
p_est = lsqnonlin( ...
    @(p) objective_hybrid(p,timeData,Gdata,Gb,Ib), ...
    p0,lb,ub);

%% Stability Check

% Local stability of the physiological equilibrium requires
%  (Routh-Hurwitz condition)

A = p_est(1) + p_est(2) + p_est(4);

B = p_est(1)*p_est(2) + p_est(1)*p_est(4) + p_est(2)*p_est(4);

C = p_est(1)*p_est(2)*p_est(4) + p_est(3)*p_est(5)*Gb*Ib;

if A*B > C
    disp('Stability condition satisfied')
else
    disp('Stability condition violated')
end


%% Estimated Parameters

fprintf('\nEstimated Parameters\n')
fprintf('--------------------\n')

fprintf('p1     = %.6g\n',p_est(1))
fprintf('p2     = %.6g\n',p_est(2))
fprintf('p3     = %.6g\n',p_est(3))
fprintf('gamma  = %.6g\n',p_est(4))
fprintf('delta  = %.6g\n',p_est(5))

%% Residual Norm Comparison

initialNorm = norm( ...
    objective_hybrid(p0,timeData,Gdata,Gb,Ib));

finalNorm = norm( ...
    objective_hybrid(p_est,timeData,Gdata,Gb,Ib));

fprintf('\nResidual Norms\n')
fprintf('--------------\n')

fprintf('Initial residual norm = %.4f\n',initialNorm)
fprintf('Final residual norm   = %.4f\n',finalNorm)

%% Simulate Fitted Hybrid Model

[tFit,Yfit] = ode45( ...
    @(t,Y) hybrid_model(t,Y,p_est,Ib,Gb), ...
    timeData,Y0);

%% Extract Fitted States

Gfit = Yfit(:,1);
Xfit = Yfit(:,2);
Ifit = Yfit(:,3);

%% Fitted Glucose Profile

figure

plot(timeData,Gfit,'b','LineWidth',2)
hold on

scatter(timeData,Gdata,50,'filled')

yline(Gb,'k--')

xlabel('Time (min)')
ylabel('Glucose (mg/dL)')
title('Fitted Hybrid Model')
legend('Model Fit','Experimental Data','Basal Glucose')
grid on

%% Clinical Indices

SI = p_est(3)/p_est(2);
SG = p_est(1);

fprintf('\nClinical Indices\n')
fprintf('----------------\n')

fprintf('Insulin sensitivity (SI) = %g\n',SI)
fprintf('Glucose effectiveness (SG) = %g\n',SG)

%% Goodness-of-Fit Metrics

RMSE = sqrt(mean((Gfit - Gdata).^2));

SSres = sum((Gdata - Gfit).^2);
SStot = sum((Gdata - mean(Gdata)).^2);

R2 = 1 - SSres/SStot;

fprintf('\nGoodness-of-Fit Metrics\n')
fprintf('-----------------------\n')

fprintf('RMSE = %.3f\n',RMSE)
fprintf('R^2  = %.4f\n',R2)

%% Fitted Insulin Profile

figure

plot(tFit,Ifit,'LineWidth',2)
hold on

scatter(timeData,Idata,'filled')

xlabel('Time (min)')
ylabel('Insulin (\muU/mL)')
title('Hybrid Model Insulin Response')
legend('Model','Data')
grid on

%% Information Criteria

k = 5;                % Number of fitted parameters
n = length(Gdata);    % Number of observations

AIC_hybrid = n*log(SSres/n) + 2*k;
BIC_hybrid = n*log(SSres/n) + k*log(n);

fprintf('\nInformation Criteria\n')
fprintf('--------------------\n')

fprintf('Hybrid AIC = %.3f\n',AIC_hybrid)
fprintf('Hybrid BIC = %.3f\n',BIC_hybrid)


%% Distribution and Residual Analysis 

residuals = Gdata - Gfit;


[h,p,W] = swtest(residuals)


figure
plot(timeData,residuals,'o-')
yline(0,'k--')
xlabel('Time')
ylabel('Residual')
title('Residuals vs Time (Hybrid Model)')
grid on