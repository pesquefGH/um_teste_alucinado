% LNCC SURGE Model 26  Chikungunya- Validation 2
% Tasks of this script (for a given BR state):
% 1) Estimate a surge model (logistic model) from the average surge behavior across the
% years, by means of a nonlinear estimator. We can consider the surge model
% as a typical (expected) surge in a given year.

% 2) From the typical surge we estimate the parameters of an incidence model, 
% from which we generate a modeled surge. 

% 3) Obtain a set of gains across the years that, when applied to the surge model, 
% match its scale to that of a given observed surge 
% (time-synchronized with the surge model). A forecast realization n is
% then calculated  as: forecast (n) = g(n) * modeled_surge. We generate 
% 10000 values of g(n) drawn from a log normal distribution. 

% 4) Provide a 'Validation 2' Forecast based on the model of the average surge

SUNDAYS = readtable('sunday_dates.csv');  % reads CSV with Sunday dates
S_dates=SUNDAYS{:,1};   % Sunday dates format YYYY-MM-DD

M = readtable(['IMDC2026_AggregatedData-Chikungunya_',UF,'.csv']); % reads aggregated
% data from the chosen state, from EW 1 of 2014 to EW 52 of 2025

SS=52;  % seasonality of 52 Epidemic Weeks (EW)

% Validation 2 

ind_EW25_2023=441+52; % time index of EW 25 of 2023 


dcases=M{:,3};   % time-series of the surges (incidence curves) across the EWs 

dcases_orig=dcases;  % stores the original time-series up to EW 52 of 2025

dcases=dcases(1:ind_EW25_2023);  % crops time-series at EW 25 of 2023 (including)

th=7*abs(median(dcases)); % threshold for discarding surges below it
% Heuristically chosen


L=26; % window length for the SSA filter
nsv=6; % number of selected eigenvalues (ordered)
[dcasesf]=ssa_modPE(dcases,L,nsv); % filtered time series


dcasesfs=dcasesf(41:ind_EW25_2023);  % selects filtered data from 
% EW 41/2014 to EW 25/2023

% From the time-series 'dcasesfs' of the Chikungunya surges, we obtain an average surge curve.
% The hypothesis is that the surges have a seasonality of 52 Epidemic Week
% (EW), and that a properly scaled average surge is sufficient to represent
% an observed surge in a given year.

DC=buffer(dcasesfs,SS); % organizes vector dcases in a matrix with PP rows
% each column of DC has cases for 52 EWs. 

[DCalign, ind_max]=alignDC_CG(DC,th); % aligns surges so that they all peak
% in the middle of the window between EW 41 year to EW 40 subsequent year

typ_DC=mean(DCalign',"omitnan")'; % typical surge waveform in 52 EWs:
% from EW 41 to EW 40 (subsequent year)  

% Surge Model Estimation (from typ_DC)

% Initial parameters of the logistic model
L=120000;  
k=0.3;
n0=26;  % time-shift
P=[L,k,n0];  % initial parameter vector
n=0:length(typ_DC)-1; n=n(:);  % time basis

% Surge Model Estimation (nonlinear) 
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective');
fun = @(P,n) (P(1).*P(2).*exp(P(2).*(n-(P(3)))))./(1+exp(P(2).*(n-(P(3)))).^2); % surge model
P_est = lsqcurvefit(fun,P,n,typ_DC,[100;0.1;24],[370000;0.5;28],options); % find estimated model parameters 
% P_est is the optimal parameter vector

P_est(3)=round(P_est(3)); % rounds estimated n0, since it is supposed to be an integer 

Model_surge=fun(P_est,n);  % Synthesize a surge from estimated model 

% Note: here, we are estimating just one model for the average surge, which
% we supposed would be enough to represent all other surges, up to a scale factor.   
% Thus, for the purpose of surge forecasting, we only need a surge model 
% and some statistics related to the set of scale gains. 
% We could also estimate several surge models, one for each year. This way,
% we could also carry out a surge forecast using statistics
% for the surge model parameters (L, k, and n0).

% Now, we calculate a set of gains g so that
% an observed surge is given by g*Model_surge

ns=size(DCalign,2); % number of observed surges (number of columns of DCalign)

g=zeros(ns,1);  % initializes with zeros vector to store the set of gains
x=Model_surge;   % surge template: attributes Model_surge to x (for notation clarity)  
for kk=1:ns
    
    a=DCalign(:,kk);   % observed surge at column kk
    cross_corr_coef=corrcoef(x,a); cross_corr_coef=cross_corr_coef(2,1);
    if isnan(cross_corr_coef)
        g(kk)=NaN; 
    else
        g(kk)=a'*x./(x'*x);   % calculates amplitude gain for each observed surge at column kk
    end
end

g(isnan(g))=[];  % removes uncertain gains (when cross correlation between
% observed and modeled surge is too low)

% A simple predictor to forecast the cases in Validation 1 is to generate
% a set of values of g for the next season based on the mean and variance
% of the previous values of g. For simplicity a log normal distribution can
% be assumed for generating g.

[param]=lognfit(g);

mg=param(1);  % mean of value of g (from 2014 to 2023)
sigma=param(2); % standard deviation of value of g (from 2014 to 2023)

MC=10000;  % number of montecarlo runs for forecast surges in 2024

g_MC=lognrnd(mg,sigma,MC,1);  % set of randomly generated gains

forecast_cases_v1=zeros(length(Model_surge),MC);  % matrix to store
% realizations of the forecast surges 
% from EW 41 of 2023 to EW 40 of 2024

for kk=1:MC
    forecast_cases_v1(:,kk)=g_MC(kk)*Model_surge; % set of MC realizations
    % of the surge forecast
end

% Computes the median forecast and related prediction intervals
set_prctile=[2.5 5 10 25 50 75 90 95 97.5]; % 2.5 to 97.5% percentiles
PP=prctile(forecast_cases_v1',set_prctile); % calculates the percentiles and stores in PP
lower_95 = PP(1,:)';  % 2.5% percentile 
lower_90 = PP(2,:)';  % 5% percentile
lower_80 = PP(3,:)';  % 10% percentile
lower_50 = PP(4,:)';  % 25% percentile
pred = PP(5,:)';  % 50% percentile - median prediction
upper_50 = PP(6,:)'; % 75% percentile
upper_80 = PP(7,:)'; % 90% percentile 
upper_90 = PP(8,:)'; % 95% percentile 
upper_95 = PP(9,:)'; % 97.5% percentile 

indf_ini=457+52; % time index of the EW 41 2023
indf_end=508+52; % time index of the EW 40 2024

date=S_dates(719:770);  % Sunday days for validation 2 

cases=dcases_orig(indf_ini:indf_end);  % known cases 

state_code=M{1,2}*ones(size(pred));  % state code vector

T = table(date,pred,lower_50,upper_50,lower_80,upper_80,lower_90,upper_90,lower_95,upper_95);

writetable(T,['..\spreadsheets\v2_SModel26_CG_',num2str(state_code(1)),'.csv'],'Delimiter',',')

% Generates plots and save as PDF files
max75=max(upper_50);
max90=max(upper_80);
max95=max(upper_90);
max97p5=max(upper_95);

EW_index=indf_ini:indf_end;

figure
subplot(221)
plot(EW_index,cases,'linewidth',2); % plot known sequence of cases 
% from EW 41 2023 up to EW 40 of 2024
hold on;
plot(EW_index,pred,'r','linewidth',2)
plot(EW_index,lower_50,'k','linewidth',2)
plot(EW_index,upper_50,'g','linewidth',2)
legend('observed','median forecast','LB50','UB50')  
xlabel('Time (EW index)')
ylabel('Number of Cases')
title(['Forecast 50% -  ',UF])
axis([EW_index(1) EW_index(end) 0 2*max75])

subplot(222)
plot(EW_index,cases,'linewidth',2); % plot known sequence of cases
% from EW 41 2023 up to EW 40 of 2024
hold on;
plot(EW_index,pred,'r','linewidth',2)
plot(EW_index,lower_80,'k','linewidth',2)
plot(EW_index,upper_80,'g','linewidth',2)
legend('observed','median forecast','LB80','UB80')  
xlabel('Time (EW index)')
ylabel('Number of Cases')
title(['Forecast 80% -  ',UF])
axis([EW_index(1) EW_index(end) 0 2*max90])

subplot(223)
plot(EW_index,cases,'linewidth',2); % plot known sequence of cases
% from EW 41 2023 up to EW 40 of 2024
hold on;
plot(EW_index,pred,'r','linewidth',2)
plot(EW_index,lower_90,'k','linewidth',2)
plot(EW_index,upper_90,'g','linewidth',2)
legend('observed','median forecast','LB90','UB90')  
xlabel('Time (EW index)')
ylabel('Number of Cases')
title(['Forecast 90% -  ',UF])
axis([EW_index(1) EW_index(end) 0 2*max95])

subplot(224)
plot(EW_index,cases,'linewidth',2); % plot known sequence of cases
% from EW 41 2023 up to EW 40 of 2024
hold on;
plot(EW_index,pred,'r','linewidth',2)
plot(EW_index,lower_95,'k','linewidth',2)
plot(EW_index,upper_95,'g','linewidth',2)
legend('observed','median forecast','LB95','UB95')  
xlabel('Time (EW index)')
ylabel('Number of Cases')
title(['Forecast 95% -  ',UF])
axis([EW_index(1) EW_index(end) 0 2*max97p5])

print(['..\plots\v2_SModel26_CG_',UF],'-dpdf')

close all



