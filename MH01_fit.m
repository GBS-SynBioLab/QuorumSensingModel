clear all

% Load the name.sbproj project, which includes the variable m1, a model object
sbioloadproject 'Quorum Sensing Tra System07'

configsetObj = getconfigset(m1);

% Set the error tolerances 
set(configsetObj.SolverOptions, 'AbsoluteTolerance', 1.0e-12);
set(configsetObj.SolverOptions, 'RelativeTolerance', 1.0e-10);

%create a vector of 24 evenly spaced values for [AHL]
ahl_values = logspace(-15,-4, 23);
a = transpose(ahl_values);

% Create a SimFunction object, where AHL is the input parameter to scan, and X is the observed species. Pass in an empty array [] as the last input argument to denote there are no dosed species.
simfunc = createSimFunction(m1,{'AHL'},{'GFP','Ptra','Ptra_activated','TraR','TraR:AHL'},[]);

% Simulate the model multiple times with different AHL values. Set the stop time to 1000.
sd = simfunc(a,1000);

% selects time data and 
[t, x] = select(sd, {'Type','species'});

%code selects the endpoint values for protein expresion data
N = length(x) ;
vector = zeros(N,1) ;
for i = 1:N
    vector(i) = x{i}(end,1) ;
end


%% plot curve

loglog(a, vector, 'kx');
ylabel('[Protein] (M)');
xlabel('[AHL] (10^x M)');
title('Dose Response of Tra Quorum Sensing System');
% xlim([10^-16 10^-3])
hold on

%% curve fitting
hill=@(q,a)q(1)+((q(2)-q(1))./(1+10.^((log10(q(3))-log10(a)).*q(4))));
b0 = [vector(end) vector(1), 1E-8, 1];
h = nlinfit(a,vector,hill, b0);
q = h;

qdata=logspace(-15, -4, 100000);
for i = 1:length(qdata)
    y(i) = q(1)+((q(2)-q(1))./(1+10.^((log10(q(3))-log10(qdata(i)))*q(4))));
end
loglog(qdata,y);
% xlim([1E-16 1E-3])
% ylim([1E-5 2E-5])

hold off

% Plot the simulation results to see how varying the level of AHL affects the level of Ga
sbioplot(sd);