%% models parameters
r =1; % signifies original model
m = []; % matrix that stores model parameters from simulation
q = fit(0,0,r);
m = [m; q];

r =0; % signifies model variant

% ref value = 1.0
pN = 'ktranslation GFP_mRNA';

pV = 10;
q = fit(pN,pV,r);
m = [m; q];

pV = 0.1;
q = fit(pN,pV,r);
m = [m; q];


% ref value = 0.208
pN = 'ktranscription Ptra_basal';

pV = 2.08;
q = fit(pN,pV,r);
m = [m; q];

pV = 0.0208;
q = fit(pN,pV,r);
m = [m; q];


% ref value = 10000
pN = 'kon Ptra/TraR:AHL';

pV = 100000;
q = fit(pN,pV,r);
m = [m; q];

pV = 1000;
q = fit(pN,pV,r);
m = [m; q];


% ref value = 5.51
pN = 'ktranscription Ptra_activated';

pV = 11.02;
q = fit(pN,pV,r);
m = [m; q];

pV = 2.255;
q = fit(pN,pV,r);
m = [m; q];


% ref value = 1.0
pN = 'ktranscription TraR';

pV = 5;
q = fit(pN,pV,r);
m = [m; q];

pV = 0.2;
q = fit(pN,pV,r);
m = [m; q];

r=2;
%ref value = 1.66E-8
pN = 'Ptra';

pV = 1.66E-7;
q = fit(pN,pV,r);
m = [m; q];

pV = 1.66E-9;
q = fit(pN,pV,r);
m = [m; q];

pN = 'decoy';

pV = 1.66E-7;
q = fit(pN,pV,r);
m = [m; q];

pV = 1.66E-9;
q = fit(pN,pV,r);
m = [m; q];



%% metrics analysis

fold_mat = [];
fold_ch_m = m(1,1)/m(1,2);

scal_mat = [];
scal_m = (m(1,2)/2)*(1+(m(1,1)/m(1,2)));

resp_mat = [];
resp_m = ((m(1,4)*(m(1,1)-m(1,2)))/(m(1,3)*((1+1))^2));

for i = 2:length(m)

fold_ch = m(i,1)/m(i,2);
per_fold_cha = ((fold_ch - fold_ch_m)/fold_ch_m)*100;
fold_mat = [fold_mat per_fold_cha];

scal = (m(i,2)/2)*(1+(m(i,1)/m(i,2)));
per_scal = ((scal - scal_m)/scal_m)*100;
scal_mat = [scal_mat per_scal];

resp =((m(i,4)*(m(i,1)-m(i,2)))/(m(i,3)*((1+1))^2));
per_resp = ((resp-resp_m)/resp_m)*100;
resp_mat = [resp_mat per_resp];

end
colour = rgb('RosyBrown');
fig = figure;
xnames = categorical({'ktranslation GFP mRNA (x10)','ktranslation GFP mRNA (:10)','ktranscription Ptra basal (x10)', 'ktranscription Ptra basal (:10)', 'kon Pta/TraR:AHL(x10)', 'kon Ptra/TraR:AHL(:10)', 'ktranscription Ptra activated (x2)', 'ktranscription Ptra activated (:2)', 'ktranscription TraR (x5)', 'ktranscription TraR (:5)', 'Ptra copies (x10)', 'Ptra copies (:10)', 'decoy copies (x10)', 'decoy copies (:10)'});
subplot(2,2,1);
barh(xnames, fold_mat, 'FaceColor',colour);
line([5 5], ylim, 'Color', 'black', 'LineStyle', '--');
line([-5 -5], ylim, 'Color', 'black', 'LineStyle', '--');
title('Fold change');
xlabel('Percentage change (%)');
hold on

subplot(2,2,2);
barh(xnames, scal_mat, 'FaceColor',colour);
line([5 5], ylim, 'Color', 'black', 'LineStyle', '--');
line([-5 -5], ylim, 'Color', 'black', 'LineStyle', '--');
title('Scaling');
xlabel('Percentage change (%)');
hold on

subplot(2,2,3);
barh(xnames, resp_mat, 'FaceColor',colour);
line([5 5], ylim, 'Color', 'black', 'LineStyle', '--');
line([-5 -5], ylim, 'Color', 'black', 'LineStyle', '--');
title('Responsiveness');
xlabel('Percentage change (%)');
hold off

%% simulation and fitting
function q=fit(pN,pV,r)

    % Load the name.sbproj project, which includes the variable m1, a model object
    sbioloadproject 'Quorum Sensing Tra Systemdecoy'

    if r ==0
        variantObj = addvariant(m1, 'v2');
        addcontent(variantObj, {'parameter', pN, 'Value', pV});
        variantObj.Active = true;
    elseif r==2
        variantObj = addvariant(m1, 'v3');
        addcontent(variantObj, {'species', pN, 'InitialAmount', pV});
        variantObj.Active = true;
    end


    %create a vector of 24 evenly spaced values for [AHL]
    ahl_values = logspace(-15,-4, 35);
    c = transpose(ahl_values);

    % Create a SimFunction object, where AHL is the input parameter to scan, and 
    % X is the observed species. Pass in an empty array [] as the last input argument to denote there are no dosed species.
    if r == 1
        simfunc = createSimFunction(m1,{'AHL'},{'GFP','Ptra','Ptra_activated','TraR','TraR:AHL'},[]);
    elseif r ==0 || r==2
        simfunc = createSimFunction(m1,{'AHL'},{'GFP','Ptra','Ptra_activated','TraR','TraR:AHL'},[], variantObj);
    end
    % Simulate the model multiple times with different AHL values. Set the stop time to 1000.
    sd = simfunc(c,1000);

    % selects time data and outputs
    [t, x] = select(sd, {'Type','species'});

    %code selects the endpoint values for protein expresion data
    N = length(x) ;
    vector = zeros(N,1) ;
    for i = 1:N
        vector(i) = x{i}(end,1) ; % x{i} = simulation number i.e AHL conc.; (end,1) = species (GFP = 1)
    end


    %% plot curve

    %semilogx(c, vector, 'kx');
    %ylabel('[GFP] (M)');
    %xlabel('[AHL] (10^x M)');
    %title('Dose Response of Tra Quorum Sensing System');
    %xlim([10^-16 10^-3])
    %hold on

    %% curve fitting
    hill=@(q,c)q(2)+((q(1)-q(2))./(1+10.^((log10(q(3))-log10(c)).*q(4))));
    b0 = [vector(end) vector(1), 1E-8, 1];
    h = nlinfit(c,vector,hill, b0);
    q = h;

    %qdata=logspace(-15, -4, 10000);
    %for i = 1:length(qdata)
    %    y(i) = q(2)+((q(1)-q(2))./(1+10.^((log10(q(3))-log10(qdata(i)))*q(4))));
    %end
    %semilogx(qdata,y);
    %hold off

end