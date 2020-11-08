%% 4. Plant model identification
%% 4.0 Initialization
clear;

RUN = 21;

folder = sprintf("./DATA/4_PlantModelIdentification%02d/",RUN);
% Create directory
mkdir(folder);
mkdir(folder+"polesZerosTrue");
% Save directory path
save(folder+"RUN.mat",...
    'RUN','folder');

%% 4.1 Generate input signals ------------------------------------------ %
tic
clear

RUN = 21;
T = 300;                        % (s) Duration of the stimulus
Nu = 10;                        % number of input signals (has to be even)
uB_RBS_span = [0.01 0.05];      % (Hz) Define PRBS bandwidth 
uf_square_span = [0.01 0.05];   % (Hz) Define square wave frequency

% Load directory path
load(sprintf("./DATA/4_PlantModelIdentification%02d/"+"RUN.mat",RUN),...
    'folder');

fs = 50;
t = (0:1/fs:T)';

% Generate open loop input signals
uB_PRBS = (uB_RBS_span(1):...
    (uB_RBS_span(2)-uB_RBS_span(1))/(Nu/2-1):uB_RBS_span(2))'; 
uf_square = (uf_square_span(1):...
    (uf_square_span(2)-uf_square_span(1))/(Nu/2-1):uf_square_span(2))'; 

u = cell(1,Nu);
for i = 1:Nu
    if i<= Nu/2
        u{i} = idinput(length(t),'prbs',[0 uB_PRBS(i)]); 
    else
        u{i} = square(2*pi*uf_square(i-Nu/2)*t);
    end
end

% Save input signals
save(folder+"inputSignals.mat",...
    'T','fs','t','uB_PRBS','uf_square','u','Nu');

fprintf('Input signals obtention duration: %03f\n',toc)
%% 4.2 Simulate systems with selected input signals
tic
clear;

RUN = 21;

% Load directory path
load(sprintf("./DATA/4_PlantModelIdentification%02d/"+"RUN.mat",RUN),...
    'folder');
% Load input signals
load(folder+"inputSignals.mat",...
    'T','t','fs','u','Nu');

% Simulate system response 
y = cell(1,Nu);

for i = 1:Nu
    y{i} = simOpenLoop(t,u{i},1/fs);
end

% Save raw responses
save(folder+"inputSignalsResponse.mat",...
    'y');

fprintf('Response signals obtention duration: %03f\n',toc)
%% 4.3 Preprocessing --------------------------------------------------- %
tic
clear;

RUN = 21;
lambda = 0.6;
Tdel = 10;
Afilter = [1 -lambda];
Bfilter = (1-lambda)*[1 -1];

% Load directory path
load(sprintf("./DATA/4_PlantModelIdentification%02d/"+"RUN.mat",RUN),...
    'folder');
% Load number of input signals
load(folder+"inputSignals.mat",...
    't','u','fs','Nu');
% Load raw responses
load(folder+"inputSignalsResponse.mat",...
    'y');

ndel = round(Tdel*fs);

yf = cell(1,Nu);
uf = cell(1,Nu);

for i = 1:Nu
    tf = t(ndel+1:end);
    
    uf{i} = u{i}(ndel+1:end); 
    uf{i} = detrend(uf{i});
    
    yf{i} = y{i}(ndel+1:end);
    yf{i} = filter(Bfilter,Afilter,yf{i});
end

% Save filtered responses
save(folder+"inputSignalsResponseProcessed.mat",...
    'tf','yf','uf','lambda');

fprintf('Preprocessing duration: %03f\n',toc)
%% 4.4 Identification
tic
clear;

RUN = 21;

nASpan = [3,6];

% Load directory path
load(sprintf("./DATA/4_PlantModelIdentification%02d/"+"RUN.mat",RUN),...
    'folder');
% Load number of input signals
load(folder+"inputSignals.mat",...
    'Nu');
% Load filtered responses
load(folder+"inputSignalsResponseProcessed.mat",...
        'yf','uf');

Nidentification = 0;
for nA = nASpan(1):nASpan(2)
    for nB = 1:nA-1
        for nK = 1:nA-nB+1
            Nidentification = Nidentification+1;
        end
    end
end

% Find armax model for each model order, for each 
M = cell(Nu,Nidentification);
nABCK = zeros(Nidentification,4);

for i = 1:Nu
    count = 0;
    
    for nA = nASpan(1):nASpan(2)
        for nB = 1:nA-1
            for nK = 1:nA-nB+1
                nC = nA;
                count = count+1;
                M{i,count} = armax([yf{i} uf{i}],[nA nB nC nK]);
                nABCK(count,:) = [nA nB nC nK];
            end
        end
    end
end

% Save filtered responses
save(folder+"identificationCellArray.mat",...
    'M','nABCK','Nidentification','nASpan');

fprintf('Identification duration: %03f\n',toc)

%% 4.5 Model validation - Average fits calculation              
tic
clear

RUN = 21;

% Load directory path
load(sprintf("./DATA/4_PlantModelIdentification%02d/"+"RUN.mat",RUN),...
    'folder');
% Load number of input signals
load(folder+"inputSignals.mat",...
    'Nu');
% Load filtered responses
load(folder+"inputSignalsResponseProcessed.mat",...
    'yf','uf');
% Load filtered responses
load(folder+"identificationModels.mat",...
    'M','nABCK','Nidentification','nASpan');

% Get fitness of each model
fitM = zeros(Nu,Nidentification);

for i = 1:Nu
    for count = 1:Nidentification
        for k = 1:Nu
            if k ~= i
                [~,fit] = compare([yf{i} uf{i}],M{i,count});
                fitM(i,count) = fitM(i,count)+fit/(Nu-1);
            end
        end  
    end
end

ordersNumber = nASpan(2)-nASpan(1)+1;
Nbest = 10;

fitMcopy = fitM;
bestFitFit = zeros(ordersNumber,Nbest);
bestFitCount = zeros(ordersNumber,Nbest);
bestFitU = zeros(ordersNumber,Nbest);
count = 0;

for nA = nASpan(1):nASpan(2)
    for b = 1:(min(Nbest,sum(nABCK(:,1)==nA)))
        idx = find(nABCK(:,1)==nA);
        
        bestFitFit(nA-2,b) = mean(max(max(fitMcopy(:,idx(1):idx(end)))));
        [aux1,aux2] = find(fitMcopy(:,idx(1):idx(end))==bestFitFit(nA-2,b));
        bestFitU(nA-2,b) = aux1(1);
        bestFitCount(nA-2,b) = aux2(1);
        bestFitCount(nA-2,b) = bestFitCount(nA-2,b)+idx(1)-1;
        fitMcopy(bestFitU(nA-2,b),bestFitCount(nA-2,b)) = 0;
    end
end

% Save identification fitness data
save(folder+"identificationFitness.mat",...
    'fitM','bestFitValue','bestFitCount','bestFitU');

fprintf('Best fit values\n')
disp(bestFitFit)
fprintf('Best fit input waves\n')
disp(bestFitU)
fprintf('Best fit orders combination numbers\n')
disp(bestFitCount)

fprintf('Fits ranking duration: %03f\n',toc)

%% 4.5.1 Check results
tic
clear

RUN = 21;

% Load directory path
load(sprintf("./DATA/4_PlantModelIdentification%02d/"+"RUN.mat",RUN),...
    'folder');
% Load identified models
load(folder+"identificationCellArray.mat",...
    'M','nABCK');
% Load best fits
load(folder+"identificationFitness.mat",...
    'bestFitCount','bestFitU','bestFitFit');
% Load input signals information
load(folder+"inputSignals.mat",...
    'fs');

fprintf('Check results loading duration: %03f\n', toc)

%% Model choice
chooseIdentification = true;
nA = 3;
b = 2;  % b-th best

% only works for nK=1
if chooseIdentification
    load("./DATA/bodeData.mat",...
        'trueMag','truePhase','w');
    
    [den,num] = polydata(M{bestFitU(nA-2,b),bestFitCount(nA-2,b)});
    
    figure()
    num1 = [num(3:4) 0 0];
    zplane(num1,den);
    saveFigAsPDF(gcf,sprintf(...
        "./DATA/4_PlantModelIdentification%01d/NA%01d_ZPlane",RUN,nA))
    
    den = conv(den,[1 -1]); % Add integrator
    
    [mag,phase,wout] = bode(tf(num,den,1/fs),w);
    figure();
    subplot(211);
    semilogx(w,20*log10(trueMag));
    hold on;
    semilogx(w,20*log10(squeeze(mag)));
    subplot(212);
    semilogx(w,truePhase);
    hold on;
    semilogx(w,squeeze(phase));
    hold off;
    hold off;
    legend('Plant','Identified system','location','best')
    saveFigAsPDF(gcf,sprintf(...
        "./DATA/4_PlantModelIdentification%01d/NA%01d_Bode",RUN,nA))
end

%% 4.6 Postprocessing
choosenIdentification = true;
if choosenIdentification
    clear;
    RUN = 21;
    % Load directory path
    load(sprintf("./DATA/4_PlantModelIdentification%02d/"+"RUN.mat",RUN),...
        'folder');
    % Load identification
    load(folder+"identification.mat",...
        'MId','fsId','TsId','orderId','fitId','polesId_','zerosId_',...
        'denId_', 'numId_');
    % Add integrator
    numId = numId_;
    denId = conv(denId_,[1 -1]);
    [AId,BId,CId,DId] = tf2ss(numId,denId);
    
    polesId = pole(tf(numId,denId,TsId));
    zerosId = zero(tf(numId,denId,TsId));

    save(folder+"identification.mat",...
        'MId','fsId', 'TsId','orderId','fitId','polesId_','zerosId_',...
        'denId_', 'numId_','denId','numId','AId','BId','CId','DId',...
        'polesId','zerosId');
end

%% Batota
clear

RUN = 21;

% Load directory path
load(sprintf("./DATA/4_PlantModelIdentification%02d/"+"RUN.mat",RUN),...
    'folder');
load("./model/barrassmodel.mat",...
    'Atrue','Btrue','Ctrue','Dtrue');
load(folder+"inputSignals.mat",...
    'fs');

SysC = ss(Atrue,Btrue,Ctrue,Dtrue);
SysD = c2d(SysC,1/fs);
polesTrue = pole(SysD);
zerosTrue = zero(SysD);
[numTrue,denTrue] = tfdata(SysD,'v');

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
set(gca,'FontSize',35);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';

p = zplane(numTrue,denTrue);
title("True poles and zeros for "+sprintf("$T_s = %g$",1/fs)+" [s]",...
    'Interpreter','latex');
saveas(gcf,folder+sprintf("polesZerosTrue/polesZerosTrue.fig"));
hold off;

% Save true poles and zeros
save(folder+"polesZerosTrue/polesZerosTrue.mat",...
    'polesTrue','zerosTrue','numTrue','denTrue');