%% 4. Plant model identification
%% 4.0 Initialization
clear;
% ---------------------------------------------------------------------- %
% --------------------------- SET PARAMETERS --------------------------- %
RUN = 3;
% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %
folder = sprintf("./DATA/4_PlantModelIdentification%02d/",RUN);
% Create directory
mkdir(folder);
mkdir(folder+"filterTuning");
% Save directory path
save(folder+"RUN.mat",...
    'RUN','folder');
%% 4.1 Generate input signals ------------------------------------------ %
clear;
% ---------------------------------------------------------------------- %
% --------------------------- SET PARAMETERS --------------------------- %
RUN = 3;
T = 120; %(s) Duration of the stimulus
Nfs = 1; % number of frequencies to test
Nu = 10; % number of input signals (has to be even)
uB_RBS_span = [0.05 0.2]; % (Hz) % Define PRBS bandwidth 
uf_square_span = [0.1 0.5]; % (Hz) % Define square wave frequency
% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %
% Load directory path
load(sprintf("./DATA/4_PlantModelIdentification%02d/"+"RUN.mat",RUN),...
    'folder');
% Linear
% fs = 500:-:50 %(Hz)
% Logarithmic vector of sampling frequencies to test
%fs = 10.^(log10(1e3):(log10(1e1)-log10(1e3))/(Nfs-1):log10(1e1)); %(Hz)
fs = 100;
Ts = 1./fs; % (s) vector of sampling interval frequencies
t = cell(Nfs,1);
 % Generate time spans for different sampling frequencies
for i = 1:Nfs
    t{i,1} = (0:1/fs(i):T)';
end
% Generate open loop input signals
tic;
uB_PRBS = (uB_RBS_span(1):...
    (uB_RBS_span(2)-uB_RBS_span(1))/(Nu/2-1):uB_RBS_span(2))'; 
uf_square = (uf_square_span(1):...
    (uf_square_span(2)-uf_square_span(1))/(Nu/2-1):uf_square_span(2))'; 
u = cell(Nu,Nfs);
for j = 1:Nfs
    for i = 1:Nu
        if i<= Nu/2
            u{i,j} = idinput(length(t{j,1}),'prbs',[0 uB_PRBS(i)]); 
        else
            u{i,j} = square(2*pi*uf_square(i-Nu/2)*t{j,1});
        end
    end
end
toc;
% Save input signals
save(folder+"inputSignals.mat",...
    'T','Nfs','fs','Ts','t','Nu','uB_PRBS','uf_square','u');
%% 4.2 Simulate systems with selected input signals
clear;
% ---------------------------------------------------------------------- %
% --------------------------- SET PARAMETERS --------------------------- %
RUN = 3;
% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %
% Load directory path
load(sprintf("./DATA/4_PlantModelIdentification%02d/"+"RUN.mat",RUN),...
    'folder');
% Load input signals
load(folder+"inputSignals.mat",...
    'T','Nfs','Ts','t','Nu','u');
% Simulate system response 
tic;
y = cell(Nu,Nfs);
for j = 1:Nfs
    for i = 1:Nu
        y{i,j} = simOpenLoop(t{j,1},u{i,j},Ts(j));
    end
end
toc;
% Save raw responses
save(folder+"inputSignalsResponse.mat",...
    'y');
%% 4.3 Preprocessing --------------------------------------------------- %
clear;
% ---------------------------------------------------------------------- %
% --------------------------- SET PARAMETERS --------------------------- %
RUN = 3;
isTuning = false;
% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %
% Load directory path
load(sprintf("./DATA/4_PlantModelIdentification%02d/"+"RUN.mat",RUN),...
    'folder');
% Load number of input signals
load(folder+"inputSignals.mat",...
    'Nfs','Nu','t','u','fs');
% Load raw responses
load(folder+"inputSignalsResponse.mat",...
    'y');
if isTuning
    % Filter tuning tune lambda
    % Select input signal to tune lambda
    i = 6; j = 1;
    lambda = 0.7:0.05:0.95;
    lambda = [0 lambda]; % Unfiltered
    yf = cell(length(lambda),1);
    for k = 1:length(lambda)
        Afilter = [1 -lambda(k)];
        Bfilter = (1-lambda(k))*[1 -1];
        yf{k,1} = filter(Bfilter,Afilter,y{i,j});
    end
    yl = [inf inf]; % conatnt y axis limits
    for k = 0:length(lambda)
        figure('units','normalized','outerposition',[0 0 1 1]);
        hold on;
        set(gca,'FontSize',35);
        yyaxis left;
        plot(t{j,1}(round(length(t{j,1})/14):2*round(length(t{j,1})/14)),...
            u{i,j}(round(length(t{j,1})/14):2*round(length(t{j,1})/14)),...
            'Linewidth',4);
        if k
            title("$\lambda = $ "+sprintf("%g",lambda(k))+...
            "$\;|\;i =$ "+sprintf("%02d",i)+...
            "$\;|\;j =$ "+sprintf("%02d",j),'Interpreter','latex');
        else
            title("w/o detrend"+...
            "$\;|\;i =$ "+sprintf("%02d",i)+...
            "$\;|\;j =$ "+sprintf("%02d",j),'Interpreter','latex');
        end
        ylabel('$u$ (V)','Interpreter','latex');
        yyaxis right;
        ax = gca;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        if k 
            plot(t{j,1}(round(length(t{j,1})/14):2*round(length(t{j,1})/14)),...
            yf{k,1}(round(length(t{j,1})/14):2*round(length(t{j,1})/14)),...
            'Linewidth',4);
            ylabel('$y_f$ (rad)','Interpreter','latex');
        else
            plot(t{j,1}(round(length(t{j,1})/14):2*round(length(t{j,1})/14)),...
            y{i,j}(round(length(t{j,1})/14):2*round(length(t{j,1})/14)),...
            'Linewidth',4);
            ylabel('$y_f$ (rad)','Interpreter','latex');
        end
        if k == 1
            yl = ylim();
        elseif  k ~= 0
            ylim(yl);
        end
        xlabel('$t$ (s)','Interpreter','latex');
        saveas(gcf,folder+...
            sprintf("filterTuning/responseLambda_i%02d_j%02d_lambda%02d.fig",...
            i,j,k));
        saveas(gcf,folder+...
            sprintf("filterTuning/responseLambda_i%02d_j%02d_lambda%02d.fig",...
            i,j,k));
        hold off;
    end
    % Select best lambda 
    lambdaTuned = 0.8;
    % Save filter tuning data
    save(folder+sprintf("filterTuning/tuningData_i%02d_j%02d.mat",i,j),...
        'lambda','yf','lambdaTuned');
    % Save last tuned lambda
    save(folder+"filterTuning/lambdaTuned.mat",...
        'lambdaTuned');
else
   % Filtering + detrending
    lambda = 0.8;
    Afilter = [1 -lambda];
    Bfilter = (1-lambda)*[1 -1];
    tic;
    yf = cell(Nu,Nfs);
    ubi = cell(Nu,Nfs);
    for j = 1:Nfs
        for i = 1:Nu
            % do not consider first 10 seconds
            uf{i,j} = u{i,j}(round(10*fs(j)+1):end); 
            uf{i,j} = detrend(uf{i,j});
            yf{i,j} = filter(Bfilter,Afilter,y{i,j}(round(10*fs(j))+1:end));
        end
    end
    toc;
    % Save filtered responses
    save(folder+"inputSignalsResponseProcessed.mat",...
        'yf','uf');
end
%% 4.4 Identification
clear;
% ---------------------------------------------------------------------- %
% --------------------------- SET PARAMETERS --------------------------- %
RUN = 3;
% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %
% Load directory path
load(sprintf("./DATA/4_PlantModelIdentification%02d/"+"RUN.mat",RUN),...
    'folder');
% Load number of input signals
load(folder+"inputSignals.mat",...
    'Nfs','Nu');
% Load filtered responses
load(folder+"inputSignalsResponseProcessed.mat",...
        'yf','uf');
% ------------------------------------- %
% Get dimension of cell array :) porque posso e n�o me apetece fazer contas
Nidentification = 0;
for nA = 3:6
    for nB = 1:nA-1
        for nK = 1:nA-nB+1
            Nidentification = Nidentification +1;
        end
    end
end
% Find armax model for each model order, for each 
tic;
M = cell(Nu,Nfs,Nidentification);

j = 1;
% for j = 1:Nfs
    for i = 1:Nu
        count = 0;
        for nA = 3:6
            for nB = 1:nA-1
                for nK = 1:nA-nB+1
                    nC = nA;
                    count = count+1;
                    M{i,j,count} = armax([yf{i,j} uf{i,j}],[nA nB nC nK]);
                    nABCK(count,:) = [nA nB nC nK];
                end
            end
        end
    end
% end
toc;
% ------------- Save Data ------------- %
% Save filtered responses
save(folder+"identificationCellArray.mat",...
    'M','nABCK','Nidentification');
% ------------------------------------- %
%% 4.5 Model validation                
clear;
RUN = 3;
% Load directory path
load(sprintf("./DATA/4_PlantModelIdentification%02d/"+"RUN.mat",RUN),...
    'folder');
% Load number of input signals
load(folder+"inputSignals.mat",...
    'Nfs','Nu');
% Load filtered responses
load(folder+"inputSignalsResponseProcessed.mat",...
        'yf','uf');
% Load filtered responses
load(folder+"identificationCellArray.mat",...
    'M','nABCK','Nidentification');
% Get fitness of each model
tic;
fitM = zeros(Nu,Nfs,Nidentification);
j = 1;
% for j = 1:Nfs
    for i = 1:Nu
        for count = 1:Nidentification
            for k = 1:Nu
                if k == i
                    % Do not use own test as validation
                    continue;
                end
                [~,fit] = compare([yf{i,j} uf{i,j}],M{i,j,count});
                fitM(i,j,count) = fitM(i,j,count)+fit/(Nu-1);
            end  
        end
    end
    % end
toc; 

tic;
fitMcopy = fitM;
Nbest = 4;
bestFitFit = zeros(4,Nfs,Nbest);
bestFitCount = zeros(4,Nfs,Nbest);
bestFitU = zeros(4,Nfs,Nbest);
count = 0;
j = 1;
% for j = 1:Nfs
    for nA = 3:6
        for b = 1:Nbest
            idx = find(nABCK(:,1)==nA);
            bestFitFit(nA-2,j,b) = max(max(fitMcopy(:,j,idx(1):idx(end))));
            [bestFitU(nA-2,j,b),bestFitCount(nA-2,j,b)]=...
                find(fitMcopy(:,j,idx(1):idx(end))==bestFitFit(nA-2,j,b));
            bestFitCount(nA-2,j,b) = bestFitCount(nA-2,j,b)+idx(1)-1;
            fitMcopy(bestFitU(nA-2,j,b),j,bestFitCount(nA-2,j,b)) = 0;
        end
    end
% end
toc;

% ------------- Save Data ------------- %
% Save filtered responses
save(folder+"identificationFitness.mat",...
    'fitM','bestFitFit','bestFitCount','bestFitU');
% ------------------------------------- %

%% 4.5.1 Check results
clear;
RUN = 3;
% Load directory path
load(sprintf("./DATA/4_PlantModelIdentification%02d/"+"RUN.mat",RUN),...
    'folder');
% Load filtered responses
load(folder+"identificationCellArray.mat",...
    'M');
load(folder+"identificationFitness.mat",...
    'bestFitCount','bestFitU');

poleValidation = true;
nA = 5;
b = 1;
j = 1;
if poleValidation
    [den, num] = polydata(M{bestFitU(nA-2,j,b),j,bestFitCount(nA-2,j,b)});
    figure;
    zplane(num,den);  
end

% ---------------------------------------------------------------------- %
% ------------------------ CHOOSE IDENTIFICATION ----------------------- %
nA = 5; % order
b = 1; % b-th best
% ---------------------------------------------------------------------- %
% ---------------------------------------------------------------------- %

%% 4.6 Postprocessing
clear;
RUN = 3;
% Load directory path
load(sprintf("./DATA/4_PlantModelIdentification%02d/"+"RUN.mat",RUN),...
    'folder');
% Load number of input signals
load(folder+"inputSignals.mat",...
    'Nfs','Nu');
% Load filtered responses
load(folder+"inputSignalsResponseProcessed.mat",...
        'yf','uf');
% Load filtered responses
load(folder+"identificationCellArray.mat",...
    'M','nABCK','Nidentification');



%% 4.$\infty$ Batota
load("./model/barrassmodel.mat",...
    'Atrue','Btrue','Ctrue','Dtrue');
load(folder+"inputSignals.mat",...
    'T','Nfs','Ts','t','Nu','u');

j = 1;
discreteSys = ss(Atrue,Btrue,Ctrue,Dtrue,Ts(j));

