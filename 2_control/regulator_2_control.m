%% Design and validate the regulator
clear

parent = "regulator_results";
modelName = "nA3";
RUN = 1;

load("identifications/identification_"+modelName,'AId','BId','CId',...
    'DId','TsId')

CPrime = [CId; eye(size(AId,1))];
DPrime = [DId; zeros(size(AId,1),1)];
ySelector = [1 zeros(1,size(AId,1))];
xSelector = eye(size(AId,1)+1);
xSelector = xSelector(2:end,:);

pulseStart = 5;
pulseAmplitude = 1;
pulseDuration = 5;
simT = 20;
K = ones(1,size(AId,1));

out = sim("design_regulator",simT);

figure()
plot(out.tout,out.y)