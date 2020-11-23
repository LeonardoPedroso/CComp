%% Design and validate the observer
clear

parent = "observer_results";
modelName = "nA3";
RUN = 1;

load("identifications/identification_"+modelName,'AId','BId','CId')

pulseStart = 5;
pulseAmplitude = 1;
pulseDuration = 5;
simT = 20;

out = sim("design_regulator",simT);

plot(out.tout,out.y)