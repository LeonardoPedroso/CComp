function y = simOpenLoop(t,u,Ts)
    barrassmodel = load("./model/barrassmodel.mat");
    sim('./model/barra1.slx',t(end),...
        simset('DstWorkspace','current','SrcWorkspace','current'));
end

