function [y,opt5biomass,minFlux5,maxFlux5,rxnList,usedList] = ...
    analyzeResult(model,targetMet,Brange,Trange,...
    GUR,OUR,minGrowth,glucoseRxnID,oxygenRxnID,biomassRxnID,TMY,MB,P)
% 'analyzeResult' is a function of GridProd.
% 'analyzeResult' outputs various information corresponding to specified grids. 
%
%INPUTS
% model       COBRA model structure
% targetMet   target metabolite
% Brange      ID for the vertical axis of the grid that achieves the highest
%             target metabolite production
% Trange      ID for the horizontal axis of the grid that achieves the highest
%             target metabolite production
% GUR         glucose uptake ratio
% OUR         oxygen uptake ratio
% minGrowth   The minimum value of biomass objective function that the
%             designed strain must achieve. 
% glucoseRxn  ID of the reaction representing glucose uptake in the model
% oxygenRxn   ID of the reaction representing oxygen uptake in the model
% biomassRxn  ID of the reaction representing biomass objective function in
%             the model
% TMY         Theoretical Maximum yield (TMPR: Theoretical Maximum Production rate)
% MB          Maximum Biomass (TMGR: Theoretical Maximum Growth Rate)
% P           paramters for the grid size (described as $P^{-1}$ in the manuscript, Default:25)
%                       
%OUTPUTS
% y            the value of target metabolite production in the second LP. If
%              the LP is not feasible, -99999 is assigned.
% opt5biomass  the value of biomass objective function in the second LP.
% minFlux5     the minimum values of the target metabolite production
%              obtained by FVA.
% maxFlux5     the maximum values of the target metabolite production
%              obtained by FVA.
% rxnList      a set of reactions that are not used in the second LP 
% usedList     a set of reactions that are used in the second LP
%
% Jul. 20, 2017   Takeyuki TAMURA
%
target=findMetIDs(model,targetMet);
filename=sprintf('results/BTconditionSearch_%d.mat',target);
load(filename);
model.lb(glucoseRxnID)=-GUR;
model.lb(oxygenRxnID)=-OUR;
[model2,rxnIDexists]=addReaction(model,'Transport',{targetMet},[-1]);
m=size(model2.mets,1);
n=size(model2.rxns,1);
model2.S(target,n)=-1;
model2.ub(n)=999999;
model2.lb(n)=0;
model2.rev(n)=0;
targetRID=n;
if isempty(rxnIDexists)==0
    model2=model;
    m=size(model2.mets,1);
    n=size(model2.rxns,1);
    targetRID=rxnIDexists;
end

model3=model2;
model3.c(biomassRxnID)=0;

for i=1:n
    if model3.rev(i)==0
        model3.c(i)=-1;
    elseif opt2.x(i)>0
        model3.c(i)=-1;
        model3.lb(i)=0;
    elseif opt2.x(i)<0
        model3.c(i)=1;
        model3.ub(i)=0;
    end
end



%%%%%%%%%%%%%%%%
biomassLB=(MB/P)*(Brange-1);
biomassUB=(MB/P)*Brange;
targetLB=(TMY/P)*(Trange-1);
targetUB=(TMY/P)*Trange;
%%%%%%%%%%%%
model4=model3;
model4.lb(biomassRxnID)=biomassLB;
model4.ub(biomassRxnID)=biomassUB;
model4.lb(targetRID)=targetLB;
model4.ub(targetRID)=targetUB;
opt4=optimizeCbModel(model4);
if opt4.stat~=1
    y=-99999;
    opt4biomass=-99999;
    opt4target=-99999;
    opt4f=-99999;
    opt5biomass=-99999;
    minFlux5=0;
    maxFlux5=0;
    usedRxns=[1:n];
    blockedRxns=[];
    rxnList=model4.rxns(blockedRxns);
    usedList=model4.rxns(usedRxns);
    return
end
opt4biomass(i,1)=opt4.x(biomassRxnID);
opt4target(i,1)=opt4.x(targetRID);
opt4f=opt4.f;
usedRxns=find(abs(opt4.x)>=0.0000001);
blockedRxns=setdiff(([1:n])',usedRxns);
rxnList=model4.rxns(blockedRxns);
usedList=model4.rxns(usedRxns);
model5=changeRxnBounds(model2,rxnList,0,'b');
opt5=optimizeCbModel(model5);
[minFlux5,maxFlux5]=fluxVariability(model5,100,'max',model5.rxns(targetRID));
switch opt5.stat
    case 1
        y=opt5.x(targetRID);
        opt5biomass=opt5.x(biomassRxnID);
    otherwise
        y=-999999;
        opt5biomass=-999999;
end
filename=sprintf('results/analyzeResult_%d.mat',target);
save(filename);
save('analyzeResult.mat');
end

