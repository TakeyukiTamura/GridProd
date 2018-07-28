function [M,Brange,Trange,B,TMY,MB]=...
    BTconstraintSearch(model,targetMet,GUR,OUR,minGrowth,glucoseRxnID,oxygenRxnID,biomassRxnID,P)
%BTconstraintSearch is a function of GridProd
%BTconstraintSearch obtains the constraints for the target production and
%biomass objective function, and returns the value of achieved target
%production.
%
%INPUTS
% model    COBRA model structure
% targetMet   target metabolite 
% GUR     glucose uptake ratio
% OUR     oxygen uptake ratio
% minGrowth    The minimum value of biomass objective fucnction that
%              GridProd must satisfy.
% glucoseRxnID     ID of glucose reaction in the COBRA model
% oxygenRxnID      ID of oxygen reaction in the COBRA model
% biomassRxnID     ID of biomass reaction in the COBRA model
% P           paramters for the grid size (described as $P^{-1}$ in the manuscript)
% 
%OUTPUTS
% M     achieved target metabolite production
% Brange     ID for the vertical axis of the grid that achieves the highest
%            target metabolite production
% Trange     ID for the horizontal axis of the grid that achieves the highest
%            target metabolite production
% B          The value of the biomass objective function when the highest
%            target metabolite production is achieved
% TMY         Theoretical Maximum yield (TMPR: Theoretical Maximum Production rate)
% MB          Maximum Biomass (TMGR: Theoretical Maximum Growth Rate)
%
%  Jul. 20, 2017   Takeyuki TAMURA

model.lb(glucoseRxnID)=-GUR;
model.lb(oxygenRxnID)=-OUR;
target=findMetIDs(model,targetMet);
[model2,rxnIDexists]=addReaction(model,'Transport',{targetMet},[-1])
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
opt2=optimizeCbModel(model2);
MB=opt2.f;
model3=model2;
model3.c(biomassRxnID)=0;
model3.c(targetRID)=1;
model3.lb(biomassRxnID)=minGrowth;
opt2=optimizeCbModel(model3);
TMY=opt2.f;

model3=model2;
model3.lb(biomassRxnID)=minGrowth;
lp.f=[zeros(n,1); ones(n,1)];
lp.A=[-eye(n) -eye(n);eye(n) -eye(n)];
lp.b=zeros(2*n,1);
lp.Aeq=[model3.S zeros(m,n)];
lp.beq=zeros(m,1);
lp.lb=[model3.lb; zeros(n,1)];
lp.ub=[model3.ub; 999999*ones(n,1)];

for i=1:P
    i
    biomassLB=(MB/P)*(i-1);
    biomassUB=(MB/P)*i;
    for j=1:P
        targetLB=(TMY/P)*(j-1);
        targetUB=(TMY/P)*j;
        [table(i,j),opt4biomass(i,j),opt4target(i,j),opt4f(i,j),opt5biomass(i,j)]=...
            integrate(model2,model3,targetRID,biomassLB,biomassUB,targetLB,targetUB,n,biomassRxnID,lp);
    end
end
table2=table;
table2(opt5biomass<minGrowth)=0;
[M,I] = max(table2(:));
[Brange,Trange] = ind2sub(size(table2),I);
B=opt5biomass(Brange,Trange);


filename=sprintf('results/BTconditionSearch_%d.mat',target);
save(filename);
end

