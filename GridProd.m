function [targetProduction, minFlux, maxFlux, blockedRxns, usedRxns, biomass]...
    = GridProd(varargin)
%GridProd is a function of GridProd that identifies a set of unused
%reactions for production of target metabolites.
%
%function [targetProduction, minFlux, maxFlux, blockedRxns, usedRxns, biomass]...
%    = GridProd(model, targetMet,  glucoseRxn, oxygenRxn, biomassRxn, options)
%
%INPUTS
% model     COBRA model structure containing the following required fields to perform GridProd.
%   rxns                    Rxns in the model
%   mets                    Metabolites in the model
%   S                       Stoichiometric matrix (sparse)
%   b                       RHS of Sv = b (usually zeros)
%   c                       Objective coefficients
%   lb                      Lower bounds for fluxes
%   ub                      Upper bounds for fluxes
%   rev                     Reversibility of fluxes
%
% targetMet   target metabolites
%             (e.g., {'succ_c','glu__D_c'} )
% glucoseRxn  Reaction representing glucose uptake (e.g., EX_glc__D_e)
% oxygenRxn   Reaction representing oxygen uptake (e.g., EX_o2_e)
% biomassRxn  Reaction representing biomass objective function
%                       (e.g., BIOMASS_Ec_iAF1260_core_59p81M)
%
%OPTIONAL INPUTS
% GUR    Glucose uptake ratio (Default: 10)
% OUR    Oxygen uptake ratio (Default: 5)
% minGrowth   The minimum value of biomass objective function that the
%             designed strain must achieve. (Default: 0.05)
% P           paramters for the grid size (described as $P^{-1}$ in the manuscript, Default:25)
%
%OUTPUTS
% targetProduction   The production rate of target metabolites achieved by the
%                    designed metabolic network.
% minFlux   The minimum values of the target metabolite production
%           obtained by FVA.
% maxFlux   The maximum values of the target metabolite production
%           obtained by FVA.
% blockedRxns   A set of unused reaction that achieves the value 
%               of targetProduction 
% usedRxns     A set of reactions that is not included in blockedRxns
% biomass      The value of biomass objective function when blockedRxns
%              is not used.
%
%
% Jul. 20, 2017   Takeyuki TAMURA
%

s=size(varargin,2);
if size(varargin,2)<5
    error('''model'',''target'',''glucoseRxn'',''oxygenRxn''.''biomassRxn'' must be specified.')
end
model=varargin{1};
targetMet=varargin{2};
glucoseRxnID=findRxnIDs(model,varargin{3});
if glucoseRxnID==0
    error('invalid glucoseRxn name')
end
oxygenRxnID=findRxnIDs(model,varargin{4});
if oxygenRxnID==0
    error('invalid oxygenRxn name')
end
biomassRxnID=findRxnIDs(model,varargin{5});
if biomassRxnID==0
    error('invalid biomassRxn name')
end
GUR=10;
OUR=5;
minGrowth=0.05;
P=25;
for i=3:floor(s/2)
    if strcmp(varargin{2*i},'GUR')==1
        GUR=varargin{2*i+1};
    elseif strcmp(varargin{2*i},'OUR')==1
        OUR=varargin{2*i+1};
    elseif strcmp(varargin{2*i},'minGrowth')==1
        minGrowth=varargin{2*i+1};
    elseif strcmp(varargin{2*i},'P')==1
        P=varargin{2*i+1};
    else
        error('Options must be a subset of {type, GUR, OUR, minGrowth}')
    end
end


for i=1:size(targetMet,2)
       [T(i),Brange(i),Trange(i),B(i),TMY,MB]=...
           BTconstraintSearch(model,targetMet{i},...
           GUR,OUR,minGrowth,glucoseRxnID,oxygenRxnID,biomassRxnID,P); 
    [M(i),I(i)]=max(T(i,:));
    [targetProduction(i),biomass(i),minFlux(i),maxFlux(i),blockedRxns{i},usedRxns{i}]=...
        analyzeResult(model,targetMet{i},Brange(i,I(i)),Trange(i,I(i)),...
        GUR,OUR,minGrowth,glucoseRxnID,oxygenRxnID,biomassRxnID,TMY,MB,P);
    targetProduction(i)
    blockedRxns{i}
end

save('GridProd.mat');
end

