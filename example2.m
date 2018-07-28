function [] = example2()
%  load the ecoli core model.
load('ecoli_core_model.mat')
%
% GridProd is employed.
% The target metabolite is "succ[c]".
% The glucose reaction is 'EX_glc(e)'.
% The oxygen reaction is 'EX_o2(e)'.
% The biomass objective fucntion reaction is
% 'Biomass_Ecoli_core_w_GAM'.
%
% The options are specified as follows.
% the glucose uptake ratio is 8,
% the oxygen uptake ratio is 18.5,
% the minimum growth ratio is 0.01.
% P is 10.
%
[targetProduction, minFlux, maxFlux, blockedRxns, usedRxns, biomass]=...
GridProd(model,{'succ[c]'},'EX_glc(e)','EX_o2(e)','Biomass_Ecoli_core_w_GAM',...
    'GUR',8,'OUR',18.5,'minGrowth',0.01,'P',10);
minFlux
maxFlux
%
%  9.9904 is obtained as the minimum and maximum production rates for succinate.
% 
save('example2.mat');
paintGrid('results/analyzeResult_69.mat')
%
% 'paintGrid' shows the heatmap that shows the levels of succinate
%  production for each grid.
% Note that 69 is ID of succinate in the ecoli core model.
save('example2.mat');
%
% Jul. 20, 2017, Takeyuki TAMURA
%
end

