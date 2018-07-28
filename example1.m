function [] = example1()
%  load the genome-scale COBRA model.
load('iAF1260.mat')
%
% GridProd is employed.
% The target metabolite is "acser_e", which corresponds to "EX_acser(e)".
% The glucose reaction is 'EX_glc__D_e'.
% The oxygen reaction is 'EX_o2_e'.
% The biomass objective fucntion reaction is
% 'BIOMASS_Ec_iAF1260_core_59p81M'.
%
% Since the options are not specified, 
% the glucose uptake ratio is 10,
% the oxygen uptake ratio is 5,
% the minimum growth ratio is 0.05 as defaults,
% the size of each grid is (TMGR/25)(TMPR/25). 
%
[targetProduction, minFlux, maxFlux, blockedRxns, usedRxns, biomass]=...
GridProd(iAF1260,{'acser_e'},'EX_glc__D_e','EX_o2_e','BIOMASS_Ec_iAF1260_core_59p81M');
minFlux
maxFlux
%
% 7.1031 is obtained as the minimum and maximum production rates for O-Acetyl-L-serine  production.
% The corresponding growth rate is 0.1496. 
%
paintGrid('results/analyzeResult_339.mat')
%
% 'paintGrid' shows the heatmap that shows the levels of O-Acetyl-L-serine
%  production rate for each grid.
% Note that 339 is ID of O-Acetyl-L-serine in iAF1260.
save('example1.mat');
%
% Jul. 20, 2017, Takeyuki TAMURA
%
end

