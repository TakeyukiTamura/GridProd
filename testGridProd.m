function [] = testGridProd()
%'testGridProd' tests whether GridProd appropriately works in user's
% environment.
%
% Jul. 20, 2017   Takeyuki TAMURA
%
load('ecoli_core_model.mat');
GridProd(model,{'succ[c]'},'EX_glc(e)','EX_o2(e)','Biomass_Ecoli_core_w_GAM','P',5);
load('analyzeResult.mat');
model6=changeRxnBounds(model2,rxnList,0,'b');
opt6=optimizeCbModel(model6);
if opt6.x(targetRID)==M
    display('testGridProd worked successfully.')
else
    error('testGridProd did not work successfully') 
end
save('testGridProd.mat');
end

