%% evaluate the results
function [nAuc, vTpr, vTnr] = fnEvaluate(tDat, tLabel, optW)

vPy = tDat * optW;
[vTpr, vTnr, sInfo] = vl_roc(tLabel', vPy');
%vl_roc(tel', py');

nAuc = sInfo.auc;
