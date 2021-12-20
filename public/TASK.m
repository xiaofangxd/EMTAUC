classdef TASK    
    %This class contains all task information and needs to be initialized
    %by initTASK.
    properties
        M;% Number of tasks
        Tdims;% Dimension of tasks
        D_multitask;% Unified search space
        Lb;% Task lower bounds
        Ub;% Task upper bounds
        fun;% TASK functions
        p;
        datFeat;
        lambda;
    end    
    methods        
        function object = initTASK(object,gData,p,lambda)
                object.p = p;
                object.M = length(p);
                object.Tdims = zeros(object.M,1);
                datFeat1 = gData.datFeat(gData.datLabel==1,:);% Positive samples
                datFeat2 = gData.datFeat(gData.datLabel==-1,:);% Negative samples
                object.lambda = lambda;
                for i = 1:object.M
                    datNum1 = round(size(datFeat1,1)*p(i));
                    datNum2 = round(size(datFeat2,1)*p(i));
                    object.datFeat(i).datFeat11 = datFeat1(randperm(size(datFeat1,1),datNum1),:);
                    object.datFeat(i).datFeat22 = datFeat2(randperm(size(datFeat2,1),datNum2),:);
                    object.Tdims(i) = gData.datDim;
                    object.fun(i).fnc=@(x)calAUC(x,object.datFeat(i).datFeat11,object.datFeat(i).datFeat22,object.lambda);
                end
                object.D_multitask = max(object.Tdims);
                object.Lb = -ones(object.M,object.D_multitask);
                object.Ub = ones(object.M,object.D_multitask);
        end  
    end
end