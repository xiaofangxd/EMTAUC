%------------------------------------------------------------------------
% We provided MATLAB implementations for
% Evolutionary Multitasking AUC optimization (SBGA, for example).
%
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021-2022 Jing Group. You are free to use the EMTAUC for 
% research purposes. All publications which use this code should acknowledge the use of "EMTAUC" 
% and reference C. Wang, K. Wu, and J. Liu, "Evolutionary Multitasking AUC optimization,"
% in IEEE Computational Intelligence Magazine.
%
% If you have any questions, please contact this Email:xiaofengxd@126.com
%--------------------------------------------------------------------------

clc,clear
addpath(genpath(pwd)) 
% warning('off');

%% algorithm_parameter setting of EMTAUC
load ('dataset.mat');                        % Load data
for name = 1                                    % Test suites£º1-2('diabetes_scale' and 'fourclass_scale', for example) All datasets can get from https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/.
    gdata = gData(name);
    gCV = 5;                                    % Cross validation
    vIndices = crossvalind('Kfold', gdata.datNum, gCV);
    lambda = 0.001;
    times = 10;                                 % Running times
    calAUC1 = zeros(times,gCV);                 % Store the AUC value obtained by SBGA

    for g = 1:gCV
        %% Get the training samples
        Train.datDim = gdata.datDim;
        Train.datFeat = gdata.datFeat(vIndices~=g, :);
        Train.datLabel = gdata.datLabel(vIndices~=g);
        %% Get the testing samples
        datTest = gdata.datFeat(vIndices==g, :);
        labTest = gdata.datLabel(vIndices==g);

        %% Task initialization
        p = [0.1,1];                           % Sampling rate
        Task = TASK();
        Task = initTASK(Task,Train,p,lambda);

        %% Algorithm_parameter setting of MTO
        N = 10;                                   % Pop size
        gen = 100;                                % Maxgen
        maxfes = sum(gen*N*(p/p(1)).^2);          % Maximum number of function evaluations
        costexp = (p(end)/p(1)).^2;               
        proC = 1;                                 % Pc
        disC = 15;                                % the index of Pc
        proM = 1;                                 % Pm
        disM = 15;                                % the index of Pm
        selection_process = 'elitist';            % selection process: elitist¡¢roulette wheel¡¢Tournament
        select = 2;                               % 1:Unified search space£¬2:Independent search space
        Ben = 0.25;                               % Beneficial factor
        Harm = 0.5;                               % HaMTMAUCOrmful factor
        utf = 20;                                 % Interval of Dynamic adjustment strategy
        
         %% SBGA
        EvBestFitness_SBGA = zeros(gen+1,Task.M);       % Record the best fitness values for all tasks
        TotalEvaluations=zeros(Task.M,1);               
        BestFitness = zeros(times,Task.M);              % Store the optimal solution fitness values
        timesSBGA = zeros(times,1);                     % Store the running time for each run
        RIJ = zeros(Task.M,Task.M,times);           
        eval_cost_SBGA = zeros(floor(maxfes/(costexp*N)),1);
        for i = 1:times
            disp(['Times = ', num2str(i)]);
            data_SBGA(i) = SBGA(proC,disC,proM,disM,selection_process,Task,Ben,Harm,N,gen,select,maxfes,costexp);
            EvBestFitness_SBGA = EvBestFitness_SBGA + data_SBGA(i).EvBestFitness;
            BestFitness(i,:) = data_SBGA(i).EvBestFitness_evn;
            TotalEvaluations = max(TotalEvaluations,data_SBGA(i).Evaluations);
            RIJ(:,:,i) = data_SBGA(i).RIJ;
            timesSBGA(i) = data_SBGA(i).wall_clock_time;
            eval_cost_SBGA = eval_cost_SBGA + data_SBGA(i).eval_cost;
            if select == 1
                minrange = Task.Lb(end,:);
                maxrange = Task.Ub(end,:);
                x = (maxrange-minrange).*data_SBGA(i).bestSolution(Task.M,:) + minrange;
            else
                x = data_SBGA(i).bestSolution(Task.M,:);
            end
            [calAUC1(i,g), ~, ~] = fnEvaluate(datTest, labTest, x');
        end
        EvBestFitness_SBGA = EvBestFitness_SBGA./times;
        %% Record
        data.EvBestFitness_SBGA=EvBestFitness_SBGA;
        data.timesSBGA = timesSBGA;
        save(['Data/',num2str(name),'_SBGA_',num2str(g),'.mat'],'data');
    end
    dataa.AUC_SBGA = calAUC1;
    disp([ num2str(mean(mean(calAUC1,1))),'(',num2str(std(std(calAUC1,1))),')']);
    save(['Data/',num2str(name),'_SBGA_',num2str(g),'_all.mat'],'dataa');
end
