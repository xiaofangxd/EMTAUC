function data_SBGA = SBGA(proC,disC,proM,disM,selection_process,Task, Ben, Harm,N,gen,select,maxfes,costexp)
    % We provided MATLAB implementations for
    % Evolutionary Multitasking AUC optimization (SBGA, for example).
    tic
    %% 0.Record the optimal solution matrix
    EvBestFitness = zeros(gen+1,Task.M);          % The best fitness value after each evaluation
    EvBestFitness_evn = zeros(Task.M,1);          
    Evaluations=zeros(Task.M,1);                  % Number of individual evaluations on each task
    bestSolution = zeros(Task.M,Task.D_multitask);% Best individuals
    eval_cost = zeros(floor(maxfes/(costexp*N)),1);
    RIJ = 0.5*ones(Task.M,Task.M); % transfer rates
    MIJ = ones(Task.M,Task.M); % benefit and benefit
    NIJ = ones(Task.M,Task.M); % neutral and neutral
    CIJ = ones(Task.M,Task.M); % harm and harm
    OIJ = ones(Task.M,Task.M); % neutral and benefit
    PIJ = ones(Task.M,Task.M); % benefit and harm
    AIJ = ones(Task.M,Task.M); % harm and neutral

    %% 1.Initial population
    Population = INDIVIDUAL();                    
    Population = initPOP(Population,N,Task, select);

    %% 2.Evaluate the objective function value of each individual
    for j = 1:Population.P
            [Population.factorial_costs(Population.flag == j,1),Evaluations,EvBestFitness_evn,eval_cost]=CalObj(Task,Population.rnvec(Population.flag == j,:),j,Evaluations,EvBestFitness_evn,select,maxfes,eval_cost,costexp,N);
            EvBestFitness(1,j) = EvBestFitness_evn(j);
            [~,y] = sort(Population.factorial_costs(Population.flag == j));
            [~,ranking] = sort(y);
            Population.ranking_c(Population.flag == j,:) = ranking;
            Population.ranking_o(Population.flag == j,:) = ranking;
            trnvec = Population.rnvec(Population.flag == j,:);
            bestSolution(j,:) = trnvec(y(1),:);
    end
    %% 3.Optimization process
    g = 1;
    while( sum(Evaluations'.*((Task.p/Task.p(1)).^2))<maxfes ) 
        % 3.1 Generate offspring
        Offspring = Population;
        for j = 1:Population.P
            Offspring.rnvec(Offspring.flag == j,:) = GA(Population.rnvec(Population.flag == j,:),Task,j,select,{proC,disC,proM,disM});
        end
        for j = 1:Population.P
            % 3.2 Knowledge transfer
            [~,atindex] = max(RIJ(j,[1:j-1,j+1:end]));% find transferred task
            if atindex >= j
                atindex = atindex + 1;
            end
            if rand() < RIJ(j,atindex)
                Si = floor(N*RIJ(j,atindex));% transfer quantity
                ind1 = randperm(N,Si);
                ind2 = randperm(N,Si);
                atipos = find(Offspring.flag == atindex);
                jpos = find(Offspring.flag == j);
                Offspring.rnvec(jpos(ind1),:) = Offspring.rnvec(atipos(ind2),:);
                Offspring.belongingtask(jpos(ind1),:) = atindex;
            end
            % 3.3 Evaluation of individual offspring
            [Offspring.factorial_costs(Offspring.flag == j,1),Evaluations,EvBestFitness_evn,eval_cost]=CalObj(Task,Offspring.rnvec(Offspring.flag == j,:),j,Evaluations,EvBestFitness_evn,select,maxfes,eval_cost,costexp,N);
            [~,y] = sort(Offspring.factorial_costs(Offspring.flag == j));
            [~,ranking] = sort(y);
            Offspring.ranking_c(Offspring.flag == j,:) = ranking;
            % 3.4 Merge parent,offspring and candidate transferred solutions
            Population = EnvironmentalSelection(Population,Offspring,selection_process,N,j);
            % 3.5 Record
            [~,y] = sort(Population.factorial_costs(Population.flag == j));
            [~,ranking] = sort(y);
            Population.ranking_c(Population.flag == j,:) = ranking;
            Population.ranking_o(Population.flag == j,:) = ranking;
            trnvec = Population.rnvec(Population.flag == j,:);
            bestSolution(j,:) = trnvec(y(1),:);
            EvBestFitness(g+1,j) = EvBestFitness_evn(j);
        end
        % 3.6 Update of symbiosis
        for j = 1:Population.P
            indexr = find(Offspring.flag == j & Offspring.belongingtask ~= j);
            ranking_c = Offspring.ranking_c(indexr);
            ranking_o = Offspring.ranking_o(indexr);
            for k = 1:length(indexr)
                if ranking_c(k) < N*Ben
                    if ranking_o(k) < N*Ben
                        MIJ(j,Offspring.belongingtask(indexr(k))) = MIJ(j,Offspring.belongingtask(indexr(k))) + 1;
                    elseif ranking_o(k) > N*(1-Harm)
                        PIJ(j,Offspring.belongingtask(indexr(k))) = PIJ(j,Offspring.belongingtask(indexr(k))) + 1;
                    else
                        OIJ(j,Offspring.belongingtask(indexr(k))) = OIJ(j,Offspring.belongingtask(indexr(k))) + 1;
                    end
                elseif ranking_c(k) > N*(1-Harm)
                    if ranking_o(k) > N*(1-Harm)
                        CIJ(j,Offspring.belongingtask(indexr(k))) = CIJ(j,Offspring.belongingtask(indexr(k))) + 1;
                    end
                else
                    if ranking_o(k) > N*(1-Harm)
                        AIJ(j,Offspring.belongingtask(indexr(k))) = AIJ(j,Offspring.belongingtask(indexr(k))) + 1;
                    elseif ranking_o(k) >= N*Ben && ranking_o(k) <= N*(1-Harm)
                        NIJ(j,Offspring.belongingtask(indexr(k))) = NIJ(j,Offspring.belongingtask(indexr(k))) + 1;
                    end
                end
            end
        end
        % 3.7 Update of transfer rates
        RIJ = (MIJ + OIJ + PIJ)./(MIJ + OIJ + PIJ + AIJ + CIJ + NIJ);
        disp(['SBGA Gen = ', num2str(g), ' EvBF = ', num2str(EvBestFitness(g+1,:))]);
        g = g+1;
    end

    %% Record algorithm results
    data_SBGA.wall_clock_time=toc;
    data_SBGA.EvBestFitness=EvBestFitness;
    data_SBGA.EvBestFitness_evn=EvBestFitness_evn;
    data_SBGA.Evaluations=Evaluations;
    data_SBGA.bestSolution = bestSolution;
    data_SBGA.RIJ = RIJ;
    data_SBGA.eval_cost = eval_cost;
end