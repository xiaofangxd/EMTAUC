function [Population,num,numm] = tfES(Population,Offspring,tfsol,tfflag,factorial_costs,selection_process,N,M,NUM)
% EnvironmentalSelection
% Input:population chromes and objective values, offspring chromes and
% objective values, candidate transferred solutions' chromes and objective
% values, selection_process (option: elitist,roulette wheel), task number
% M.
% Output: population chromes and objective values, the number of effect
% transfer solutions.
%--------------------------------------------------------------------------
    nvec = [Population.rnvec(Population.flag == M,:);Offspring.rnvec(Offspring.flag == M,:);tfsol];
    fitness = [Population.factorial_costs(Population.flag == M,1);Offspring.factorial_costs(Offspring.flag == M,1);factorial_costs];
    if strcmp(selection_process,'elitist')
        [~,index]=sort(fitness);
    elseif strcmp(selection_process,'roulette wheel')
        index = RouletteWheelSelection(N,fitness);
    elseif strcmp(selection_process,'Tournament')
        index = TournamentSelection(2,N,fitness);
    end
    Population.rnvec(Population.flag == M,:) = nvec(index(1:N),:);
    Population.factorial_costs(Population.flag == M,1) = fitness(index(1:N),1);
    
    num = sum(index(1:N) > 2*N);
    ttfflag = tfflag(index(index(1:N) > 2*N)-2*N);
    numm = zeros(1,NUM);
    for i = 1:NUM
        numm(i) = numel(find(ttfflag==i));
    end
%     disp(['the number of effect transfer solutions: ' num2str(num)]);
end