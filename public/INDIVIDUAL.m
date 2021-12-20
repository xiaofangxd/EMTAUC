classdef INDIVIDUAL
    % This class contains P populations.
    properties
        P; % the number of pop
        flag; % the index of pop
        rnvec; % (genotype)--> decode to find design variables --> (phenotype) 
        factorial_costs;
        belongingtask;
        ranking_o; 
        ranking_c;
    end    
    methods        
        function object = initPOP(object,N,Task,select)
            object.P = Task.M;
            object.rnvec = rand(object.P*N,Task.D_multitask);
            object.factorial_costs = inf*ones(object.P*N,1);
            object.flag = zeros(object.P*N,1);
            for i=1:object.P
                object.flag((i-1)*N+1:i*N,1) = i;
                if select == 2
                    minrange = Task.Lb(i,1:Task.Tdims(i));
                    maxrange = Task.Ub(i,1:Task.Tdims(i));
                    object.rnvec(object.flag == i,1:Task.Tdims(i)) = repmat(maxrange-minrange,[N,1]).*object.rnvec(object.flag == i,1:Task.Tdims(i))+repmat(minrange,[N,1]);
                end
            end
            object.belongingtask = object.flag;
            object.ranking_o = repmat((1:N)',object.P,1);
            object.ranking_c = repmat((1:N)',object.P,1);
        end
    end
end