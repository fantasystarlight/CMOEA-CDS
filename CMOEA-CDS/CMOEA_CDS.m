classdef CMOEA_CDS < ALGORITHM
% <2026> <multi/many> <real/binary/permutation><constrained/none>
% CMOEA with constraint-diversity dominance sorting
% t ---  2 --- Mutation distribution control parameter
% k --- 10 --- Diversity control parameter

%------------------------------- Reference --------------------------------

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhe Liu

    methods
        function main(Algorithm,Problem)
        %% Generate random population
            [t,k] = Algorithm.ParameterSet(2,10);
            Archive          = Problem.Initialization();
            CF = Archive(1:floor(end/2));
            CR = Archive(floor(end/2)+1:2*floor(end/2)); 
        %% External archive initialization
            while Algorithm.NotTerminated(Archive)
        %% Optimization
                MatingPool_UCP = [randperm(length(CF)),randperm(length(CF))];
                MatingPool_CP = [randperm(length(CR)),randperm(length(CR))];
                Offspring1 = OperatorDE(Problem, CF, CF(MatingPool_UCP(1:end/2)), CF(MatingPool_UCP(end/2+1:end)),{1,0.5,1,20*rand()^t});
                Offspring2 = OperatorDE(Problem, CR, CR(MatingPool_CP(1:end/2)), CR(MatingPool_CP(end/2+1:end)),{1,0.5,1,20});
                CF = PopulationUpdateCDS([CF, Offspring1, Offspring2], floor(Problem.N/2), k, false);  
                CR = PopulationUpdateCDS([CR, Offspring1, Offspring2], floor(Problem.N/2), k, true);
                Archive = ArchiveUpdate([Archive, CR], Problem.N);  
            end
        end
    end
end

