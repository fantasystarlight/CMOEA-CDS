function Population = PopulationUpdateCDS(Population, N, k, constrained)
% Population update of CMOEA-CDS

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhe Liu
    %% Fitness calculation
    if constrained
        PenaltyValue = CDS(Population.objs, Population.cons, k);
    else
        PenaltyValue = CDS(Population.objs, k);
    end
    TotalPenaltyValue = sum(PenaltyValue, 2);    
    [~, Rank] = sort(TotalPenaltyValue);
    Next = false(1, length(Population));
    Next(Rank(1:N)) = true;    
    Population = Population(Next);    
end