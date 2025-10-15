function Population = ArchiveUpdate(Population, N)
% The Constrained based search of CCMO-CDS

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
    PopObj = Population.objs;
    [P, ~] = size(PopObj);
    PopCon = Population.cons;
    CV = sum(max(0,PopCon),2);  
    CV(CV < 0) = 0;
    Dominated = zeros(P);
    for i = 1: P-1
        for j = i+1: P
            Domination = any(PopObj(i,:)>PopObj(j,:))-any(PopObj(i,:)<PopObj(j,:));
            if Domination == 1 && CV(i) >= CV(j)
                Dominated(i,j) = P;
            elseif Domination == -1 && CV(i) <= CV(j)
                Dominated(j,i) = P;
            elseif CV(i) > CV(j)
                Dominated(i,j) = 1;
            elseif CV(i) < CV(j)
                Dominated(j,i) = 1;            
            end
        end
    end
    PenaltyValue = sum(Dominated, 2);
    Next = PenaltyValue == 0;
    %% Ranking and environmental selection
    if sum(Next) <= N
        [~,Rank] = sort(CV);
        Next(Rank(1:N)) = true;
    elseif sum(Next) > N
        Del  = Truncation(Population(Next).objs,sum(Next)-N);
        Temp = find(Next);
        Next(Temp(Del)) = false;
    end
    Population = Population(Next);
end

function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

    %% Truncation
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end