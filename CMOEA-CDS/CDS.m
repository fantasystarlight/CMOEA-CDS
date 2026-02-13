function PenaltyValue = CDS(varargin)
% Constraint-diversity dominance sorting

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhe Liu
    PopObj = varargin{1};
    [N, M] = size(PopObj);
    if nargin == 2
        k  = varargin{2};
        CV = zeros(N, 1);
    else
        PopCon = varargin{2};
        k  = varargin{3};
        CV = sum(max(0,PopCon),2);
    end
    PopObj = roundn(PopObj, -10);    
    PenaltyValue = zeros(N);    
    crowdingdistance = pdist2(PopObj, PopObj);
    maxdistance = max(crowdingdistance,[],"all");
    for i = 1 : N-1
        err = PopObj(i,:) - PopObj;
        max_err = max(err, [],  2);
        min_err = min(err, [],  2);
        for j = i + 1: N
            eq = 0;
            for m = 1 : M
                if err(j, m) ~= 0
                    break;
                end
                if m == M
                    eq = 1;
                end
            end           
            if eq == 1
                if CV(i) <= CV(j)
                    PenaltyValue(j,i) = N;
                else
                    PenaltyValue(i,j) = N;
                end
            elseif min_err(j) >= 0 && CV(j) <= CV(i)
                PenaltyValue(i,j) = N;
            elseif max_err(j) <= 0 && CV(i) <= CV(j)
                PenaltyValue(j,i) = N;
            elseif CV(j) < CV(i)
                    PenaltyValue(i,j) = 1 - (crowdingdistance(i,j) / maxdistance)^(1/k);
            elseif CV(j) > CV(i)
                    PenaltyValue(j,i) = 1 - (crowdingdistance(j,i) / maxdistance)^(1/k);
            else
                    PenaltyValue(i,j) = 1 - (crowdingdistance(i,j) / maxdistance)^(1/k);
                    PenaltyValue(j,i) = 1 - (crowdingdistance(j,i) / maxdistance)^(1/k);
            end
        end
    end
end