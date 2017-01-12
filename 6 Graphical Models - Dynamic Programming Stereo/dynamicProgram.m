%the goal of this routine is to return the minimum cost dynamic programming
%solution given a set of unary and pairwise costs

function bestPath = dynamicProgram(unaryCosts, pairwiseCosts)


%count number of positions (i.e. pixels in the scanline), and nodes at each
%position (i.e. the number of distinct possible disparities at each position)
[nNodesPerPosition nPosition] = size(unaryCosts);

%define minimum cost matrix - each element will eventually contain
%the minimum cost to reach this node from the left hand side.
%We will update it as we move from left to right
minimumCost = zeros(nNodesPerPosition, nPosition);

%define parent matrix - each element will contain the (vertical) index of
%the node that preceded it on the path.  Since the first column has no
%parents, we will leave it set to zeros.
parents = zeros(nNodesPerPosition, nPosition);

%FORWARD PASS

%TO DO:  fill in first column of minimum cost matrix
minimumCost(:,1) = unaryCosts(:,1);

%Now run through each position (column)
for (cPosition = 2:nPosition)
    %run through each node (element of column)
    for (cNode = 1:nNodesPerPosition)
        %now we find the costs of all paths from the previous column to this node
        possPathCosts = zeros(nNodesPerPosition,1);
        
        for (cPrevNode = 1:nNodesPerPosition)
            %TO DO  - fill in elements of possPathCosts
            possPathCosts(cPrevNode) = unaryCosts(cNode,cPosition) + pairwiseCosts(cPrevNode,cNode) + minimumCost(cPrevNode,cPosition-1);
        end;
        
        %TO DO - find the minimum of the possible paths 
        [minPathCost, indexMin] = min(possPathCosts);
        %TO DO - store the minimum cost in the minimumCost matrix
        minimumCost(cNode,cPosition) = minPathCost;
        %TO DO - store the parent index in the parents matrix
        parents(cNode,cPosition) = indexMin;
    end;
end;


%BACKWARD PASS

%we will now fill in the bestPath vector
bestPath = zeros(1,nPosition);

%TO DO  - find the index of the overall minimum cost from the last column and put this
%into the last entry of best path
[temp,indexMinLast] = min(minimumCost(:,nPosition));
bestPath(nPosition) = indexMinLast;

%TO DO - find the parent of the node you just found
% parent = parents(bestPath(nPosition),nPosition); %% not used

%run backwards through the cost matrix tracing the best patch
for (cPosition = nPosition-1:-1:1)
    %TO DO - work through matrix backwards, updating bestPath by tracing
    %parents.
    bestPath(cPosition) = parents(bestPath(cPosition+1),cPosition+1);    
end