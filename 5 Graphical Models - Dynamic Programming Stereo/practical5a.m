function r=practical5a
%The goal of this practical is to investigate dynamic progamming along
%scanlines in stereo vision.  In this part of the practical, we build a
%basic dynamic programming algorithm, which will solve the same problem as
%presented in the notes.  
%In practical 5b we will apply this to the stereo problem.  


%Use the template below, filling in parts marked "TO DO".
% Note: Optionally, you may want to move the function-skeleton
%       function bestPath = dynamicProgram(unaryCosts, pairwiseCosts)...
% from below into its own new and separate dynamicProgram.m file, for 
% easier debugging and re-use.

%close any previous figures;
close all;

%define unary costs : unaryCosts(d,j) represents the cost for having
%disparity d at pixel j.
unaryCosts = [2.0 1.1 5.7 1.5 6.0 3.1 ;...
              0.8 4.8 1.0 3.0 6.9 3.3 ;...
              4.3 2.3 2.4 2.4 6.6 6.2 ;...
              6.4 0.0 6.1 0.8 7.1 2.1 ;...
              2.3 2.2 4.9 8.9 1.0 9.8 ];

%define pairwise costs:  pairwiseCosts(a,b) represents the cost for changing from 
%disparity level A at pixel j to disparity level B at pixel j+1;
pairwiseCosts = [   0   2 100 100 100;...
                2   0   2 100 100;...
              100   2   0   2 100;...
              100 100   2   0   2;...
              100 100 100   2   0];
          
%Now, use dynamic programming to find best solution.
%TO DO - fill in this routine (further down in this file).
bestPath = dynamicProgram(unaryCosts,pairwiseCosts)

%display the nodes
figure; set(gcf,'Color',[1 1 1]);
[nY nX] = size(unaryCosts);
for(cY = 1:nY)
    for (cX = 1:nX)
        plot(cX,cY,'r.','MarkerSize',20); hold on;
        text(cX+0.1,cY+0.1,sprintf('%1.1f',unaryCosts(cY,cX)));
        axis off;set(gca,'YDir','reverse');
    end;
end;

%display the best path
plot(1:nX,bestPath,'b-');

%calculate the cost of the path computed
bestPathCostUnary = sum(unaryCosts(bestPath+(0:nX-1)*nY));
bestPathCostPairwise = sum(pairwiseCosts(bestPath(1:end-1)+nY*(bestPath(2:end)-1)));
bestCost = bestPathCostUnary+bestPathCostPairwise;
fprintf('Path Cost = %3.3f\n',bestCost);


%==========================================================================
%==========================================================================

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


%Now run through each position (column)
for (cPosition = 2:nPosition)
    %run through each node (element of column)
    for (cNode = 1:nNodesPerPosition)
        %now we find the costs of all paths from the previous column to this node
        possPathCosts = zeros(nNodesPerPosition,1);
        for (cPrevNode = 1:nNodesPerPosition)
            %TO DO  - fill in elements of possPathCosts
            
        end;
        %TO DO - find the minimum of the possible paths 
        
        %TO DO - store the minimum cost in the minimumCost matrix
        
        %TO DO - store the parent index in the parents matrix
        
    end;
end;


%BACKWARD PASS

%we will now fill in the bestPath vector
bestPath = zeros(1,nPosition);

%TO DO  - find the index of the overall minimum cost from the last column and put this
%into the last entry of best path

%TO DO - find the parent of the node you just found


%run backwards through the cost matrix tracing the best patch
for (cPosition = nPosition-1:-1:1)
    %TO DO - work through matrix backwards, updating bestPath by tracing
    %parents.
    
    
end


%TO DO: REMOVE THIS WHEN YOU ARE DONE
bestPath = ceil(rand(1,nPosition)*nNodesPerPosition);






