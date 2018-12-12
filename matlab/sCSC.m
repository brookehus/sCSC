%% MATLAB M-FILE TO GENERATE sCSC DENDROGRAM 
%% Brooke E. Husic, Kristy L. Schlueter-Kuck, and John O. Dabiri
%% Stanford University
%% 2017-2018

%% Instructions: Execute this m-file in a folder containing the adjacency matrix to be analyzed. To use the default code, the adjacency matrix should be in a comma-delimited ASCII file named 'adjacency.dat'



clear all
close all

%% Load adjacency matrix %%

A = dlmread('adjacency.dat',',');  



%% Compute color field %%
              
Adegree = sum(A,2);     % create column vector of row-sums of adjacency matrix

DD = diag(Adegree);     % create diagonal degree matrix with row-sums of adjacency matrix

L = DD-A;               % create graph Laplacian

[eigvec, eigval] = eig(L,DD);   % compute eigenvectors and eigenvalues of generalized eigenvalue problem. 
                                % NOTE: for large matrices, it may be
                                % more efficient to use the 'eigs' command
                                % and sort for the largest eigenvalues/vectors
              
[sortedeigs, lambdaindex] = sort(diag(eigval),'descend');  % sort eigenvalues in decending order

CSC_fullset=eigvec(:,lambdaindex);  % sort eigenvectors in order corresponding to sorted eigenvalues


%% Generate sCSC Dendrogram %%

numeigvecs = 7;         % select number of eigenvectors to include in dendrogram analysis


evec = CSC_fullset(:,1:numeigvecs); % create matrix with selected subset of eigenvectors

for col = 1:length(evec(1,:))       % loop through eigenvectors to create binary code for each state in dataset
    lnk = linkage(evec(:,col),'average');       % create agglomerative hierarchical tree using average linkage
    T = cluster(lnk,'maxclust',2);              % assign value of 1 or 2 to each state state depending on membership in bifurcated tree
    binary(:,col) = T-1;                        % shift assigment values to 0 and 1
    
end


groupdist = zeros(2^(numeigvecs-1),numeigvecs);     % initialize matrix containing length of each dendrogram branch
plotcoordsx = zeros(2^(numeigvecs-1),numeigvecs);   % initialize matrix containing x coordinates of each dendrogram node
plotcoordsy = zeros(2^(numeigvecs-1),numeigvecs);   % initialize matrix containing y coordinates of each dendrogram node
lineangle = pi/4;                                   % define angle of each dendrogram branch in plot
maxlinewidth = 15;                                  % define maximum line width of dendrogram branches

for split = 1:numeigvecs                            % loop through eigenvectors, from largest to last one included in analysis
    
    split                                           % output current dendrogram level
    row = 1;                                        % initialize level of dendrogram
    
    for combos = 0:2:2^(split)-2                    % loop through binary codes
        group1 = fliplr(de2bi(combos,split));       % create binary code for group 1 of a pair of branches
        group2 = fliplr(de2bi(combos+1,split));     % create binary code for group 2 of a pair of branches
        
        g1label = dec2bin(combos,split);            % create string-compatible version of binary code for group 1 for plotting
        g2label = dec2bin(combos+1,split);          % create strong-compatible version of binary code for group 2 for plotting
        
        
        Xgroup = evec(:,split);                     % initialize eigenvector for current level of dendrogram
        Xgroup(find(~ismember(binary(:,1:split),group1,'rows') & ~ismember(binary(:,1:split),group2,'rows'))) = [];     % remove eigenvector elements not in in group 1 or group 2
             
        Lgroup = L;                                 % initialize graph Laplacian    
        Lgroup(:,find(~ismember(binary(:,1:split),group1,'rows') & ~ismember(binary(:,1:split),group2,'rows'))) = [];   % remove columns not in group 1 or group 2
        Lgroup(find(~ismember(binary(:,1:split),group1,'rows') & ~ismember(binary(:,1:split),group2,'rows')),:) = [];   % remove rows not in group 1 or group 2
        
        if ~isempty(find(ismember(binary(:,1:split),group1,'rows')))  % if group 1 is occupied by any states
            Z1 = Xgroup'*Lgroup*Xgroup;                                             % calculate dissimilarity metric for states in groups 1 and 2
            Z1count = numel(find(ismember(binary(:,1:split),group1,'rows')));       % calculate number of states in group 1
        else
            Z1 = NaN;                                                               % dissimilarity metric is undefined if group 1 is unoccupied
        end
        
        if ~isempty(find(ismember(binary(:,1:split),group2,'rows')))  % if group 2 is occupied by any states
            Z2 = Xgroup'*Lgroup*Xgroup;                                             % calculate dissimilarity metric for states in groups 1 and 2
            Z2count = numel(find(ismember(binary(:,1:split),group2,'rows')));       % calculate number of states in group 2
        else
            Z2 = NaN;                                                               % dissimilarity metric is undefined if group 2 is unoccupied
        end
        
                
        groupdist(row,split) = Z1;                  % record dissimilarity metric for group 1
        groupdist(row+1,split) = Z2;                % record dissimilarity metric for group 2
        
        
        
        
        figure(1)                                   % open figure
        hold on                                     % retain current plot when adding new features
        
        xrel1 = -Z1*sin(lineangle);                 % relative x-coordinate for endpoint of dendrogram branch of group 1
        xrel2 = Z2*sin(lineangle);                  % relative x-coordinate for endpoint of dendrogram branch of group 2
        
        yrel1 = -Z1*cos(lineangle);                 % relative y-coordinate for endpoint of dendrogram branch of group 1
        yrel2 = -Z2*cos(lineangle);                 % relative y-coordinate for endpoint of dendrogram branch of group 2
        
        if split == 1                               % for first dendrogram level
            plotcoordsx(row,split) = xrel1;         % absolute x-coordinate for endpoint of dendrogram branch of group 1
            plotcoordsx(row+1,split) = xrel2;       % absolute x-coordinate for endpoint of dendrogram branch of group 2
            
            plotcoordsy(row,split) = yrel1;         % absolute y-coordinate for endpoint of dendrogram branch of group 1
            plotcoordsy(row+1,split) = yrel2;       % absolute y-coordinate for endpoint of dendrogram branch of group 2
            
            
            plot([0 plotcoordsx(row,split)],[0 plotcoordsy(row,split)],'-o','Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',maxlinewidth*(Z1count/length(A(:,1))))          % plot dendrogram branch of group 1
            plot([0 plotcoordsx(row+1,split)],[0 plotcoordsy(row+1,split)],'-o','Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',maxlinewidth*(Z2count/length(A(:,1))))      % plot dendrogram branch of group 2
            
            text(plotcoordsx(row,split),mean([0 plotcoordsy(row,split)]),g1label,'FontSize',12)                                 % label branch of group 1 with corresponding binary code
            text(plotcoordsx(row+1,split),mean([0 plotcoordsy(row+1,split)]),g2label,'FontSize',12)                             % label branch of group 2 with corresponding binary code
            
            text(plotcoordsx(row,split),plotcoordsy(row,split),['     ' num2str(Z1count)],'Color','r','FontSize',12)            % label end node of group 1 with corresponding number of states
            text(plotcoordsx(row+1,split),plotcoordsy(row+1,split),['     ' num2str(Z2count)],'Color','r','FontSize',12)        % label enb node of group 2 with corresponding number of states
            
        else
            plotcoordsx(row,split) = xrel1+plotcoordsx((row+1)/2,split-1);      % absolute x-coordinate for endpoint of dendrogram branch of group 1
            plotcoordsx(row+1,split) = xrel2+plotcoordsx((row+1)/2,split-1);    % absolute x-coordinate for endpoint of dendrogram branch of group 2
            
            plotcoordsy(row,split) = yrel1+plotcoordsy((row+1)/2,split-1);      % absolute y-coordinate for endpoint of dendrogram branch of group 1
            plotcoordsy(row+1,split) = yrel2+plotcoordsy((row+1)/2,split-1);    % absolute y-coordinate for endpoint of dendrogram branch of group 2
                
            plot([plotcoordsx((row+1)/2,split-1) plotcoordsx(row,split)],[plotcoordsy((row+1)/2,split-1) plotcoordsy(row,split)],'-o','Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',maxlinewidth*(Z1count/length(A(:,1))))        % plot dendrogram branch of group 1
            plot([plotcoordsx((row+1)/2,split-1) plotcoordsx(row+1,split)],[plotcoordsy((row+1)/2,split-1) plotcoordsy(row+1,split)],'-o','Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',maxlinewidth*(Z2count/length(A(:,1))))    % plot dendrogram branch of group 2
            
            text(plotcoordsx(row,split),mean([plotcoordsy((row+1)/2,split-1) plotcoordsy(row,split)]),g1label,'FontSize',12)        % label branch of group 1 with corresponding binary code
            text(plotcoordsx(row+1,split),mean([plotcoordsy((row+1)/2,split-1) plotcoordsy(row+1,split)]),g2label,'FontSize',12)    % label branch of group 2 with corresponding binary code
            
            text(plotcoordsx(row,split),plotcoordsy(row,split),['     ' num2str(Z1count)],'Color','r','FontSize',12)                % label end node of group 1 with corresponding number of states
            text(plotcoordsx(row+1,split),plotcoordsy(row+1,split),['     ' num2str(Z2count)],'Color','r','FontSize',12)            % label enb node of group 2 with corresponding number of states
              
        end
        
        
        
        row = row+2;        % advance dendrogram level counter
        
    end
end

axis equal                  % set plot axes equal


%% Export sCSC Dendrogram distance matrix %%

dlmwrite('sCSC_Dendrogram.dat',groupdist,',')


%% Export indices of states in a selected group in sCSC Dendrogram %%

group = [1 1 1 1]               % select binary code of group of interest. Note: include spaces between digits of binary code

split = numel(group);           % determine level of binary code in dendrogram

groupidx = find(ismember(binary(:,1:split),group,'rows'));      % find indices corresponding to selected binary code

dlmwrite(['sCSC_GroupIndices_' fliplr(dec2bin(bi2de(group),split)) '.dat'],groupidx,',')   % save indices to file


