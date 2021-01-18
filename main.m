
clearvars -except

%% Load Required Data Files
load('Network_AMS_MTB_Revised1.mat'); %loads as amsterdam
load('NetworkDemand_AMS_MTB.mat'); %loads as amsterdamdemand

%% Nodes

%Gets a list of nodes in cell array "nodes"
%Input used is cell array of lines and their respective stops in sequence

numnetworks=2;

nodes={};
nodes_network={};
nodes_network=cell(numnetworks,1);

transferpoints=[];

numlines=zeros(1,numnetworks);
maxnodes=zeros(1,numnetworks);
%Store lines and find maximum nodes
for i=1:numnetworks
    numlines(i)=size(amsterdam{i}{1},1);
    maxnodes(i)=size(amsterdam{i}{1},2)-1;
    nodes_network{i}={};
end

%Create list of nodes and transfer points
for n=1:numnetworks
    for i=1:numlines(n)
        for j=1:maxnodes(n)
            % checks if node being inserted into "nodes" 
            % is not already present
            if all(strcmp(amsterdam{n}{1}{i,j+1},nodes)==0)
                nodes=[nodes amsterdam{n}{1}{i,j+1}];
                nodes_network{n,1}=[nodes_network{n,1} ...
                    amsterdam{n}{1}{i,j+1}];
            else
                % checks if node being inserted into "nodes" 
                % is not already present
                if mod(i,2)~=0 && ...
                        any(transferpoints==...
                        find(strcmp(nodes,amsterdam{n}{1}{i,j+1})))==0 
                %if a node is present in two lines it is a transfer point
                    transferpoints=[transferpoints ...
                        find(strcmp(nodes,amsterdam{n}{1}{i,j+1}))];
                end
            end
        end
    end
end

numnodes=size(nodes,2); %finding the number of nodes
odmat=zeros(numnodes);%arranging imported odmatrix as per nodes order

for i=1:size(amsterdamdemand,1)
    o=strcmp(nodes,amsterdamdemand{i,1});
    d=strcmp(nodes,amsterdamdemand{i,2});
    if sum(o+d)==2
        odmat(o,d)=odmat(o,d)+amsterdamdemand{i,3}/2;%/\/\/\/\\/\/\/\/\//\
    end
end

for i=1:numel(nodes)
   odmat(i,i)=0; 
end

%% Lines

%Gets a list of lines in cell array "lines"
%Input used is cell array of lines and their respective stops in sequence

% global lines;

linescoded=zeros(sum(numlines),max(maxnodes));
for n=1:numnetworks
    if n==1
        lineshift=0;
    else
        lineshift=sum(numlines(1:n-1));
    end
    for i=1:numlines(n)
        for j=1:maxnodes(n)
            if isequal(amsterdam{n}{1}(i,j+1),{''})
                continue;
            end
            linescoded(i+lineshift,j)=...
                find(strcmp(nodes,amsterdam{n}{1}{i,j+1}));
        end
    end
end
% which lines is a node part of
nodelines=zeros(numnodes,sum(numlines));
for n=1:numnetworks
    if n==1
        lineshift=0;
    else
        lineshift=sum(numlines(1:n-1));
    end
    for i=1:numnodes
        for j=1:numlines(n)
            if any(linescoded(j+lineshift,:)==i)
                nodelines(i,j+lineshift)=1;
            end
        end
    end
end

freq=zeros(1,sum(numlines));
for n=1:numnetworks
    if n==1
        lineshift=0;
    else
        lineshift=sum(numlines(1:n-1));
    end
    for i=1:numlines(n)
        freq(i+lineshift)=amsterdam{n}{3}{i,2};
    end
end
freq(11:20)=3*freq(11:20);%**********************************/\/\/\/\/\/\/\
freq(23:26)=3*freq(23:26);%**********************************/\/\/\/\/\/\/\

%line lengths
for i=1:sum(numlines)
    j=find(linescoded(i,:)==0,1);
    if isempty(j)==1
        j=size(linescoded,2)+1;
    end
    linelength(i)=j-1;
end

%% L' Space

%Gets an adjacency matrix representation of L' space
%Input used is cell array of lines and their respective stops in sequence

% L' space adjacency matrix with representation of 
% different lines as 3rd dimension

ldash=zeros(numnodes,numnodes,sum(numlines)/2);
for n=1:numnetworks
    if n==1
        lineshift=0;
    else
        lineshift=sum(numlines(1:n-1))/2;
    end
    for i=(1:2:numlines(n))
        for j=2:maxnodes(n)
            %origin to destination: ldash(origin,destination,linenum)
            % finds matching node position in "nodes" 
            % (find removed for better speed)
            ldash(strcmp(nodes,amsterdam{n}{1}{i,j}),...
                strcmp(nodes,amsterdam{n}{1}{i,j+1}),((i+1)/2)+...
                lineshift)=amsterdam{n}{2}{i,j};
            ldash(strcmp(nodes,amsterdam{n}{1}{i+1,j}),...
                strcmp(nodes,amsterdam{n}{1}{i+1,j+1}),((i+1)/2)+...
                lineshift)=amsterdam{n}{2}{i+1,j};
        end %uses the order of the stops given in the input
    end
end

% Ldash2 for PSL (ldash uses sum(numlines)/2 system whereas the nature of
% link codes requires sum(numlines) system so making ldash2 for that
ldash2 = zeros(numnodes,numnodes,sum(numlines));
for i = 1:sum(numlines)/2
ldash2(:,:,(i-1)*2+1) = triu(ldash(:,:,i))';
ldash2(:,:,i*2) = tril(ldash(:,:,i));
end

%% P Space

pspace0_freq=zeros(numnodes,numnodes,sum(numlines)/2);
pspaceivt=zeros(numnodes,numnodes);
pspacefreq=zeros(numnodes,numnodes);
pspace0_IVT=zeros(numnodes,numnodes,sum(numlines)/2);

for i=1:numnodes
    connectedlines=mod(find(linescoded==i),sum(numlines));
    connectedlines(connectedlines==0)=sum(numlines);
    for j=1:numel(connectedlines)
        numvisit=0;
        endofthisline=find(linescoded(connectedlines(j),:)==0,1)-1;
        if size(endofthisline,2)==0
            endofthisline=max(maxnodes);
        end
        visited=linescoded(connectedlines(j),...
            find(linescoded(connectedlines(j),:)==i):endofthisline);
        numvisit=numvisit+numel(visited);
        if mod(connectedlines(j),2)~=0
            wholeline=(connectedlines(j)+1)/2;
        else
            wholeline=connectedlines(j)/2;
        end
        visitval=...
            getting_shortest_paths(i,wholeline,visited,numvisit,ldash);
        for k=1:numvisit
            %origin to destination format: (origin,destionation)
            pspace0_freq(i,visited(k),connectedlines(j))=...
                freq(connectedlines(j)); %assigning frequencies
            pspace0_IVT(i,visited(k),connectedlines(j))=...
                visitval(k); %assigning IVT
        end
        pspace0_freq(i,i,connectedlines(j))=0; 
        % cannot travel from j to j, therefore 
        % making those 0 in the adjacency matrix
        pspace0_IVT(i,i,connectedlines(j))=0; 
        % cannot travel from j to j, therefore 
        % making those 0 in the adjacency matrix
        pspacefreq(i,:)=pspacefreq(i,:)+pspace0_freq(i,:,connectedlines(j));
    end
end

for i=1:numnodes %origin
    for j=1:numnodes %destination
        for k=1:sum(numlines)
            pspaceivt(i,j)=pspaceivt(i,j)+...
                pspace0_IVT(i,j,k)*pspace0_freq(i,j,k);
        end
    end
end
pspaceivt=pspaceivt./pspacefreq;
pspaceivt(isnan(pspaceivt))=0;

%% Master Route Choice Set (only topological)
tic
paths=cell(numnodes);

for origin=1:numnodes

    visited=[];
    
    % Rules
    maxtransfers=3;
    
    % Initializations
    nextoriginarray=origin;
    level=0;
    
    routeinformation={};
    oldrouteinformation={};
    
    % Route information for each "last node"
    routepath=[];
    routepath=[routepath origin];
    routedirections=zeros(1,sum(numlines));
    circularnodes = [];%*****
    routeinformation{1}={routepath,routedirections,circularnodes};
    
    % Origin connections
    connections=find(pspaceivt(origin,:));
    
    while level<maxtransfers
        % Increase level
        level=level+1;
        % Initialize origin array
        originarray=nextoriginarray;
        nextoriginarray=[];
        % Add this level nodes to visited
        visited=[visited originarray];
        % Initialize routeinformation
        oldrouteinformation=routeinformation;
        routeinformation={};
        routepath=[];
        
        for i=1:numel(originarray)
            
            % Find all connections to "origin/transfer point"
            connections=find(pspaceivt(originarray(i),:));
            
            % CONDITIONS
            
            % Not same line (destination connection allowed)
            
            if level>1
                clines1=nodelines(connections,:)+...
                    repmat(nodelines(originarray(i),:),...
                    numel(connections),1);
                clines2=nodelines(oldrouteinformation{i}{1}(end-1),:)+...
                    nodelines(originarray(i),:);
                positionsclines2=clines2==2;
%                 find(clines2==2);
                samelinecheck2=clines1==2;
                notsameline=any(bsxfun(@minus,positionsclines2,...
                    samelinecheck2),2);
% sum(clines1(:,positionsclines2)==2,2)...
%                     ~=numel(positionsclines2);
                connections=connections(notsameline);
            end
            
%             connectionlines=nodelines(connections,:);
%             linesvisited=find(oldrouteinformation{i}{2});
%             notonvisitedlines=sum(connectionlines(:,linesvisited),2)==0;
%             connections=connections(notonvisitedlines);

            % No Stop Repititions
            %*****
            visited = oldrouteinformation{i}{1};
            %*****
            visitedconnections=~ismember(connections,visited);
            connections=connections(visitedconnections);
%             routedirections=routedirections(visitedconnections,:);
            
            % No circular paths (to destination)
            if level > 1
                for j = 1:numel(oldrouteinformation{i}{3})
                    connections(connections==...
                        oldrouteinformation{i}{3}(j)) = [];
                end
            end


            % No Reverses (passed to another function)
            [reverselineconnections,routedirections,circularnodes]=...
                getting_noreverses3(originarray(i),connections,...
                oldrouteinformation(i),nodelines,linescoded);
            connections(reverselineconnections)=[];
            routedirections(reverselineconnections,:)=[];
            circularnodes(reverselineconnections) = [];
            
            % No circular paths (conflicting with nodes actually in the path)
            circularwithvisited = [];
            if level > 1 && isempty(circularnodes) == 0
                for j = 1:size(circularnodes,1)
                    for k = 1:numel(visited)
                        if sum(circularnodes{j} == visited(k)) > 0
                            circularwithvisited = ...
                                [circularwithvisited j];
                            break;
                        end
                    end
                end
            end
            connections(circularwithvisited) = [];
            routedirections(circularwithvisited,:) = [];
            circularnodes(circularwithvisited) = [];
                            
            
            % Assign Paths to "Paths"
            for j=1:numel(connections)
                paths{origin,connections(j)}{end+1}=...
                    [oldrouteinformation{i}{1} connections(j)];
            end
            
            % CONDITIONS

            % No Same Line Transfer
%             samelinecheck1=nodelines(connections,:)+...
%                 repmat(nodelines(originarray(i),:),numel(connections),1);
%             samelinecheck2=sum(samelinecheck1~=0,2);
%             notsamelineconnections=...
%                 find(samelinecheck2>nnz(nodelines(originarray(i),:)));
%             connections=...
%                 connections(notsamelineconnections);
%             routedirections=routedirections(notsamelineconnections,:);

            % Find Transfer Connections
            transferconnections=ismember(connections,transferpoints);
            connections=connections(transferconnections);
            routedirections=routedirections(transferconnections,:);
            circularnodes = circularnodes(transferconnections);
            
            
            % Adding route information
            routepath=horzcat(...
                repmat(oldrouteinformation{i}{1},numel(connections),1),...
                connections');
            for j=1:numel(connections)
                routeinformation{end+1}=...
                    {routepath(j,:),routedirections(j,:),circularnodes{j}};
            end
            nextoriginarray=[nextoriginarray connections];
        
        end
    end
end
%}
toc
%% Disrupting links
disrpIVT = zeros(sum(numlines),max(maxnodes));
disrpWT = zeros(sum(numlines),max(maxnodes));
disrptransfers = zeros(sum(numlines),max(maxnodes));
disrpGTC = zeros(sum(numlines),max(maxnodes));
disrpisolating = zeros(sum(numlines),max(maxnodes));
disrpdisconnor = cell(sum(numlines),max(maxnodes));
disrpdisconndest = cell(sum(numlines),max(maxnodes));
disrplinkloads = cell(sum(numlines),max(maxnodes));
disrplinkvc = cell(sum(numlines),max(maxnodes));

% disruptions = [26,1; 1,10; 3,10; 25,4; 10,16; 1,8]; %disruptedline,disruptor
% for disr = 1:size(disruptions,1)
tic   
for disruptedline = 0%disruptions(disr,1)%36:38
    for disruptor = 0%1:linelength(disruptedline)-1%disruptions(disr,2)%
        disruptdest = 0;%disruptor+1;

%% Final route choice set

% Dominance parameters
extratransfer=1; 
extralengthratio = 0.5;
minprob=0.05;

% Generalized cost parameters
alpha = 1; %scaling parameter: trying to see if the spread can be reduced
betaivt=0.04*alpha; %1;%
betawt=0.07*alpha; %3;%
betatransfer=0.334*alpha; %10;%

% Initialization

transferdomroutesnum=zeros(numnodes);
transferdomroutes=cell(numnodes);
linkscoded=cell(numnodes);
pathlines=cell(numnodes);


masterpath=cell(numnodes);
masterlinks=cell(numnodes);
masterlines=cell(numnodes);
masterfreq=cell(numnodes);
masterlength=cell(numnodes);
masterwt=cell(numnodes);
mastertransfers=cell(numnodes);
masterutility=cell(numnodes);
masterprob=cell(numnodes);
masterload=cell(numnodes);
isolatingnum = 0;

isolatedor = [];
isolateddest = [];

for i=1:numnodes
    for j=1:numnodes
%         i
%         j
%         choicesnum(i,j)=numel(paths{i,j});
        if numel(paths{i,j})~=0
%         if choicesnum(i,j)~=0
            
            % Finding number of transfers
            routesizes=cellfun('length',paths{i,j});
            transferdomroutesnum(i,j) = numel(routesizes);
            transferdomroutes{i,j}=...
                paths{i,j}(1:transferdomroutesnum(i,j));
            odroutepath=zeros(transferdomroutesnum(i,j),...
                maxtransfers+extratransfer+1);
            odroutetransfers=zeros(transferdomroutesnum(i,j),1);
            for k=1:transferdomroutesnum(i,j)
                odroutetransfers(k)=numel(transferdomroutes{i,j}{k});
                odroutepath(k,1:odroutetransfers(k))=...
                    transferdomroutes{i,j}{k};
            end
            
            % Finding full paths (in L')
            [odroutelines,odroutelinkscoded,linkscoded,...
                pathlines,NApaths,illogicalpaths]=...
                linking_pspace3sandbox(linkscoded,odroutepath,...
                odroutetransfers,pspace0_IVT,linescoded,...
                disruptedline,disruptor,disruptdest,pathlines);
            
            odroutepath([NApaths illogicalpaths],:) = [];
            odroutetransfers([NApaths illogicalpaths]) = [];
            odroutelines([NApaths illogicalpaths],:) = [];
            odroutelinkscoded([NApaths illogicalpaths],:) = [];
            
            % For isolating links 
            % (links that disconnect some node(s) when disrupted)
            if sum(odroutepath)==0
                isolatingnum = isolatingnum + 1;
                isolatedor = [isolatedor i];
                isolateddest = [isolateddest j];
                continue;
            end
            
            % Filtering transfer dominated routes
            transferfilter = min(odroutetransfers) + 1;
            odroutepath=odroutepath(odroutetransfers<=transferfilter,:);
            odroutelinkscoded = ...
                odroutelinkscoded(odroutetransfers<=transferfilter,:);
            odroutelines = odroutelines(odroutetransfers<=transferfilter,:);
            odroutetransfers=...
                odroutetransfers(odroutetransfers<=transferfilter);

            % Finding IVT and WT
            odroutelength = zeros(size(odroutetransfers,1),1);
            odroutewt = zeros(size(odroutetransfers,1),1);
            for k = 1:size(odroutetransfers,1) %path
                for m = 2:odroutetransfers(k) %transfer node
                    freqsum = 0;
                    avglength = 0;
                    for n = 1:numel(odroutelines{k,m-1}) %alternate lines
                        avglength = avglength + ...
                            pspace0_IVT(odroutepath(k,m-1),...
                            odroutepath(k,m),...
                            odroutelines{k,m-1}(n))*...
                            freq(odroutelines{k,m-1}(n));
                        freqsum = freqsum + freq(odroutelines{k,m-1}(n));
                    end
                    odroutelength(k) = odroutelength(k)+avglength/freqsum;
                    odroutewt(k) = odroutewt(k) + 30./freqsum;
                end
            end
            
            % Filtering length dominated routes
            lengthfilter = min(odroutelength) * (1 + extralengthratio);
            odroutepath=odroutepath(odroutelength<=lengthfilter,:);
            odroutewt=odroutewt(odroutelength<=lengthfilter);
            odroutetransfers=...
                odroutetransfers(odroutelength<=lengthfilter);
            odroutelinkscoded = ...
                odroutelinkscoded(odroutelength<=lengthfilter,:);
            odroutelines = odroutelines(odroutelength<=lengthfilter,:);
            odroutelength=odroutelength(odroutelength<=lengthfilter);
            
            % Applying Path Size Logit
%             [odpsl]=...
%                 psl(odroutepath,odroutetransfers,...
%                 odroutelinkscoded,odroutelength,ldash2);%,...
%                 odrouteutility,odrouteprobability);
%             odrouteutility=odrouteutility+ 1*log(odpsl);
%             odrouteprobability=...
%                 exp(odrouteutility)/sum(exp(odrouteutility));
            
            
            
            % Finding utilities and probabilities and filtering
            % probabilities
            % Changing seconds to minutes (WT is already in minutes)
            odroutelength=odroutelength/60;
            odrouteutility=-betaivt*odroutelength-betawt*odroutewt-...
                betatransfer*(odroutetransfers-2);
            odrouteprobability=...
                exp(odrouteutility)/sum(exp(odrouteutility));
%             odroutepath=odroutepath(odrouteprobability>=minprob,:);
%             odroutelength=odroutelength(odrouteprobability>=minprob);
%             odroutewt=odroutewt(odrouteprobability>=minprob);
%             odroutetransfers=...
%                 odroutetransfers(odrouteprobability>=minprob);
%             odrouteutility=odrouteutility(odrouteprobability>=minprob);
%             odroutelinkscoded = ...
%                 odroutelinkscoded(odrouteprobability>=minprob,:);
%             odroutelines = odroutelines(odrouteprobability>=minprob,:);
%             odrouteprobability=...
%                 exp(odrouteutility)/sum(exp(odrouteutility));
%         
            odroutefreq=odroutelines;
            for f1=1:size(odroutelines,1)
                for f2 = 1:odroutetransfers(f1)-1
                    odroutefreq{f1,f2}=freq(odroutelines{f1,f2});
                end
            end
            
            %% Adding to master choice set
            masterpath{i,j}=odroutepath;
            masterlinks{i,j}=odroutelinkscoded;
            masterlines{i,j}=odroutelines;
            masterfreq{i,j}=odroutefreq;
            masterlength{i,j}=odroutelength;
            masterwt{i,j}=odroutewt;
            mastertransfers{i,j}=odroutetransfers;
            masterutility{i,j}=odrouteutility;
            masterprob{i,j}=odrouteprobability;
        end
    end
end

for i=1:numnodes
    for j=1:numnodes
        if i~=j
            masterload{i,j} = masterprob{i,j}.*odmat(i,j);
        end
    end
end
%}
toc
%% Congestion
tic
congestion10ADA2_2dis2sandbox
toc
%% Statistics

statistics

%% Saving results
disrpIVT(disruptedline,disruptor) = totalIVT;
disrpWT(disruptedline,disruptor) = totalWT;
disrptransfers(disruptedline,disruptor) = ...
    totaltransfers;
disrpGTC(disruptedline,disruptor) = totalutility;
disrpisolating(disruptedline,disruptor) = isolatingnum;
disrpdisconnor{disruptedline,disruptor} = isolatedor;
disrpdisconndest{disruptedline,disruptor} = isolateddest;
disrplinkloads{disruptedline,disruptor} = outputLL(:,end);
disrplinkvc{disruptedline,disruptor} = outputVC(:,end);

%% For spatial effects assessment
% save(['dislines' datestr(datetime('today'),29) '_34_1.mat'],'disrp*');
    end
%     save(['basecase' datestr(datetime('today'),29) '.mat'],'master*','topo*','total*','outputLL','outputVC')
    
end

% end