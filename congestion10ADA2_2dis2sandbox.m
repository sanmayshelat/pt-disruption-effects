%% Initializations and definitions
    %% Capacity definitions
% Assuming standing capacity is equal to seating capacity
metroseat=150;
metrototal=500;
tramseat=60;
tramtotal=180;
busseat=30;
bustotal=90;
% freq=2*freq;
capacitytotal=zeros(sum(numlines),1);
capacitytotal(1:numlines(1),1)=freq(1:numlines(1))'.*metrototal;
capacitytotal(sum(numlines(1:1))+1:sum(numlines(1:2)),1)=...
    freq(sum(numlines(1:1))+1:sum(numlines(1:2)))'.*tramtotal;
% capacitytotal(sum(numlines(1:2))+1:sum(numlines(1:3)))=freq(sum(numlines(1:2))+1:sum(numlines(1:3)))'.*bustotal;

capacityseat=zeros(sum(numlines),1);
capacityseat(1:numlines(1),1)=freq(1:numlines(1))'.*metroseat;
capacityseat(sum(numlines(1:1))+1:sum(numlines(1:2)),1)=...
    freq(sum(numlines(1:1))+1:sum(numlines(1:2)))'.*tramseat;
% capacityseat(sum(numlines(1:2))+1:sum(numlines(1:3)))=...
%     freq(sum(numlines(1:2))+1:sum(numlines(1:3)))'.*busseat;
%}
    %% Defining crowding parameters and rule
% From Wardman [v/c;seat factor;stand factor]
crowded = [0.50 0.75 1.00 1.25 1.50 1.75 2.00;...
           0.86 0.95 1.05 1.16 1.27 1.40 1.55;...
           0.00 0.00 1.62 1.79 1.99 2.20 2.44];
% Defining parameter curve
cg1 = polyfit(crowded(1,:),crowded(2,:),3);
cg2 = polyfit(crowded(1,3:7),crowded(3,3:7),3);
% Defining crowding functions for IVT based on crowding factors
g1 = @(load) (load<=2).*(cg1(1)*(load).^3+cg1(2)*(load).^2+cg1(3)*(load)+cg1(4))+(load>2).*1.55;
g2 = @(load) 0.*(load<1)+(load>=1).*(load<=2).*(cg2(1)*load.^3+cg2(2)*load.^2+cg2(3)*load+cg2(4))+(load>2).*2.44;
%% First Loading
[linkload,fractionboarding]=...
    linkloading8ADA2_2(masterload,numnodes,numlines,masterfreq...
    ,masterlines,masterpath,mastertransfers,masterlinks,linescoded,...
    linelength,capacitytotal);
%masterprob,odmat
%     outputLL1 = reshape(linkload',[],1);
%     outputLL2 = [];
%     for i = 1:size(linkload,1)
%         outputLL2 = [outputLL2 (size(linkload,2)*(i-1) + linelength(i)):...
%                     (size(linkload,2)*i)];
%     end
%     outputLL1(outputLL2) = [];
%     outputLL=zeros(size(outputLL1,1),1);
%     outputLL(:,1) = outputLL1;
%     vclink=linkload./repmat(capacityseat,1,size(linkload,2));
%     outputLL1 = reshape(linkload',[],1);
%     outputLL2 = [];
%     outputVC1 = reshape(vclink,[],1);
%     outputVC2 = [];
%     for i = 1:size(linkload,1)
%         outputLL2 = [outputLL2 (size(linkload,2)*(i-1) + linelength(i)):...
%             (size(linkload,2)*i)];
%         outputVC2 = [outputVC2 (size(vclink,2)*(i-1) + linelength(i)):...
%             (size(vclink,2)*i)];
%     end
%     outputVC1(outputVC2) = [];
%     outputVC(:,iterationnum) = outputVC1;
%     outputVC1(outputVC2) = [];
%     outputVC(:,iterationnum) = outputVC1;
%% Deterministic User Equiibrium (DUE)
    %% Initializations
dualitygapdiff=0;
iterationnum=1;%first iteration already performed for non-congested case
relativedg=100;
prevrdg=100;
plotrdg=0;
plotdgdiff=0;
plotutility=0;
                            plotutillength=0;
                            plotutilwait=0;
                            
crowdlength=cell(numnodes);
crowdwt=cell(numnodes);
crowdutility=cell(numnodes);
crowdprob=cell(numnodes);
crowdload=cell(numnodes);
oldload=masterload;
assignmentutility(1)=0;     
outputWT=zeros(sum(linelength),1);
    %% DUE While Loop
while dualitygapdiff>0.01 || iterationnum==1 %iterationnum<10%
iterationnum=iterationnum+1;
    
        %% V/C Ratios
vclink=linkload./repmat(capacityseat,1,size(linkload,2));
vcldash=zeros(size(ldash));
for i=1:size(linescoded,1)
    for j=1:size(linescoded,2)-1
        if mod(i,2)==1
            vcldash(j,j+1,(i+1)/2)=vclink(i,j);
        else
            vcldash(j+1,j,i/2)=vclink(i,j);
        end
    end
end
crowdseat=g1(vcldash);
crowdstand=g2(vcldash);    
        %% Getting L', P space with Crowd Factors
% Crowding factor for ldash is weighted avg of sitting and standing
crowdldash=ldash.*...
    ((crowdseat.*1+crowdstand.*max(vcldash-1,0))./...
    (vcldash.*(vcldash>=1)+1.*(vcldash<1)));
% crowdldash=ldash.*...
%     ((crowdseat.*1+crowdstand.*min(vcldash-1,1))./...
%     (min(vcldash,2).*(vcldash>=1)+1.*(vcldash<1)));
[~,crowdpspaceivt,~,~]=...
    pspace(numnodes,numlines,linescoded,crowdldash,freq,maxnodes);
        %% Calculating crowded utility, probability
newload=cell(numnodes);
% WT 
maxwait = ceil(1./fractionboarding);
wtwithdenied = zeros(size(fractionboarding));

for i=1:sum(numlines)
    originalwt = 30./freq(i);
    for j=1:linelength(i)-1
        wtwithdenied(i,j)=...
            [[0 originalwt*2*[1:(maxwait(i,j)-1)]]+...
            originalwt*ones(1,maxwait(i,j))]*...
            [fractionboarding(i,j)*ones(1,(maxwait(i,j)-1)) ...
            1-(maxwait(i,j)-1)*fractionboarding(i,j)]';
    end
end

for i=1:size(linescoded,1)
    if i==1
        k=0;
    else
        k=sum(linelength(1:i-1));
    end
    for j=1:linelength(i)-1
        outputWT(k+j,iterationnum-1)=wtwithdenied(i,j);
    end
end

wtdeniedeasysearch=zeros(numnodes,sum(numlines)); % for easy search by node number later on in the program
for i=1:sum(numlines)
    for j=1:linelength(i)-1
        wtdeniedeasysearch(linescoded(i,j),i)=wtwithdenied(i,j);
    end
end
%*********
for i=1:numnodes %origin
for j=1:numnodes % destination
if i~=j
crowdlength{i,j}=zeros(numel(masterlength{i,j}),1);
crowdwt{i,j}=zeros(numel(masterlength{i,j}),1);
if isempty(masterlength{i,j})~=1
for k=1:numel(masterlength{i,j}) % paths
for m=2:mastertransfers{i,j}(k) % transfernodes
pointa=masterpath{i,j}(k,m-1); % origin transfer node
pointb=masterpath{i,j}(k,m); % dest transfer node


            %% In-vehicle time 
%(convert s to min)
crowdlength{i,j}(k)=crowdlength{i,j}(k)+pspaceivt(pointa,pointb)/60;%crowdpspaceivt(pointa,pointb)/60; %
            %% Waiting time
wtalts=zeros(numel(masterlinks{i,j}{k,m-1}),1);
for alts=1:numel(masterlinks{i,j}{k,m-1})
altsline=masterlines{i,j}{k,m-1}(alts);
% originalwt=30./freq(altsline);
% load on this alt
% if (1-fractionboarding(pointa,altsline))>0
%     maxwait=ceil(1/fractionboarding(pointa,altsline))-1;
%     wtalts(alts)=((originalwt+((2*originalwt).*[0:maxwait]))*...
%         [ones(1,maxwait) mod(1/fractionboarding(pointa,altsline),1)]')/(1/fractionboarding(pointa,altsline));
% if wtwithdenied(pointa,altsline)>0
%     wtalts(alts)=wtwithdenied(pointa,altsline);
% else
%     wtalts(alts)=originalwt;
% end
    wtalts(alts)=wtdeniedeasysearch(pointa,altsline);
end
% Average waiting time is 1/sum(frequencies)
crowdwt{i,j}(k)=crowdwt{i,j}(k)+1/sum(1./wtalts);
                        
end
end
                
            %% Utility, Probability, Load with MSA
crowdutility{i,j}=-betaivt*crowdlength{i,j}-betawt*crowdwt{i,j}-betatransfer*(mastertransfers{i,j}-2);
crowdprob{i,j}=exp(crowdutility{i,j})/sum(exp(crowdutility{i,j}));
crowdload{i,j}=crowdprob{i,j}.*odmat(i,j);
newload{i,j}=oldload{i,j}+...
    ((1/(iterationnum)))*(crowdload{i,j} - oldload{i,j});

                
end
end
end
end

        %% Finding new linkload from crowded utilities
[linkload,fractionboarding]=linkloading8ADA2_2...
(newload,numnodes,numlines,...
masterfreq,masterlines,masterpath,mastertransfers,masterlinks,...
linescoded,linelength,capacitytotal);
        %% Calculating duality gap
% linkload=oldlinkload+(1/iterationnum)*(newlinkload-oldlinkload);
dgdenom=0;
dgnumer=0;
dgnumer1=0;
            dglength=0;
            dgwait=0;
for i=1:numnodes
for j=1:numnodes
if i~=j && isempty(crowdutility{i,j})~=1
dgdenom=dgdenom+max(crowdutility{i,j})*odmat(i,j);
dgnumer=dgnumer+sum(crowdutility{i,j}.*oldload{i,j}); %utility gained by old loading
% dgnumer=dgnumer+sum(crowdutility{i,j}.*newload{i,j});
dglength=dglength+sum((-betaivt*crowdlength{i,j}).*oldload{i,j});
dgwait=dgwait+sum((-betawt*crowdwt{i,j}).*oldload{i,j});
end
end
end
oldload=newload;
oldwt=crowdwt;
oldlength=crowdlength;
oldutility=crowdutility;
        %% Stopping criterion calculation
% changeload = linkload - oldlinkload;
prevrdg=relativedg;
relativedg=dgnumer/dgdenom;
dualitygapdiff=abs(prevrdg-relativedg);
plotutility(iterationnum-1)=dgnumer;
plotutillength(iterationnum-1)=dglength;
plotutilwait(iterationnum-1)=dgwait;
plotrdg(iterationnum-1)=relativedg;
plotdgdiff(iterationnum-1)=dualitygapdiff;
% plot(plotdgdiff);
% plot(plotutility);

        %% Output
vclink=linkload./repmat(capacityseat,1,size(linkload,2));
outputLL1 = reshape(linkload',[],1);
outputLL2 = [];
outputVC1 = reshape(vclink',[],1);
outputVC2 = [];
for i = 1:size(linkload,1)
    outputLL2 = [outputLL2 (size(linkload,2)*(i-1) + linelength(i)):...
                (size(linkload,2)*i)];
    outputVC2 = [outputVC2 (size(vclink,2)*(i-1) + linelength(i)):...
                (size(vclink,2)*i)];
end
outputLL1(outputVC2) = [];
outputLL(:,iterationnum) = outputLL1;
outputVC1(outputVC2) = [];
outputVC(:,iterationnum) = outputVC1;
% pause;    
end
