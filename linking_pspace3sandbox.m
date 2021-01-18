function [odroutelines,odroutelinkscoded,linkscoded,...
    pathlines,NApaths,illogicalpaths]=...
    linking_pspace3sandbox...
    (linkscoded,odroutepath,odroutetransfers,pspace0_IVT,linescoded,...
    disruptedline,disruptor,disruptdest,pathlines)

% Initialisation
odroutelinkscoded=cell(size(odroutepath,1),size(odroutepath,2)-1);
odroutelines=cell(size(odroutepath,1),size(odroutepath,2)-1);
numnodes=size(pspace0_IVT,1);
NApaths = [];
illogicalpaths = [];
% For every path b/w O-D
for i=1:size(odroutepath,1)
%     disruptedpathflag = 0;
    % For every transfer node in the path
    for j=2:odroutetransfers(i)
%         % Disrupted path check
%         if disruptedpathflag == 1
%             odroutelinkscoded{i,j-1} = 0;
%             odroutelines{i,j-1} = 0;
%             break;
%         end
        
        if isempty(linkscoded{odroutepath(i,j-1),odroutepath(i,j)})==1
            % Find connecting lines
            connectinglines=...
                pspace0_IVT(odroutepath(i,j-1),odroutepath(i,j),:);
            connectinglines=connectinglines>0;
            odroutelines{i,j-1}=find(connectinglines);
            connectinglines=linescoded(connectinglines,:);
            
            % Find position of adjacent transfer nodes on line(s)
            originposition=...
                mod(find(connectinglines'==odroutepath(i,j-1)),...
                size(linescoded,2));
            originposition(originposition==0)=size(linescoded,2);
            destposition=...
                mod(find(connectinglines'==odroutepath(i,j)),...
                size(linescoded,2));
            destposition(destposition==0)=size(linescoded,2);
            
            % Disrupted path decision
            disruptedlinepos = odroutelines{i,j-1}==disruptedline;
            %one of the lines is a disrupted lines
            if sum(disruptedlinepos)~=0 
                %problem is only if disrupted part is being used
                if originposition(disruptedlinepos)<disruptdest &&...
                        destposition(disruptedlinepos)>disruptor
                    %if alternative lines are NOT available
                    if numel(odroutelines{i,j-1})==1
                        %0 signifies that connection not available
                        linkscoded{odroutepath(i,j-1),...
                            odroutepath(i,j)} = 0;
                        odroutelinkscoded{i,j-1} = 0;
                        pathlines{odroutepath(i,j-1),...
                            odroutepath(i,j)} = 0;
                        odroutelines{i,j-1} = 0;
%                         disruptedpathflag = 1;
                        NApaths = [NApaths i];
                        break; %this path is not possible
                    %if alternative lines are available
                    else
                        odroutelines{i,j-1}(disruptedlinepos) = [];
                        originposition(disruptedlinepos) = [];
                        destposition(disruptedlinepos) = [];
                    end
                end
            end
            
            % Find coded links
            % For all lines connecting these 2 nodes
            for k=1:size(odroutelines{i,j-1},1) % previously was sum(connectinglines)???
                codedlinks=zeros(destposition(k)-originposition(k),1);
                for m=1:destposition(k)-originposition(k)
                    % for each link, its origin and destination
                    codedlinks(m)=numnodes^2*(odroutelines{i,j-1}(k)-1)+...
                        numnodes*...
                        (connectinglines(k,originposition(k)+m-1)-1)+...
                        connectinglines(k,originposition(k)+m-1+1);
                end
                linkscoded{odroutepath(i,j-1),odroutepath(i,j)}{k}=...
                    codedlinks;
                pathlines{odroutepath(i,j-1),odroutepath(i,j)}{k}=...
                    odroutelines{i,j-1};
            end
            odroutelinkscoded{i,j-1}=...
                linkscoded{odroutepath(i,j-1),odroutepath(i,j)};
            
        else
            % if path disrupted (it is not a cell anymore but integer 0)
            if iscell(linkscoded{odroutepath(i,j-1),odroutepath(i,j)}) == 0
                odroutelinkscoded{i,j-1} = 0;
                odroutelines{i,j-1} = 0;
%                 disruptedpathflag = 1;
                NApaths = [NApaths i];
                break; %this path is not possible
            else
                odroutelinkscoded{i,j-1} = ...
                    linkscoded{odroutepath(i,j-1),odroutepath(i,j)};
%                 connectinglines=...
%                     pspace0_IVT(odroutepath(i,j-1),odroutepath(i,j),:);
%                 connectinglines=connectinglines>0;
%                 odroutelines{i,j-1}=find(connectinglines);
                odroutelines{i,j-1} = ...
                    pathlines{odroutepath(i,j-1),odroutepath(i,j)}{:};
            end
        end
        
        % Illogical paths (same line transfers)
        if j>2
%             i,j
            temp = ismembc(odroutelines{i,j-1},odroutelines{i,j-2});
            l2 = numel(odroutelines{i,j-1});
            l1 = numel(odroutelines{i,j-2});
            if sum(temp) == max(l1,l2) %all lines same
                illogicalpaths = [illogicalpaths i];
                break;
            
%             elseif sum(temp) > 0 && min(l1,l2) == 1 %at least one line same and one of them has no alternative
%                 if l1 == 1%deleting common line from route with alternative
%                     odroutelinkscoded{i,j-1}...
%                         (odroutelines{i,j-1}==odroutelines{i,j-2}) = [];
%                     odroutelines{i,j-1}...
%                         (odroutelines{i,j-1}==odroutelines{i,j-2}) = [];
%                 else
%                     odroutelinkscoded{i,j-2}...
%                         (odroutelines{i,j-2}==odroutelines{i,j-1}) = [];
%                     odroutelines{i,j-2}...
%                         (odroutelines{i,j-2}==odroutelines{i,j-1}) = [];
%                 end
                
            elseif sum(temp) > 0 % at least one line same (but not all same)
                if sum (odroutelines{i,j-1}(temp) == disruptedline) == 0 %if none of overlapping lines is disrupted line
                    illogicalpaths = [illogicalpaths i];
                    break;
                end
            end
        end
        
    end
end
end
