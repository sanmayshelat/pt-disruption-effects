function [reverselineconnections,routedirections,circularnodes]=...
    getting_noreverses3(origin,connections,oldrouteinformation,nodelines,linescoded)

routedirections=repmat(oldrouteinformation{1}{2},numel(connections),1);
circularnodes = cell(numel(connections),1);
circularnodes(:) = oldrouteinformation{1}(3);

reverselineconnections=[];
for j=1:numel(connections)
    reverse=0;
    % Find all lines on which both origin and connection exist
    commonlines=find((nodelines(origin,:)+...
        nodelines(connections(j),:))==2);
    
    % Adjacent nodes in commonlines
    adjacentnodes=zeros(numel(commonlines),1);
    for k=1:numel(commonlines)
        direction_multiplier=mod(commonlines(k),2);
        % Find position of origin and connection in the line
        originnumber=find(linescoded(commonlines(k),:)==...
            origin);
        connectionnumber=find(linescoded(commonlines(k),:)==...
            connections(j));
        
        if originnumber<connectionnumber
            if oldrouteinformation{1}{2}(commonlines(k))==-1
                reverselineconnections=[reverselineconnections j];
                reverse=1;
                break;
            else
                routedirections(j,commonlines(k))=1;
                circularnodes{j,1} = [circularnodes{j,1} ...
                    linescoded(commonlines(k),...
                    originnumber+1:connectionnumber-1)];
                
                if direction_multiplier~=0
                    routedirections(j,commonlines(k)+1)=-1;
                else
                    routedirections(j,commonlines(k)-1)=-1;
                end
            end
        else
            continue;
        end
        % Node before connection in L'
%         origin
%         commonlines(k)
%         k
%         j
%         connectionnumber
        adjacentnodes(k)=...
            linescoded(commonlines(k),...
            connectionnumber-routedirections(j,commonlines(k)));
    end
    
    if reverse==1
        continue;
    end
    adjacentnodes(adjacentnodes==0)=[];
%     adjacentnodes=unique(adjacentnodes);
    donethisadjnode=zeros(1,numel(adjacentnodes));
    for k=1:numel(adjacentnodes)
        if sum(donethisadjnode(donethisadjnode==adjacentnodes(k)))==1
            continue;
        end
        donethisadjnode(k)=adjacentnodes(k);
        commonlines2=find((nodelines(connections(j),:)+...
        nodelines(adjacentnodes(k),:))==2);
        %commonlines2=setdiff(commonlines2,commonlines);
        
        for m=1:numel(commonlines2)
            if sum(commonlines(commonlines==commonlines2(m)))>1
                continue;
            end
            direction_multiplier=mod(commonlines2(m),2);
            % Find position of connection and node before connection on
            % this line
            connectionnumber=find(linescoded(commonlines2(m),:)==connections(j));
            adjconnectionnumber=find(linescoded(commonlines2(m),:)==...
                    adjacentnodes(k));
            if adjconnectionnumber<connectionnumber
                if oldrouteinformation{1}{2}(commonlines2(m))==-1
                    reverselineconnections=[reverselineconnections j];
                    break;
                else
                    routedirections(j,commonlines2(m))=1;
                    if direction_multiplier~=0
                        routedirections(j,commonlines2(m)+1)=-1;
                    else
                        routedirections(j,commonlines2(m)-1)=-1;
                    end
                end
            else
                continue;
            end
        end
    end
end
    
end