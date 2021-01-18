function p_labels=getting_shortest_paths(a,b,visited,numvisit,ldash)
p_labels=zeros(1,numvisit-1);
p_labels(1)=0;

for i=2:numvisit
    if ldash(a,visited(i),b)==0
        %origin to destination format: ldash(origin,destination,linenum)
        p_labels(i)=p_labels(i-1)+ldash(visited(i-1),visited(i),b);
    else
        p_labels(i)=ldash(a,visited(i),b);
    end
end
end