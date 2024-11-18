function [ip]=myfind(S,str)
n=length(S);
m=length(str);
for j=1:m
    i=1;
    while ~strcmp(S(i),str(j)) & i < n
        i=i+1;
    end;
    if (strcmp(S(i),str(j)))
        %   disp(['Found ' str '=' S(i)]);
        ip(j)=i;
    else
        disp([str(j) ' NOT found']);
        ip(j)=-1;
    end;
end