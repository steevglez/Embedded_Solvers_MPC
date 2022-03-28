function [xinfy,uinfy] = stationary_state_values_rt(A,B,C,d,yref)
    n=size(A,1);
    m=size(B,2);
    p=size(C,1);
    L=[A-eye(n),B;C,zeros(p,m)];              
    bl=[-B*d;yref];
    infy=L\bl; 
    xinfy=infy(1:n);         
    uinfy=infy(n+1:end);
end
