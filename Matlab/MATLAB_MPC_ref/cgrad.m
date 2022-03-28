function [x,iter] = cgrad(A,b,x0,imax,tol)
i=0;
x=x0;
r=b-A*x;
d=r;
dw=r'*r;
tce=norm(r);
iter=i;
while (i<imax && tce>tol)
    q=A*d;
    alpha=dw/(d'*q);
    x=x+alpha*d;
    if floor(i/50)==i/50
        r=b-A*x;
    else
        r=r-alpha*q;
    end
    dl=dw;
    dw=r'*r;
    beta=dw/dl;
    d=r+beta*d;
    i=i+1;
    iter=i;
    tce=norm(x-x0);
    x0=x;
end