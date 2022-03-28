function [Acal,Ocal,Omg,H,Mx]=setupMPC(A,B,OmegaN,Omega,Gamma,N)
    % ---Calculos que se tienen que hacer una sola vez --------------------
    [Acal,Ocal] = ABcal(A,B,N);
    Gmm=kron(eye(N),Gamma);
    Omg=blkdiag(kron(eye(N-1),Omega),OmegaN);
    %Gmm = mdiag(Gamma,N);                  
    %Omg = mdiag(Omega,OmegaN,N);
    % Matriz H del funcional de costo VN(x_0,vec{u})
    H=2*Gmm+2*Ocal'*Omg*Ocal;              
    H=(H+H')/2;
    % Matriz M de las restricciones lineales Mu<=c
    M=[Ocal;-Ocal];  
    % Convertir restricciones caja en restricciones de Desigualdad
    Mx=[M;eye(N);-eye(N)];
end

% Función para generar matrices diagonales por bloques
function [M] = mdiag(varargin) 
    if nargin==2
        A=varargin{1}; N=varargin{2}; M=A;    
        for i=1:1:N-1; M=blkdiag(M,A); end
    end
    if nargin==3
        A=varargin{1}; B=varargin{2}; N=varargin{3};  M=A;
        for i=1:1:N-2; M=blkdiag(M,A); end
        M=blkdiag(M,B);
    end
end
