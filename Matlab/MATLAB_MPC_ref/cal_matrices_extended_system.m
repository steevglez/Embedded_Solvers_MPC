function [Acal,AN,Ocal,CN] = cal_matrices_extended_system(A,B,N)
    Acal=[]; Ocal=[];
    for i=1:1:N
        AB=[];
        Acal=[Acal;A^i];
        for j=1:1:i; AB=[AB,A^(i-j)*B]; end
        Fl=[AB,zeros(size(AB,1),(N-i)*size(B,2))];
        Ocal=[Ocal;Fl];
        if i==N; AN=A^i; CN=Fl;  end
    end
end
