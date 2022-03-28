function writeMPCSamples(n, matDir)

    outputDir = "samples/samplesMPC_N"+n+".txt";

    disp(['Guardando: '+ outputDir]);
    load(matDir);

    I = size(Mx,1);
    nSamples = size(u,2);

    fileID = fopen(outputDir,'w');

    fprintf(fileID,"N:\t%d\tI:\t%d\tnSamples:\t%d\n", n, I, nSamples);
    fprintf(fileID,'%.30f\t',reshape(xmin',1,[]));
    fprintf(fileID,'\n'); 
    fprintf(fileID,'%.30f\t',reshape(xmax',1,[]));
    fprintf(fileID,'\n'); 
    fprintf(fileID,'%.30f\t',reshape(umin',1,[]));
    fprintf(fileID,'\n'); 
    fprintf(fileID,'%.30f\t',reshape(umax',1,[]));
    fprintf(fileID,'\n');
    fprintf(fileID,'%.30f\t',reshape(A',1,[]));
    fprintf(fileID,'\n');
    fprintf(fileID,'%.30f\t',reshape(B',1,[]));
    fprintf(fileID,'\n');
    fprintf(fileID,'%.30f\t',reshape(Acal',1,[]));
    fprintf(fileID,'\n'); 
    fprintf(fileID,'%.30f\t',reshape(AcalOmgOcal',1,[]));
    fprintf(fileID,'\n'); 
    fprintf(fileID,'%.30f\t',reshape(H',1,[]));
    fprintf(fileID,'\n');    
    fprintf(fileID,'%.30f\t',reshape(Mx',1,[]));
    fprintf(fileID,'\n'); 
    fprintf(fileID,'%.30f\t',reshape(L_invLast',1,[]));
    fprintf(fileID,'\n');

    for sample = 1:nSamples
        fprintf(fileID,'%.30f\t',reshape(x(:,sample)',1,[]));
        fprintf(fileID,'\n');
        fprintf(fileID,'%.30f\t',reshape(yref(:,sample)',1,[]));
        fprintf(fileID,'\n');
        fprintf(fileID,'%.30f\t',reshape(u(:,sample),1,[]));
        fprintf(fileID,'\n');
    end
    fclose(fileID);
end


