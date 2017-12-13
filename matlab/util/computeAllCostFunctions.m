function costs = computeAllCostFunctions(x, F, L, pitchBounds, ...
                                         validFftIndices, ...
                                         fftShiftVector,...
                                         nPitches,...
                                         crossCorrelationVectors,...
                                         Gamma1, Gamma2, ...
                                         dcIsIncluded)

X  = fft(x, F);
dftData = X(1:ceil(F/2)).*fftShiftVector;
    
costs = nan(L, nPitches);


for l = 1:L
    % Update the valid pitches and the DFT matrices
    maxFftIndex = ...
        floor(min(F*pitchBounds(2),F/(2*l)-1));
    nPitches = nPitches-(validFftIndices(end)-maxFftIndex);
    validPitchIndices = 1:nPitches;
    validFftIndices = validFftIndices(validPitchIndices);
    
    if dcIsIncluded
        gamma1 = Gamma1{l+1};
    else
        gamma1 = Gamma1{l};
    end
    gamma2 = Gamma2{l};
                
    dftData1 = real(dftData(validFftIndices*l + 1))';
    dftData2 = imag(dftData(validFftIndices*l + 1))';

    if l == 1
        if dcIsIncluded
            lsSol1 = sum(x)*Gamma1{1};
            R1 = computeRowsOfToeplitzHankelMatrix(2, 2, ...
                                                   crossCorrelationVectors(1:2*l+1, ...
                                                              validPitchIndices), ... 
                                                   true, true);
            
            lambda1 = (ones(2, 1)*(dftData1-...
                                   sum(R1(1:end-1, :).*lsSol1(:, validPitchIndices),1)));

            lsSol1 = [lsSol1; zeros(1, nPitches)]+lambda1.*gamma1;
        else
            lsSol1 = dftData1.*gamma1;
        end
        
        lsSol2 = dftData2.*gamma2; 
    else

        if dcIsIncluded
            ll = l+1;
        else
            ll = l;
        end

        R1 = computeRowsOfToeplitzHankelMatrix(ll, ll, ...
                                               crossCorrelationVectors(1:2*l+1, ...
                                                          validPitchIndices), ... 
                                               true, dcIsIncluded);
                    
        lambda1 = (ones(ll, 1)*(dftData1-...
                               sum(R1(1:end-1, :).*lsSol1(:, validPitchIndices),1)));
                
        lsSol1 = [lsSol1(:, validPitchIndices); zeros(1, nPitches)] ...
                 + lambda1.*gamma1;

        R2 = computeRowsOfToeplitzHankelMatrix(l, l, ...
                                               crossCorrelationVectors(1:2*l+1, ...
                                                          validPitchIndices), ... 
                                               false, dcIsIncluded);
        
        lambda2 = (ones(l, 1)*(dftData2-...
                               sum(R2(1:end-1, :).*lsSol2(:, validPitchIndices),1)));

        lsSol2 = [lsSol2(:, validPitchIndices); zeros(1, nPitches)] ...
                 + lambda2.*gamma2;
                
    end

    % compute the cost function
    if dcIsIncluded
        costFunction = sum(x)*lsSol1(1, :);
        for i=1:l
            costFunction = costFunction + ...
                real(dftData(validFftIndices*i+1))'.*lsSol1(i+1, :)...
                + imag(dftData(validFftIndices*i+1))'.*lsSol2(i, :);
        end        
    else
        costFunction = zeros(1, length(validPitchIndices));
        for i=1:l
            costFunction = costFunction + ...
                real(dftData(validFftIndices*i+1))'.*lsSol1(i, :)...
                + imag(dftData(validFftIndices*i+1))'.*lsSol2(i, :);
        end        

    end
    
    costs(l, validPitchIndices) = costFunction; 
end