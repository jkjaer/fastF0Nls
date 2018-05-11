function costFunctions = computeAllCostFunctions(x, L, fullPitchGrid, ...
        fftShiftVector, crossCorrelationVectors, Gamma1, Gamma2, ...
        dcIsIncluded)
% compute Zc'*y for all model orders and frequencies where Zc is the
% complex-valued sinusoidal matrix
[harmonicDfts, pitchGrids] = ...
        computeComplexHarmonicDfts(x, fullPitchGrid, L, fftShiftVector);
nPitches = length(fullPitchGrid);
costFunctions = nan(L, nPitches);

for l = 1:L
    % Update the valid pitches and the DFT matrices
    nPitches = sum(~isnan(pitchGrids(l,:)));
   
    if dcIsIncluded
        gamma1 = Gamma1{l+1};
    else
        gamma1 = Gamma1{l};
    end
    gamma2 = Gamma2{l};
                
    if l == 1
        if dcIsIncluded
            lsSol1 = sum(x)*Gamma1{1};
            R1 = computeRowsOfToeplitzHankelMatrix(2, 2, ...
            	crossCorrelationVectors(1:2*l+1, 1:nPitches), ... 
                true, true);
            
            lambda1 = (ones(2, 1)*(real(harmonicDfts(l,1:nPitches))-...
            	sum(R1(1:end-1, :).*lsSol1(:, 1:nPitches),1)));

            lsSol1 = [lsSol1; zeros(1, nPitches)]+lambda1.*gamma1;
        else
            lsSol1 = real(harmonicDfts(l,1:nPitches)).*gamma1;
        end
        
        lsSol2 = imag(harmonicDfts(l,1:nPitches)).*gamma2; 
    else

        if dcIsIncluded
            ll = l+1;
        else
            ll = l;
        end

        R1 = computeRowsOfToeplitzHankelMatrix(ll, ll, ...
        	crossCorrelationVectors(1:2*l+1, 1:nPitches), ... 
            true, dcIsIncluded);
                    
        lambda1 = (ones(ll, 1)*(real(harmonicDfts(l,1:nPitches))-...
            sum(R1(1:end-1, :).*lsSol1(:, 1:nPitches),1)));
                
        lsSol1 = [lsSol1(:, 1:nPitches); zeros(1, nPitches)] ...
                 + lambda1.*gamma1;

        R2 = computeRowsOfToeplitzHankelMatrix(l, l, ...
        	crossCorrelationVectors(1:2*l+1, 1:nPitches), ... 
            false, dcIsIncluded);
        
        lambda2 = (ones(l, 1)*(imag(harmonicDfts(l,1:nPitches))-...
        	sum(R2(1:end-1, :).*lsSol2(:, 1:nPitches),1)));

        lsSol2 = [lsSol2(:, 1:nPitches); zeros(1, nPitches)] ...
                 + lambda2.*gamma2;
                
    end

    % compute the cost function
    if dcIsIncluded
        costFunctions(l, 1:nPitches) = sum(x)*lsSol1(1, :)+...
            sum(real(harmonicDfts(1:l,1:nPitches)).*lsSol1(2:l+1,:),1)+...
            sum(imag(harmonicDfts(1:l,1:nPitches)).*lsSol2,1);
    else
        costFunctions(l, 1:nPitches) = ...
            sum(real(harmonicDfts(1:l,1:nPitches)).*lsSol1,1)+...
            sum(imag(harmonicDfts(1:l,1:nPitches)).*lsSol2,1);
    end
end