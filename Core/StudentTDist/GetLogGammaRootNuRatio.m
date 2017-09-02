function logGammaRootNuRatio = GetLogGammaRootNuRatio( nu )

    ESTNLSSassert( numel( nu ) == 1, 'ESTNLSS:GetLogGammaRootNuRatio:NuSize', 'GetLogGammaRootNuRatio only supports univariate nu.' );
    ESTNLSSassert( real( nu ) >= 0, 'ESTNLSS:GetLogGammaRootNuRatio:NuSign', 'GetLogGammaRootNuRatio requires nu to be weakly positive.' );

    if real( nu ) >= 10.0522850359735 % fzero(@(nu) eps(cgammaln( 0.5 * ( nu + 1 )))-7e-16,[8,12]) = 10.0522850359735
        inu = 1 ./ nu;
        inu2 = inu .* inu;
        inu4 = inu2 .* inu2;
        inu6 = inu4 .* inu2;
        inu8 = inu4 .* inu4;
        inu10 = inu6 .* inu4;
        inu12 = inu6 .* inu6;
        inu14 = inu8 .* inu6;
        inu16 = inu8 .* inu8;
        
        DSeq = [ -0.25;
                 +0.04166666666666666666667.*inu2;
                 -0.05.*inu4;
                 +0.1517857142857142857143.*inu6;
                 -0.8611111111111111111111.*inu8;
                 +7.852272727272727272727.*inu10;
                 -105.0192307692307692308.*inu12;
                 +1936.602083333333333333.*inu14;
                 -47092.51470588235294118.*inu16;
               ];
        
        logGammaRootNuRatio = -0.2257913526447274323631 + inu .* Shanks( DSeq );
    else
        if isreal( nu )
            logGammaRootNuRatio = gammaln( 0.5 * ( nu + 1 ) ) - gammaln( 1 + 0.5 * nu ) + 0.5 * reallog( nu ) - 0.5723649429247000870717; % 0.5723649429247000870717 = log(Pi)/2
        else
            logGammaRootNuRatio = cgammaln( 0.5 * ( nu + 1 ) ) - cgammaln( 1 + 0.5 * nu ) + 0.5 * log( nu ) - 0.5723649429247000870717; % 0.5723649429247000870717 = log(Pi)/2
        end
    end
    
    ESTNLSSassert( all( ~isnan( logGammaRootNuRatio(:) ) ), 'ESTNLSS:GetLogGammaRootNuRatio:NaNOutputLogGammaRootNuRatio', 'GetLogGammaRootNuRatio returned a NaN output log_y.' );    
    
end
