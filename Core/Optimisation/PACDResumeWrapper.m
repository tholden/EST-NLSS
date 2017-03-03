function [ x, f, PersistentState ] = PACDResumeWrapper( OptiFunction, x, lb, ub, PersistentState, varargin )

    InitialTimeOutLikelihoodEvaluation = Inf;
    
    [ x, f, PersistentState ] = PACDMinimisation( ...
        @( XV, PS, DesiredNumberOfNonTimeouts ) ParallelWrapper( @( X ) OptiFunction( X, PS ), XV, DesiredNumberOfNonTimeouts, InitialTimeOutLikelihoodEvaluation ),...
        x, lb, ub, [], [], PersistentState, true );
    
    x = max( lb, min( ub, x ) );
    
    f = -f;
    
end

