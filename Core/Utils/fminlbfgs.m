function [x,fval,exitflag,output,grad]=fminlbfgs(funfcn,x_init,optim)
%FMINLBFGS finds a local minimum of a function of several variables. 
%   This optimizer is developed for image registration methods with large 
%	amounts of unknown variables.
%
%   Optimization methods supported:
%	- Quasi Newton Broyden–Fletcher–Goldfarb–Shanno (BFGS)  
%   - Limited memory BFGS (L-BFGS)
%   - Steepest Gradient Descent optimization.
%   
%   [X,FVAL,EXITFLAG,OUTPUT,GRAD] = FMINLBFGS(FUN,X0,OPTIONS) 
%
%   Inputs,
%		FUN: Function handle or string which is minimized, returning an
%				error value and optional the error gradient. 
%		X0: Initial values of unknowns can be a scalar, vector or matrix
%	 (optional)
%		OPTIONS: Structure with optimizer options, made by a struct or
%				optimset. (optimset doesnot support all input options)
%
%   Outputs,
%		X : The found location (values) which minimize the function.
%		FVAL : The minimum found
%		EXITFLAG : Gives value, which explain why the minimizer stopt
%		OUTPUT : Structure with all important ouput values and parameters
%		GRAD : The gradient at this location 
%
%   Extended description of input/ouput variables 
%   OPTIONS,
%		OPTIONS.GoalsExactAchieve : If set to 0, a line search method is
%               used which uses a few function calls to do a good line
%               search. When set to 1 a normal line search method with Wolfe 
%				conditions is used (default).
%		OPTIONS.GradConstr, Set this variable to true if gradient calls are
%				cpu-expensive (default). If false more gradient calls are 
%				used and less function calls.
%	    OPTIONS.HessUpdate : If set to 'bfgs', Broyden–Fletcher–Goldfarb–Shanno 
%				optimization is used (default), when the number of unknowns is 
%				larger then 3000 the function will switch to Limited memory BFGS, 
%				or if you set it to 'lbfgs'. When set to 'steepdesc', steepest 
%				decent optimization is used.
%		OPTIONS.StoreN : Number of itterations used to approximate the Hessian,
%			 	in L-BFGS, 20 is default. A lower value may work better with
%				non smooth functions, because than the Hessian is only valid for
%				a specific position. A higher value is recommend with quadratic equations. 
%		OPTIONS.GradObj : Set to 'on' if gradient available otherwise finited difference
%				is used.
%     	OPTIONS.Display : Level of display. 'off' displays no output; 'plot' displays
%				all linesearch results in figures. 'iter' displays output at  each 
%               iteration; 'final' displays just the final output; 'notify' 
%				displays output only if the function does not converge; 
%	    OPTIONS.TolX : Termination tolerance on x, default 1e-6.
%	    OPTIONS.TolFun : Termination tolerance on the function value, default 1e-6.
%		OPTIONS.MaxIter : Maximum number of iterations allowed, default 400.
% 		OPTIONS.MaxFunEvals : Maximum number of function evaluations allowed, 
%				default 100 times the amount of unknowns.
%		OPTIONS.DiffMaxChange : Maximum stepsize used for finite difference gradients.
%		OPTIONS.DiffMinChange : Minimum stepsize used for finite difference gradients.
%		OPTIONS.OutputFcn : User-defined function that an optimization function calls
%				at each iteration.
%		OPTIONS.rho : Wolfe condition on gradient (c1 on wikipedia), default 0.01.
%		OPTIONS.sigma : Wolfe condition on gradient (c2 on wikipedia), default 0.9. 
%		OPTIONS.tau1 : Bracket expansion if stepsize becomes larger, default 3.
%		OPTIONS.tau2 : Left bracket reduction used in section phase,
%		default 0.1.
%		OPTIONS.tau3 : Right bracket reduction used in section phase, default 0.5.
%   FUN,
%		The speed of this optimizer can be improved by also providing
%   	the gradient at X. Write the FUN function as follows
%   	function [f,g]=FUN(X)
%       	f , value calculation at X;
%   	if ( nargout > 1 )
%       	g , gradient calculation at X;
%   	end
%	EXITFLAG,
%		Possible values of EXITFLAG, and the corresponding exit conditions
%		are
%  		1, 'Change in the objective function value was less than the specified tolerance TolFun.';
%  		2, 'Change in x was smaller than the specified tolerance TolX.'; 
%  		3, 'Magnitude of gradient smaller than the specified tolerance';
%  		4, 'Boundary fminimum reached.';
%  		0, 'Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.FunEvals.';
%  		-1, 'Algorithm was terminated by the output function.';
%  		-2, 'Line search cannot find an acceptable point along the current search';
%
%   Examples
%       options = optimset('GradObj','on');
%       X = fminlbfgs(@myfun,2,options)
%
%   	% where myfun is a MATLAB function such as:
%       function [f,g] = myfun(x)
%       f = sin(x) + 3;
%	    if ( nargout > 1 ), g = cos(x); end
%
%   See also OPTIMSET, FMINSEARCH, FMINBND, FMINCON, FMINUNC, @, INLINE.
%
%   Function is written by D.Kroon University of Twente (Updated Nov. 2010)

% Read Optimalisation Parameters
% defaultopt = struct('Display','final','HessUpdate','bfgs','GoalsExactAchieve',1,'GradConstr',true,  ...
% 			'TolX',1e-6,'TolFun',1e-6,'GradObj','off','MaxIter',400,'MaxFunEvals',100*numel(x_init)-1,  ...
% 			'DiffMaxChange',1e-1,'DiffMinChange',1e-8,'OutputFcn',[], ...
% 			'rho',0.0100,'sigma',0.900,'tau1',3,'tau2', 0.1, 'tau3', 0.5,'StoreN',20);
% 
% f = fieldnames(defaultopt);
% for i=1:length(f)
%     if (~isfield(optim,f{i})||(isempty(optim.(f{i})))), optim.(f{i})=defaultopt.(f{i}); end
% end

% Initialize the data structure
data.fval=0;
data.gradient=0;
data.fOld=NaN; 
numVarTmp = numel(x_init);
data.numberOfVariables = numVarTmp;
data.xInitial = x_init(:);
data.alpha=1;
data.beta=NaN;
data.xOld=x_init(:); 
data.iteration=0;
data.funcCount=0;
data.gradCount=0;
data.exitflag=[];
data.nStored=0;
data.timeTotal=0;
data.timeExtern=0;

data.bracket_exitflag=Inf;
data.section_exitflag=Inf;
data.exitflag=Inf;

nanvec = NaN( numVarTmp, 1 );
data.gradient = nanvec;
data.grad = nanvec;
data.dir = nanvec;
data.gOld = nanvec;

data.fInitial = NaN;
data.fPrimeInitial = NaN;
data.fPrime_alpha = NaN;
data.fPrime_beta = NaN;
data.f_alpha = NaN;
data.f_beta = NaN;
data.fPrime_a = NaN;
data.fPrime_b = NaN;
data.f_a = NaN;
data.f_b = NaN;
data.a = NaN;
data.b = NaN;
data.Hessian = NaN( numVarTmp );
data.TolFunLnS = eps;
data.fminimum = Inf;

data.initialStepLength=1;

coder.varsize( 'data.storefx', 'data.storepx', 'data.storex', [], [ false, true ] );
coder.varsize( 'data.storegx', [], [ true, true ] );

data.storefx=zeros( 1, 0 ); data.storepx=zeros( 1, 0 ); data.storex=zeros( 1, 0 ); data.storegx=zeros( numVarTmp, 0 );

% Switch to L-BFGS in case of more than 3000 unknown variables
if(optim.HessUpdate(1)=='b') 
    if(numVarTmp<3000)
        optim.HessUpdate='bfgs ';
    else
        optim.HessUpdate='lbfgs';
    end
end

data.deltaX=zeros(numVarTmp,optim.StoreN);
data.deltaG=zeros(numVarTmp,optim.StoreN);
data.saveD=zeros(numVarTmp,optim.StoreN);

exitflag=[];

% Display column headers
if(strcmp(optim.Display,'iter'))
    disp('     Iteration  Func-count   Grad-count         f(x)         Step-size');
end

% Calculate the initial error and gradient
[data,fval,grad]=gradient_function(data.xInitial,funfcn, data, optim);
data.gradient=grad;
data.dir = -data.gradient;
data.fInitial = fval;
data.fPrimeInitial= data.gradient'*data.dir(:);
data.fOld=data.fInitial;
data.xOld=data.xInitial;
data.gOld=data.gradient;
    

gNorm = norm(data.gradient,Inf);  % Norm of gradient
data.initialStepLength = min(1/gNorm,5); 

% Show the current iteration
if(strcmp(optim.Display,'iter'))
        fprintf('     %5.0f       %5.0f       %5.0f       %13.6g    \n',data.iteration,data.funcCount,data.gradCount,data.fInitial);
end
  
% Hessian intialization
if(optim.HessUpdate(1)=='b')
	data.Hessian=eye(data.numberOfVariables);
end

% Call output function
if(call_output_function(data,optim,'init')), exitflag=-1; end
    
% Start Minimizing
while(true)
    % Update number of itterations
    data.iteration=data.iteration+1; 

    % Set current lineSearch parameters
    data.TolFunLnS = eps(max(1,abs(data.fInitial )));
    data.fminimum = data.fInitial - 1e16*(1+abs(data.fInitial));
    
	% Make arrays to store linesearch results
    data.storefx=zeros( 1, 0 ); data.storepx=zeros( 1, 0 ); data.storex=zeros( 1, 0 ); data.storegx=zeros( numVarTmp, 0 );
	
    % Find a good step size in the direction of the gradient: Linesearch
    if(optim.GoalsExactAchieve==1)
		data=linesearch(funfcn, data,optim);
    else
        data=linesearch_simple(funfcn, data, optim);
    end
	
    % Check if exitflag is set
    if isfinite(data.exitflag)
        exitflag=data.exitflag;
        data.xInitial=data.xOld; 
        data.fInitial=data.fOld;
        data.gradient=data.gOld;
        break, 
    end
    
    % Update x with the alpha step
    data.xInitial = data.xInitial + data.alpha*data.dir;
    
    % Set the current error and gradient
    data.fInitial =  data.f_alpha;
	data.gradient = data.grad;
    
    % Set initial steplength to 1
    data.initialStepLength = 1;
    
    
    gNorm = norm(data.gradient,Inf);  % Norm of gradient
    
    % Set exit flags 
    if(gNorm <optim.TolFun), exitflag=1; end
    if(max(abs(data.xOld-data.xInitial)) <optim.TolX), exitflag=2; end
    if(data.iteration>=optim.MaxIter), exitflag=0; end
    
    % Check if exitflag is set
    if(~isempty(exitflag))
        break
    end

    % Update the inverse Hessian matrix
    if(optim.HessUpdate(1)~='s')
        % Do the Quasi-Neton Hessian update.
        data = updateQuasiNewtonMatrix_LBFGS(data,optim);
    else
        data.dir = -data.gradient;
    end
  
    % Derivative of direction
    data.fPrimeInitial= data.gradient'*data.dir(:);

    % Call output function
    if(call_output_function(data,optim,'iter')), exitflag=-1; end
    
    % Show the current iteration
    if(strcmp(optim.Display(1),'i')||strcmp(optim.Display(1),'p'))
        fprintf('     %5.0f       %5.0f       %5.0f       %13.6g   %13.6g\n',data.iteration,data.funcCount,data.gradCount,data.fInitial,data.alpha);
    end
    
    % Keep the variables for next iteration
    data.fOld=data.fInitial;
    data.xOld=data.xInitial;
    data.gOld=data.gradient;
end
% Set output parameters
fval=data.fInitial;
grad=data.gradient;
x = data.xInitial;

% Call output function
if(call_output_function(data,optim,'done')), exitflag=-1; end

% Make exist output structure
if(optim.HessUpdate(1)=='b')
    output.algorithm='Broyden–Fletcher–Goldfarb–Shanno (BFGS)';
elseif(optim.HessUpdate(1)=='l')
    output.algorithm='limited memory BFGS (L-BFGS)';
else
    output.algorithm='Steepest Gradient Descent';
end
output.message=getexitmessage(exitflag);
output.iteration = data.iteration;
output.funccount = data.funcCount;
output.fval = data.fInitial;
output.stepsize = data.alpha;
output.directionalderivative = data.fPrimeInitial;
output.gradient = (data.gradient);
output.searchdirection = data.dir;
output.timeTotal=0;    
output.timeExtern=data.timeExtern;
oupput.timeIntern=output.timeTotal-output.timeExtern;
% Display final results
if(~strcmp(optim.Display,'off'))
    fprintf('    Optimizer Results\n');
    fprintf('        Algorithm Used: %s\n', output.algorithm);
    fprintf('        Exit message : %s\n', output.message);
    fprintf('        iterations : %d\n', int32(data.iteration));
    fprintf('        Function Count : %d\n', int32(data.funcCount));
    fprintf('        Minimum found : %f\n', (fval));
    fprintf('        Intern Time : %f seconds\n', (oupput.timeIntern));
    fprintf('        Total Time : %f seconds\n', (output.timeTotal));
end

function message=getexitmessage(exitflag)
    switch(exitflag)
        case 1, message='Change in the objective function value was less than the specified tolerance TolFun.';
        case 2, message='Change in x was smaller than the specified tolerance TolX.'; 
        case 3, message='Magnitude of gradient smaller than the specified tolerance';
        case 4, message='Boundary fminimum reached.';
        case 0, message='Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.FunEvals.';
        case -1, message='Algorithm was terminated by the output function.';
        case -2, message='Line search cannot find an acceptable point along the current search';
        otherwise, message='Undefined exit code';
    end

    
function stopt=call_output_function(data,optim,where)
stopt=false;
if(~isempty(optim.OutputFcn))
    output.iteration = data.iteration;
    output.funccount = data.funcCount;
    output.fval = data.fInitial;
    output.stepsize = data.alpha;
    output.directionalderivative = data.fPrimeInitial;
    output.gradient = (data.gradient);
    output.searchdirection = data.dir;
    stopt=feval(optim.OutputFcn,(data.xInitial),output,where); 
end
        
	
function data=linesearch_simple(funfcn, data, optim)
% Find a bracket of acceptable points
data = bracketingPhase_simple(funfcn, data, optim);

if (data.bracket_exitflag  == 2)
  % BracketingPhase found a bracket containing acceptable points; 
  % now find acceptable point within bracket
  data = sectioningPhase_simple(funfcn, data, optim);
  data.exitflag = data.section_exitflag; 
else
  % Already acceptable point found or MaxFunEvals reached
  data.exitflag = data.bracket_exitflag; 
end

function data = bracketingPhase_simple(funfcn, data,optim)
% Number of itterations
itw=0; 

% Point with smaller value, initial
data.beta=0; 
data.f_beta=data.fInitial; 
data.fPrime_beta=data.fPrimeInitial;

% Initial step is equal to alpha of previous step.
alpha = data.initialStepLength;

% Going up hill
hill=false;

% Search for brackets
while(true)
    % Calculate the error registration gradient
    if(optim.GradConstr)
        [data,f_alpha]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
        fPrime_alpha=nan;
        grad=nan;
    else
        [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
        fPrime_alpha = grad'*data.dir(:);
    end
    
	% Store values linesearch
    data.storefx = [ data.storefx f_alpha ]; 
    data.storepx = [ data.storepx fPrime_alpha ];
    data.storex  = [ data.storex  alpha ];
    data.storegx = [ data.storegx grad(:) ];
    
    % Update step value
    if(data.f_beta<f_alpha)
        % Go to smaller stepsize
        alpha=alpha*optim.tau3;
        
        % Set hill variable
        hill=true;
    else
        % Save current minium point
        data.beta=alpha; data.f_beta=f_alpha; data.fPrime_beta=fPrime_alpha; data.grad=grad;
        if(~hill)
            alpha=alpha*optim.tau1;  
        end
    end
                        
    % Update number of loop iterations
    itw=itw+1; 
		
    if(itw>(log(optim.TolFun)/log(optim.tau3)))
      % No new optium found, linesearch failed.
      data.bracket_exitflag=-2; break; 
    end
    
    if(data.beta>0&&hill)
            % Get the brackets around minimum point
            % Pick bracket A from stored trials
            [~,i]=sort(data.storex,'ascend');
            storefx=data.storefx(i);storepx=data.storepx(i); storex=data.storex(i);
            [~,i]=find(storex>data.beta,1);
            if(isempty(i)), [~,i]=find(storex==data.beta,1); end
            alpha=storex(i(1)); f_alpha=storefx(i(1)); fPrime_alpha=storepx(i(1));
            
            % Pick bracket B from stored trials
            [~,i]=sort(data.storex,'descend');
            storefx=data.storefx(i);storepx=data.storepx(i); storex=data.storex(i);
            [~,i]=find(storex<data.beta,1);
            if(isempty(i)), [~,i]=find(storex==data.beta,1); end
            beta=storex(i(1)); f_beta=storefx(i(1)); fPrime_beta=storepx(i(1));
            
            % Calculate derivatives if not already calculated
            if(optim.GradConstr)
                gstep=data.initialStepLength/1e6; 
                if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
                if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
                [data,f_alpha2]=gradient_function(data.xInitial(:)+(alpha+gstep)*data.dir(:),funfcn, data, optim);
                [data,f_beta2]=gradient_function(data.xInitial(:)+(beta+gstep)*data.dir(:),funfcn, data, optim);
                fPrime_alpha=(f_alpha2-f_alpha)/gstep;
                fPrime_beta=(f_beta2-f_beta)/gstep;
            end

            % Set the brackets A and B
            data.a=alpha; data.f_a=f_alpha; data.fPrime_a=fPrime_alpha;
            data.b=beta; data.f_b=f_beta; data.fPrime_b=fPrime_beta;
  
            % Finished bracketing phase
            data.bracket_exitflag  = 2; return
    end

	% Reached max function evaluations
	if(data.funcCount>=optim.MaxFunEvals), data.bracket_exitflag=0; return; end
end
    

function data = sectioningPhase_simple(funfcn, data, optim)
% Get the brackets
brcktEndpntA=data.a; brcktEndpntB=data.b;

% Calculate minimum between brackets
[alpha,f_alpha_estimated] = pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,data.a,data.b,data.f_a,data.fPrime_a,data.f_b,data.fPrime_b,optim);  
if(isfield(data,'beta')&&(data.f_beta<f_alpha_estimated)), alpha=data.beta; end


[~,i]=find(data.storex==alpha,1);
if((~isempty(i))&&(~isnan(data.storegx(i(1)))))
    f_alpha=data.storefx(i(1)); grad=data.storegx(:,i(1));
else
    % Calculate the error and gradient for the next minimizer itteration
    [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
    if(isfield(data,'beta')&&(data.f_beta<f_alpha)) 
        alpha=data.beta; 
        if((~isempty(i))&&(~isnan(data.storegx(i(1)))))
            f_alpha=data.storefx(i(1)); grad=data.storegx(:,i(1));
        else
            [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
        end
    end
end

% Store values linesearch
data.storefx=[data.storefx f_alpha]; data.storex=[data.storex alpha];

fPrime_alpha = grad'*data.dir(:);
data.alpha=alpha; 
data.fPrime_alpha= fPrime_alpha; 
data.f_alpha= f_alpha;
data.grad=grad;

% Set the exit flag to succes   
data.section_exitflag=Inf;


function data=linesearch(funfcn, data, optim)

% Find a bracket of acceptable points
data = bracketingPhase(funfcn, data,optim);

if (data.bracket_exitflag  == 2)
  % BracketingPhase found a bracket containing acceptable points; 
  % now find acceptable point within bracket
  data = sectioningPhase(funfcn, data, optim);
  data.exitflag = data.section_exitflag; 
else
  % Already acceptable point found or MaxFunEvals reached
  data.exitflag = data.bracket_exitflag; 
end

function data = sectioningPhase(funfcn, data, optim)
%
% sectioningPhase finds an acceptable point alpha within a given bracket [a,b] 
% containing acceptable points. Notice that funcCount counts the total number of 
% function evaluations including those of the bracketing phase. 

while(true)
    
    % Pick alpha in reduced bracket
    brcktEndpntA = data.a + min(optim.tau2,optim.sigma)*(data.b - data.a); 
    brcktEndpntB = data.b - optim.tau3*(data.b - data.a);
    
    % Find global minimizer in bracket [brcktEndpntA,brcktEndpntB] of 3rd-degree 
    % polynomial that interpolates f() and f'() at "a" and at "b".
    alpha = pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,data.a,data.b,data.f_a,data.fPrime_a,data.f_b,data.fPrime_b,optim);  

    % No acceptable point could be found
    if (abs( (alpha - data.a)*data.fPrime_a ) <= data.TolFunLnS), data.section_exitflag = -2; return; end
    
    % Calculate value (and gradient if no extra time cost) of current alpha
    grad = NaN( size( data.grad ) );
    if(~optim.GradConstr)
        [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
        fPrime_alpha = grad'*data.dir(:);
    else
        gstep=data.initialStepLength/1e6; 
        if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
        if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
        [data,f_alpha]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
        [data,f_alpha2]=gradient_function(data.xInitial(:)+(alpha+gstep)*data.dir(:),funfcn, data, optim);
        fPrime_alpha=(f_alpha2-f_alpha)/gstep;
    end

	% Store values linesearch 
	data.storefx=[data.storefx f_alpha]; data.storex=[data.storex alpha]; 
	
    % Store current bracket position of A
    aPrev = data.a; 
    f_aPrev = data.f_a; 
    fPrime_aPrev = data.fPrime_a; 

    % Update the current brackets
    if ((f_alpha > data.fInitial + alpha*optim.rho*data.fPrimeInitial) || (f_alpha >= data.f_a))
        % Update bracket B to current alpha
        data.b = alpha; data.f_b = f_alpha; data.fPrime_b = fPrime_alpha;
    else
        % Wolfe conditions, if true then acceptable point found 
        if (abs(fPrime_alpha) <= -optim.sigma*data.fPrimeInitial) 
            if(optim.GradConstr)
                % Gradient was not yet calculated because of time costs
                [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
                fPrime_alpha = grad'*data.dir(:);
            end
            % Store the found alpha values
            data.alpha=alpha; data.fPrime_alpha= fPrime_alpha; data.f_alpha= f_alpha;
            data.grad=grad;
            data.section_exitflag = Inf; return, 
        end
        
        % Update bracket A
        data.a = alpha; data.f_a = f_alpha;  data.fPrime_a = fPrime_alpha;
        
        if (data.b - data.a)*fPrime_alpha >= 0
            % B becomes old bracket A;
            data.b = aPrev; data.f_b = f_aPrev;  data.fPrime_b = fPrime_aPrev;
        end
    end
    
    % No acceptable point could be found
    if (abs(data.b-data.a) < eps), data.section_exitflag = -2; return, end

    % maxFunEvals reached
    if(data.funcCount >optim.MaxFunEvals), data.section_exitflag = -1; return, end
end

function data = bracketingPhase(funfcn, data, optim)
% bracketingPhase finds a bracket [a,b] that contains acceptable points; a bracket 
% is the same as a closed interval, except that a > b is allowed.
%
% The outputs f_a and fPrime_a are the values of the function and the derivative 
% evaluated at the bracket endpoint 'a'. Similar notation applies to the endpoint 
% 'b'. 

% Parameters of bracket A
data.a = NaN; 
data.f_a = NaN; 
data.fPrime_a = NaN; 

% Parameters of bracket B
data.b = NaN; 
data.f_b = NaN; 
data.fPrime_b = NaN;

% First trial alpha is user-supplied
% f_alpha will contain f(alpha) for all trial points alpha
% fPrime_alpha will contain f'(alpha) for all trial points alpha
alpha = data.initialStepLength;
f_alpha = data.fInitial;              
fPrime_alpha = data.fPrimeInitial;    

% Set maximum value of alpha (determined by fminimum)
alphaMax = (data.fminimum - data.fInitial)/(optim.rho*data.fPrimeInitial); 
alphaPrev = 0;

while(true) 
  % Evaluate f(alpha) and f'(alpha)
  fPrev = f_alpha;
  fPrimePrev = fPrime_alpha;
  
  grad = NaN( size( data.grad ) );
  % Calculate value (and gradient if no extra time cost) of current alpha
  if(~optim.GradConstr)
      [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
      fPrime_alpha = grad'*data.dir(:);
  else
      gstep=data.initialStepLength/1e6;
      if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
      if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
      [data,f_alpha]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
      [data,f_alpha2]=gradient_function(data.xInitial(:)+(alpha+gstep)*data.dir(:),funfcn, data, optim);
      fPrime_alpha=(f_alpha2-f_alpha)/gstep;
  end
  
  % Store values linesearch 
  data.storefx=[data.storefx f_alpha]; data.storex=[data.storex alpha]; 
	
  % Terminate if f < fminimum
  if (f_alpha <= data.fminimum), data.bracket_exitflag = 4; return; end
  
  % Bracket located - case 1 (Wolfe conditions)
  if (f_alpha > (data.fInitial + alpha*optim.rho*data.fPrimeInitial)) || (f_alpha >= fPrev)
    % Set the bracket values
    data.a = alphaPrev; data.f_a = fPrev;  data.fPrime_a = fPrimePrev;
    data.b = alpha; data.f_b = f_alpha;  data.fPrime_b = fPrime_alpha;
    % Finished bracketing phase
    data.bracket_exitflag  = 2; return 
  end

  % Acceptable steplength found
  if (abs(fPrime_alpha) <= -optim.sigma*data.fPrimeInitial) 
      if(optim.GradConstr)
          % Gradient was not yet calculated because of time costs
          [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
          fPrime_alpha = grad'*data.dir(:);
      end
      % Store the found alpha values
      data.alpha=alpha;
      data.fPrime_alpha= fPrime_alpha; data.f_alpha= f_alpha;
      data.grad=grad;
      % Finished bracketing phase, and no need to call sectioning phase
      data.bracket_exitflag = Inf;  return 
  end
  
  % Bracket located - case 2  
  if (fPrime_alpha >= 0)
    % Set the bracket values
    data.a = alpha; data.f_a = f_alpha;  data.fPrime_a = fPrime_alpha;
    data.b = alphaPrev; data.f_b = fPrev; data.fPrime_b = fPrimePrev;
    % Finished bracketing phase
    data.bracket_exitflag  = 2; return
  end
 
  % Update alpha
  if (2*alpha - alphaPrev < alphaMax )
      brcktEndpntA = 2*alpha-alphaPrev; 
      brcktEndpntB = min(alphaMax,alpha+optim.tau1*(alpha-alphaPrev));
      % Find global minimizer in bracket [brcktEndpntA,brcktEndpntB] of 3rd-degree polynomial 
      % that interpolates f() and f'() at alphaPrev and at alpha
      alphaNew = pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,alphaPrev,alpha,fPrev, ...
                                         fPrimePrev,f_alpha,fPrime_alpha,optim);
      alphaPrev = alpha;
      alpha = alphaNew;
  else
      alpha = alphaMax;
  end

  % maxFunEvals reached
  if(data.funcCount >optim.MaxFunEvals), data.bracket_exitflag = -1; return, end
end

function [alpha,f_alpha]= pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,alpha1,alpha2,f1,fPrime1,f2,fPrime2,~)
% finds a global minimizer alpha within the bracket [brcktEndpntA,brcktEndpntB] of the cubic polynomial 
% that interpolates f() and f'() at alpha1 and alpha2. Here f(alpha1) = f1, f'(alpha1) = fPrime1, 
% f(alpha2) = f2, f'(alpha2) = fPrime2.

% determines the coefficients of the cubic polynomial with c(alpha1) = f1, 
% c'(alpha1) = fPrime1, c(alpha2) = f2, c'(alpha2) = fPrime2.
coeff = [(fPrime1+fPrime2)*(alpha2-alpha1)-2*(f2-f1) ...
    3*(f2-f1)-(2*fPrime1+fPrime2)*(alpha2-alpha1) (alpha2-alpha1)*fPrime1 f1];

% Convert bounds to the z-space
lowerBound = (brcktEndpntA - alpha1)/(alpha2 - alpha1);
upperBound = (brcktEndpntB - alpha1)/(alpha2 - alpha1);

% Swap if lowerbound is higher than the upperbound
if (lowerBound  > upperBound), t=upperBound; upperBound=lowerBound; lowerBound=t; end 

% Find minima and maxima from the roots of the derivative of the polynomial.
if all( isfinite( coeff ) )
    sPoints = roots([3*coeff(1) 2*coeff(2) coeff(3)]); 
else
    sPoints = [];
end

% Remove imaginaire and points outside range

sPoints(imag(sPoints)~=0)=[]; 
sPoints(sPoints<lowerBound)=[]; sPoints(sPoints>upperBound)=[];

% Make vector with all possible solutions
sPoints=real([lowerBound sPoints(:)' upperBound]);

% Select the global minimum point
[f_alpha,index]=min(polyval(coeff,sPoints)); z=sPoints(index);

% Add the offset and scale back from [0..1] to the alpha domain
alpha = alpha1 + z*(alpha2 - alpha1);

function [data,fval,grad]=gradient_function(x,funfcn, data, optim)
    % Call the error function for error (and gradient)
    if ( nargout <3 )
        fval=funfcn((x));
        fval = fval( 1 );
        data.timeExtern=data.timeExtern;
        data.funcCount=data.funcCount+1;
    else
        if(strcmp(optim.GradObj,'on'))
            [fval, grad]=funfcn((x));
            fval = fval( 1 );
            data.timeExtern=data.timeExtern;
            data.funcCount=data.funcCount+1;
            data.gradCount=data.gradCount+1;
        else
            % Calculate gradient with forward difference if not provided by the function
            grad=zeros(length(x),1);
            fval=funfcn((x));
            fval = fval( 1 );
            gstep=data.initialStepLength/1e6; 
            if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
            if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
            mag_x = max( 1, abs( x ) );
            for i=1:length(x)
                current_gstep = gstep * mag_x( i );
                x_temp=x; x_temp(i)=x_temp(i)+gstep;
                [fval_g]=funfcn((x_temp));
                fval_g = fval_g( 1 );
                data.funcCount=data.funcCount+1;
                data.timeExtern=data.timeExtern;
                grad(i)=(fval_g-fval)/current_gstep;
            end
        end
        grad=grad(:);
    end
    
function data = updateQuasiNewtonMatrix_LBFGS(data,optim)
% updates the quasi-Newton matrix that approximates the inverse to the Hessian.
% Two methods are support BFGS and L-BFGS, in L-BFGS the hessian is not
% constructed or stored.
% Calculate position, and gradient diference between the
% itterations
deltaX=data.alpha* data.dir;
deltaG=data.gradient-data.gOld;
        
if ((deltaX'*deltaG) >= sqrt(eps)*max( eps,norm(deltaX)*norm(deltaG) ))

    if(optim.HessUpdate(1)=='b')
        % Default BFGS as described by Nocedal
        p_k = 1 / (deltaG'*deltaX);
        Vk = eye(data.numberOfVariables) - p_k*deltaG*deltaX';
        % Set Hessian
        data.Hessian = Vk'*data.Hessian *Vk + p_k * (deltaX*deltaX');
        % Set new Direction
        data.dir = -data.Hessian*data.gradient;
    else
        % L-BFGS with scaling as described by Nocedal
       
        % Update a list with the history of deltaX and deltaG
        data.deltaX(:,2:optim.StoreN)=data.deltaX(:,1:optim.StoreN-1); data.deltaX(:,1)=deltaX;
        data.deltaG(:,2:optim.StoreN)=data.deltaG(:,1:optim.StoreN-1); data.deltaG(:,1)=deltaG;
    
        data.nStored=data.nStored+1; if(data.nStored>optim.StoreN), data.nStored=optim.StoreN; end

        % Initialize variables
        a=zeros(1,data.nStored);
        p=zeros(1,data.nStored);

        q = data.gradient;
        for i=1:data.nStored
            p(i)= 1 / (data.deltaG(:,i)'*data.deltaX(:,i));
            a(i) = p(i)* data.deltaX(:,i)' * q;
            q = q - a(i) * data.deltaG(:,i);
        end
        % Scaling of initial Hessian (identity matrix)
        p_k = data.deltaG(:,1)'*data.deltaX(:,1) / sum(data.deltaG(:,1).^2); 
        
        % Make r = - Hessian * gradient
        r = p_k * q;
        for i=data.nStored:-1:1
            b = p(i) * data.deltaG(:,i)' * r;
            r = r + data.deltaX(:,i)*(a(i)-b);
        end
        
        % Set new direction
        data.dir = -r;
    end
end





    
    