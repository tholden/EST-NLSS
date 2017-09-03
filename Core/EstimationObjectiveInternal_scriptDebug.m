% ESTIMATIONOBJECTIVEINTERNAL_SCRIPT   Generate MEX-function
%  EstimationObjectiveInternal_mex from EstimationObjectiveInternal.
% 
% Script generated from project 'EstimationObjectiveInternal.prj' on 03-Sep-2017.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.EnableMemcpy = true;
cfg.InitFltsAndDblsToZero = true;
cfg.EnableOpenMP = false;
cfg.MATLABSourceComments = true;
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = true;
cfg.ConstantFoldingTimeout = 0;
cfg.CompileTimeRecursionLimit = 2147483647;
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
cfg.SaturateOnIntegerOverflow = false;
cfg.EnableAutoExtrinsicCalls = false;
cfg.InlineThreshold = 0;
cfg.InlineThresholdMax = 0;
cfg.InlineStackLimit = 2147483647;
cfg.StackUsageMax = 16777216;
cfg.IntegrityChecks = true;
cfg.ResponsivenessChecks = true;
cfg.ExtrinsicCalls = false;
cfg.EchoExpressions = false;
cfg.GlobalDataSyncMethod = 'NoSync';
cfg.EnableJIT = true;
cfg.EnableDebugging = true;
cfg.LaunchReport = true;

%% Define argument types for entry-point 'EstimationObjectiveInternal'.
ARGS = cell(1,1);
ARGS{1} = cell(8,1);
ARGS{1}{1} = coder.typeof(0,[Inf Inf],[1 1]);
ARGS{1}{2} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{3} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{4} = coder.typeof(false);
ARGS{1}{5} = coder.typeof(false);
ARGS{1}{6} = coder.typeof(0);
ARGS{1}{7} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{8} = coder.typeof(0,[Inf Inf],[1 1]);

%% Invoke MATLAB Coder.
codegen -config cfg EstimationObjectiveInternal -args ARGS{1} -g -jit -launchreport -O disable:inline

