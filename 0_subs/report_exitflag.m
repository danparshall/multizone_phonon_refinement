function report_exitflag(exitflag);

% http://www.mathworks.com/help/optim/ug/lsqnonlin.html#outputarg_exitflag
	disp(['  exitflag = ', num2str(exitflag)]);
	if exitflag==0
	disp('Number of iterations exceeded options.MaxIterations or number of function evaluations exceeded options.MaxFunctionEvaluations.');
	elseif exitflag == 2
		disp('Change in x was less than the specified tolerance.');
	elseif exitflag == 3
		disp('Change in the residual was less than the specified tolerance.');
	elseif exitflag == 4
		disp('Magnitude of search direction was smaller than the specified tolerance.');
	elseif exitflag == -1
		disp('Output function terminated the algorithm.');
	elseif exitflag == -2
		disp('Problem is infeasible: the bounds lb and ub are inconsistent.');
	else
		assert(exitflag == 1);
		disp('Function converged to a solution x.');
	end

