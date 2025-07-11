function [ EventTimes ] = genNHPP(rate_fh,T,n,dt)
%GENNHPP Generate a NHPP
%   Generates a NHPP, N(t): # of events by time t
% INPUTS
%   rate_fh : function handle for rate (vectorized)
%   T : end of time horizon
%   n : number of sample paths to generate (default n = 1)
% OUTPUTS
%   EventTimes : if n = 1, EventTimes is a vector (length is number of events)
%                if n > 1, EventTimes is a cell array with n rows (# columns is number of events)
%   NumArrived : Number Arrived by time t (time is row, columns are sample paths)

% Input Error Checking ****************************************************
narginchk(2, 4)
if nargin < 3 || isempty(n)
	n = 1; 
end % Default is 1 sample path

if nargin < 4 || isempty(dt)
	dt = 1; 
end % Default is 1 day

% if ~isa(rate_fh,'function_handle')
%     error('rate_fh must be a valid function handle')
% end

if T <= 0 
	error('T must be a positive real number')
end
% End (Input Error Checking) **********************************************



MaxLambda = max(rate_fh);
ph=rate_fh/MaxLambda;
all_t = (0:0.25:300)/1000; %linspace(0, 300, 1200);
all_t = all_t(1:end-1);

% Generate a NHPP 
if n == 1 
	
    % Single Sample Path
    t = 0; Nevents = 0; EventTimes = []; done = false;
    while ~done
        t = t + (-1/MaxLambda)*log(rand());
        if t > T
            done = true;
		else
			% Get ph closest to t
			[~, closest_ind] = min(abs(all_t-t));
            if rand()<= ph(closest_ind)
                Nevents = Nevents+1;
                EventTimes(Nevents) = t;
            end
        end
    end

else 
    % Multiple Sample Paths
    NumPaths = n; EventTimes = {};
    for path = 1:NumPaths
        t = 0; Nevents = 0; done = false; 
        while ~done
            t = t + (-1/MaxLambda)*log(rand());
            if t > T
                done = true;
			else
				[~, closest_ind] = min(abs(all_t-t));
                if rand() <= ph(closest_ind)
                    Nevents = Nevents+1;
                    EventTimes{path,Nevents} = t;
                end
            end
        end
    end
end