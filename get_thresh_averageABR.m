function thresh = get_thresh_averageABR(M_sorted, A_descending)
% Uses cross covariance between adjacent thresholds. The threshold level is
% where cross correlation between average ABR of that level and level above
% is more than a multiple above cross correlation between the noise and the level above it. 
%
% Input: M_sorted - sorted average ABR traces. m x n matrix, each row is average ABR trace at level, in
% descending intensity level. All same frequency.
%
% output: value of row that is threshold level
%
% Method based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6726521/
%
% 8/16/21 George Liu
% Dependencies: none

%% Constant
CRITERION = 0.35;
TYPICAL_THRESHOLD = 35;
TYPICAL_WIDTH = 35;

%% function

x = M_sorted;

% % RMS
% ave_rms = rms(x, 2);
% ave_rms_max = max(ave_rms, [], 2);

% % Max minus min values
% ave_maxmin = max(x, [], 2) -  min(x, [], 2);

% % dot product with max average ABR
% loudestABR = repmat(x(1,:), size(x,1), 1);
% dot_prod = dot(x, loudestABR, 2);

% dot product with previous ABR normalized
prevABR = x([1, 1:end-1],:);
prevABR_norm = sqrt(dot(prevABR, prevABR, 2));
prevABR_unit = prevABR./prevABR_norm;
prev_dotprod = dot(x, prevABR_unit, 2);
prev_dotprod_norm = dot(x./diag(sqrt(x*x')), prevABR_unit, 2);

% % Cross covariance, normalized
% [ave_xcov, ~] = xcov(x', 'normalized');
% [m, n] = size(ave_xcov);
% row = ceil(m/2);
% cols = sqrt(n);
% ave_xcov_mat = reshape(ave_xcov(row,:), [cols, cols]);
% r = ave_xcov_mat([1, 2:(cols+1):end]);
% r = reshape(r, size(ave_rms_max)); 

% Calculate threshold based on noise floor
% metric = rms([r, ave_rms_max], 2); % RMS of covariance and RMS
% metric = ave_rms_max; % RMS
% metric = dot_prod; %dot product with maximum ABR trace
% metric = prev_dotprod; % dot product with previous ABR unit trace
metric = prev_dotprod_norm; % dot product of normalized trace with previous ABR unit trace. Same as xcov normalized

% Fit to sigmoid function to obtain threshold estimate
ft = fittype('0 + (1-0)/(1 + exp(-4*log(3)*(x-xmid)/xscale80))','indep','x'); % general form: ft = fittype('L + (U-L)/(1 + exp(-4*log(3)*(x-xmid)/xscale80))','indep','x');
mdl = fit(A_descending, metric,ft, 'start', [TYPICAL_THRESHOLD, TYPICAL_WIDTH]);
all_db = -20:1:95;
mdl_out = mdl(all_db);
diff = abs(mdl_out - CRITERION);
[~, i_min] = min(diff);
thresh = all_db(i_min(1));

end
