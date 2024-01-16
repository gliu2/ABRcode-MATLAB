function [d_prime_cache, delta_wave1amp_cache, ste_delta_wave1amp_cache, sortedX_real_cache, sortedX_real_significantchange_cache, delta_wave1amp_significantchange_cache] = get_wave1amp_change_singlemouse(finalTable_metadata, is_mouse, is_side, isfreq, metric_selector, timepoints)
% Calculate D' and wave 1 amplitude change and standard error for an
% individual mouse using single trial wave 1 amplitude distributions.
%
% Input: 
%       ***
%
% Output: 
%       ***cache has index for each time point (pre, post24h, post1w)
%
% Dependencies: none
%
% Last edit George Liu 9/21/2022

n_timepoints = 3;
NUM_CHANGE_METRICS = 2;
sortedX_real_cache = [];
d_prime_cache = cell(n_timepoints, 1);
delta_wave1amp_cache = cell(n_timepoints, 1);
ste_delta_wave1amp_cache = cell(n_timepoints, 1);
sortedX_real_significantchange_cache = cell(n_timepoints, 1);
delta_wave1amp_significantchange_cache = cell(n_timepoints, 1);

for i=1:NUM_CHANGE_METRICS % Plot D' then wave 1 amp change
    for ii = 1:n_timepoints
        combined_filter = is_mouse & is_side & isfreq & metric_selector & timepoints{ii};
        data_selected = finalTable_metadata(combined_filter, :);
        x = finalTable_metadata(1, COLS_METRICVALS);
        y = data_selected(:, COLS_METRICVALS);

        n_pts = sum(isfreq);
        ymean = nanmean(y{:,:}, 1);

        % sort points in order of ascending x value to avoid criss
        % crossing of plot lines.
%             x_double = cellfun(@str2double, x); % convert cell array of char vectors to double matrix
        x_double = x.(1); % convert table 1x1 to double matrix
        [sortedX, sortIndex] = sort(x_double);
        sortedYmean = ymean(sortIndex);

        % Remove nan values to plot uninterrupted curve
        idx = ~(isnan(sortedYmean));
        sortedX_real = sortedX(idx);

        % Calculate change in wave 1 amplitude from baseline to
        % post-noise date
        d_prime = zeros(num_levels, 1);
        delta_wave1amp = zeros(num_levels, 1);
        ste_delta_wave1amp = zeros(num_levels, 1);
        
        num_levels = numel(distribution);

        if ii==1
            distribution_pre = data_selected.('Distribution'){1};
        else                        
            for aa = 1:num_levels % ascending levels
                ind = num_levels - aa + 1;
                distribution_post = data_selected.('Distribution'){1};

                this_pre_wave1amp_dist = distribution_pre{ind};
                this_post_wave1amp_dist = distribution_post{ind};

                % Calculate D' for this stimulus level
                n1 = numel(this_pre_wave1amp_dist);
                mean1 = mean(this_pre_wave1amp_dist);
                std1 = std(this_pre_wave1amp_dist);
                n2 = numel(this_post_wave1amp_dist);
                mean2 = mean(this_post_wave1amp_dist);
                std2 = std(this_post_wave1amp_dist);
                [npool,meanpool,stdpool] = pooledmeanstd(n1,mean1,std1,n2,mean2,std2);
                d_prime(aa) = (mean2 - mean1)/stdpool;

                delta_wave1amp(aa) = mean2 - mean1;
                ste1 = std1/sqrt(n1);
                ste2 = std2/sqrt(n2);
                ste_delta_wave1amp(aa) = sqrt(ste1^2 + ste2^2);
            end

        end

        % Remove stimulus levels that are not divisible by 10
        isdiv10 = mod(sortedX_real, 10)==0;
        sortedX_real_isdiv10 = sortedX_real(isdiv10);
        delta_wave1amp_isdiv10 = delta_wave1amp(isdiv10);
        ste_delta_wave1amp_isdiv10 = ste_delta_wave1amp(isdiv10);

        if i==1 % Plot D'
%             plot(sortedX_real(isdiv10), d_prime(isdiv10), colors(ii), 'LineWidth', LINEWIDTH)
            
            sortedX_real_cache = sortedX_real(isdiv10);
            d_prime_cache{ii} = d_prime(isdiv10);
        else

%             % Plot change in wave 1 amplitude from baseline to after
%             % noise exposure
%             errorbar(sortedX_real_isdiv10, delta_wave1amp_isdiv10, ste_delta_wave1amp_isdiv10, colors(ii), 'LineWidth', LINEWIDTH)

            delta_wave1amp_cache{ii} = delta_wave1amp_isdiv10;
            ste_delta_wave1amp_cache{ii} = ste_delta_wave1amp_isdiv10;
            
            % Check statistical significance from zero
            z_vector = delta_wave1amp_isdiv10./ste_delta_wave1amp_isdiv10;
            p_twoTailed = 2*(1 - normcdf(abs(z_vector)));
            is_significant = p_twoTailed < P_CRIT;
            zz = find(is_significant);
%             text(sortedX_real_isdiv10(zz), delta_wave1amp_isdiv10(zz), '*', ...
%                 'Color', colors(ii), 'FontSize', FONTSIZE*2, 'FontWeight', 'bold') 
            
            sortedX_real_significantchange_cache{ii} = sortedX_real_isdiv10(zz);
            delta_wave1amp_significantchange_cache{ii} = delta_wave1amp_isdiv10(zz);
        end
    end
end
