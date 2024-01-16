function effectsize = cohen_d_unpooled(mean1, var1, n1, mean2, var2, n2)
% Calculate unpooled effect size (Cohen's D).
%
% Why unpooled std instead of pooled, see comments on https://statisticsbyjim.com/basics/cohens-d/
%
% Last edit: 8/24/23 George Liu
% Dependencies: none

unpooledstd = sqrt(var1/n1 + var2/n2);

effectsize = (mean1  - mean2) / unpooledstd;

end