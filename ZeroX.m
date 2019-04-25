% Calculate zero crossing of the signal. 
% INput:  x and y. 
% Output: eposition of zero crossing points
% 
% 2/22/2019
% George Liu
% Dependencies: none

function ZC = ZeroX(x,y)
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);                    % Returns Approximate Zero-Crossing Indices Of Argument Vector
zxidx = zci(y);
ZC = zeros(1, numel(zxidx));
for k1 = 1:numel(zxidx)
    idxrng = max([1 zxidx(k1)-1]):min([zxidx(k1)+1 numel(y)]);
    xrng = x(idxrng);
    yrng = y(idxrng);
    
    % Remove duplicate points to avoid triggering error with interpolation
%     [yrng, index] = unique(yrng); 
    
%     ZC(k1) = interp1( yrng(:), xrng(:), 0, 'linear', 'extrap' );

    if all(size(yrng)==size(unique(yrng)))
        ZC(k1) = interp1( yrng(:), xrng(:), 0, 'linear', 'extrap' );
%         disp(ZC(k1))
%         disp(xrng(1))
    else
        ZC(k1) = mean(xrng);
%         disp(['xrng: ', num2str(xrng(1)), ' ', num2str(xrng(2))])
%         disp(['yrng: ', num2str(yrng(1)), ' ', num2str(yrng(2))])
    end
end
end

% x = linspace(0, 11*pi, 42);                                             % Create Data
% y = sin(x);                                                             % Create Data
% ZC = ZeroX(x,y);
% figure
% plot(x, y, '-r')
% hold on
% plot(ZC, zeros(size(ZC)), 'pg')
% hold off
% grid
