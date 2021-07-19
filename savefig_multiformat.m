function savefig_multiformat(figh, path, name)
%SAVEFIG_MULTIFORMAT - Save figure in multiple formats (FIG, SVG, PNG)
%   Inputs
%    figh: figure handle
%    path: save path
%    name: save name
%
%   7/15/21 George Liu
EXT = {'fig'; 'png'; 'svg'};

fig_path = fullfile(path, name);
for j = 1:length(EXT)
    saveas(figh, fig_path, EXT{j})
end

end

