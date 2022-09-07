function G = compute_grid(mgrid,mode)

d = numel(mgrid); % dimension

if strcmp(mode,'spatial')
   T = cellfun(@(m)(0:m-1)'/m, num2cell(mgrid), 'UniformOutput', false);

elseif strcmp(mode,'spectral')
   T = cellfun(@(m)(0:m)', num2cell(mgrid), 'UniformOutput', false);

elseif strcmp(mode,'spectral-sym')
   T = cellfun(@(m)(-m:m)', num2cell(mgrid), 'UniformOutput', false);

else
   error(['Unrecognized type: ', mode]);

end

G = cell(1,d); [G{:}] = ndgrid(T{:});
G = cat(d+1,G{:}); % switch back to matrix storage. TODO: check for arbitrary dim

end
