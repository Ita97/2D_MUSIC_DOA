function [m,n,p] = spectrum_peaks(A,K)
    % A - spectrum
    % N - number of peaks to select
    [M, N] = size(A);
    p = FastPeakFind(abs(A));
    x = p(1:2:end);
    y = p(2:2:end);
    
    clear p;
    p = zeros(size(x));
    for i = 1:numel(x)
        p(i) = abs(A(y(i),x(i)));
    end
    
    [p, idx]=sort(p,"descend");

    K=min(K,numel(x));
    m=x(idx(1:K));
    n=y(idx(1:K));
    p=p(1:K);

end
