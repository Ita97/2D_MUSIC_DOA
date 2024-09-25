function [az,el,P,p] = MUSIC_URA_estimator (R,M,N,az_range,el_range,K,AntennaArray,fc)    % MUSIC 2D-DoA in URA MxN 
%% INPUT
% R - covariance matrix, M - antenna elements over x,
% N - antenna elements over y,
% az_range, el_range - Azimuth and Elevation range limits [deg], 
% K - number of signals
%% OUTPUT
% az,el - 2D MUSIC DOA estimation [deg],
% P - 2D MUSIC Spatial Spectrum,
% p - Spatial Spectrum peak power
    
    c = physconst('LightSpeed');
    %AntennaArray.release();
    copied_antenna = AntennaArray.clone();
    copied_antenna.Size = [M N 1 1];
    steervec = phased.SteeringVector("SensorArray",copied_antenna, "PropagationSpeed",c);
    % az_angles = az_range(1):0.1:az_range(2);
    % el_angles = el_range(1):0.1:el_range(2);
    % az_angles = linspace(az_range(1), az_range(2), 720);
    % el_angles = linspace(el_range(1), el_range(2), 720);

    [V,Lam] = eig(R);   % Eigen decomposition
    lambdas = diag(Lam);  % Extract eigenvalues
    [lambdas,idx] = sort(abs(lambdas), 'descend');  % Sort eigenvalues in descending order
    
    % Choose number of sources (K) explicitly, based on prior knowledge
    Nn = length(lambdas) - K;  % Number of noise eigenvectors
    
    % Extract the noise subspace eigenvectors
    Vn = V(:, idx(K+1:end));  % Noise eigenvectors
    VVn = Vn*Vn';
    P = zeros(length(el_range),length(az_range));
    for m = 1:length(az_range)
        for n = 1:length(el_range)
            %angles=[];
            ph = az_range(m); 
            ps = el_range(n);      
            a = steervec(fc,[ph;ps]);  
            %a = a(:);      % Steering vector
            %P(n,m)=1./sum(abs((a'*VVn)).^2,2);
            P(n,m) = 1./(a'*VVn*a);   % Spatial spectrum for MUSIC
        end
    end
    
    [mmax,nmax] = spectrum_peaks(P,K);       
    az = az_range(mmax); el = el_range(nmax);   % Estimation in degrees of the Azimuth / Zenith angles
    p=diag(P(nmax,mmax));
end