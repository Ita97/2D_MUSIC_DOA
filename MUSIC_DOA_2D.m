function [doa, spectrum] = MUSIC_DOA_2D(x, SensorArray, fc, az_range, el_range, K, sps, fb, plot)
%% INPUT
% x - received signal SxM*N, where S is the signal number of samples, 
% M*N the number of Rx antennas (in x and y), and S the number of snapshots
% SensorArray - phased.NRRectangularPanelArray (default) of size MxN
% fc - carrier frequency
% az_range and el_range - azimuth and elevation range [deg]
% K - number of signals, 1 (default)
% sps - boolean to indicate the use of spatial smoothing or not, false (default)
% if sps=true, then M>=3K/2 and N>=3K/2
% fb - boolean to indicate forward and backward, false (default)
% plot - boolean for plotting, false (default)
%% Output
% doa - 2xK azimuth and elevation of each signals [rad]
% spectrum - 2D MUSIC spectrum [dB]
    

    narginchk(2,9);
    if nargin<9
        plot=false;
        if nargin<8
            fb=false;
            if nargin<7
                sps=false;
                if nargin<6
                    K=1;
                    if nargin<5
                        el_range = -45:45;     
                        if nargin<4
                            az_range = 60:60;    
                        end
                    end
    
                end
            end
        end
    end

    M = SensorArray.Size(1); N = SensorArray.Size(2);
    %d1 = SensorArray.Spacing(1); d2 = SensorArray.Spacing(2);
    S = size(x,1);

    if sps && (M<3*K/2 || N<3*K/2)
        error("Array size too small for DOA estimation with spatial smoothing!")
    end

    x = permute(x,[2 1]);

    %% 2D Spatial Smoothing (SPS)
    if sps
        M1 = M-K+1; N1 = N-K+1;
        Ms = M-M1+1; Ns = N-N1+1;
        R = zeros(M1*N1,M1*N1);
        
        for m = 1:Ms
            for n = 1:Ns
                J = flip(diag(ones(M1*N1,1)));
                Rf = zeros(M1*N1,M1*N1);
                Rb = zeros(M1*N1,M1*N1);

                matrix_filter = zeros(M,N);
                matrix_filter(m:(M1+m-1),n:(N1+n-1))=1;
                matrix_filter=logical(reshape(matrix_filter,[],1));
                y_mn = x(matrix_filter,:);

                if fb
                    for s = 1:S
                        Rf = Rf + y_mn(:,s)*y_mn(:,s)'; % Sub-array Covariance Matrix forward
                        y_b = J*conj(y_mn(:,s));
                        Rb = Rb + y_b*y_b'; % Sub-array Covariance Matrix backward
                    end
                    R = R+(Rf + Rb)/(2*S);
                else
                    for s = 1:S
                        Rf = Rf + y_mn(:,s)*y_mn(:,s)'; % Sub-array Covariance Matrix forward
                    end
                    R = R+Rf/S;
                end
                
            end
        end
        R = R/(Ms*Ns); 

        M=M1; N=N1;
    
    %% No SPS
    else
        R = zeros(M*N,M*N); 
        if fb
            Rf = zeros(M*N,M*N); 
            Rb = zeros(M*N,M*N); 
            J = flip(diag(ones(M*N,1)));
            for s=1:S
                Rf = Rf + x(:,s)*x(:,s)';    % Forward Covariance Matrix
                y = J*conj(x(:,s));
                Rb = Rb + y*y'; % Backward Covariance Matrix
            end
            R = (Rf + Rb)/(2*S); % Forward-Backward Covariance Matrix
        else
            for s=1:S
                R = R + x(:,s)*x(:,s)';    
            end
            R = R/S; % Covariance Matrix
        end
    end

    %% Forward and Backward
    [az,el,P,p] = MUSIC_URA_estimator(R,M,N,az_range,el_range,K,SensorArray,fc);

    doa = [az; el];
    spectrum = pow2db(abs(P));    % Spatial spectrum in dB
    
    if plot
        max_PdB = pow2db(abs(p));    % Peaks of the spatial spectrum in dB
        figure; hold on;
        mesh(az_range,el_range,spectrum,"FaceColor","flat");
        plot3(az,el,max_PdB,'ro')   % DoA peaks
        
        xlim([az_range(1),az_range(end)]); ylim([el_range(1),el_range(end)]);
        title('Custom MUSIC')
        xlabel('Azimuth Angle [deg]')
        ylabel('Elevation Angle [deg]')
        zlabel('Power [dB]')
    end

end


