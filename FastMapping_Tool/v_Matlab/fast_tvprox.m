function [out] = fast_tvprox(inp, lam)

% This function computes the proximity operator for the TV-norm using 
% Condat's direct algorithm
% See also: `Condat, L. (2013). A direct algorithm for 1-D total variation denoising. IEEE Signal Processing Letters, 20(11), 1054?1057. doi:10.1109/LSP.2013.2278339 <http://dx.doi.org/10.1109/LSP.2013.2278339>`_

% Adapted from Condat's original C code. (https://www.gipsa-lab.grenoble-inp.fr/~laurent.condat/download/condat_fast_tv.c) 


% Younes Farouj, @ MIPLAB-EPFL,
% May 30, 2018

%%

out = zeros(size(inp));
N = size(inp,1);
  

        % a)
        k = 1; k0 = 1;
        vmin = inp(1) - lam; vmax = inp(1) + lam;
        umin = lam; umax = - lam;
        kp = 1; km = 1; 
        
        while true
          % #2
            if k == N
                out(N) = vmin + umin;
                break;
            end
            
            
            if k+1 >= N
                break;
            end

          % #3
            if inp(k+1) + umin < vmin - lam
                for i = k0:km
                    out(i) = vmin;
                end
                k = km + 1; k0 = km + 1; kp = km + 1; km = km + 1;
                vmin = inp(k);
                vmax = inp(k) + 2*lam;
                umin = lam;
                umax = - lam;
            else
          % #4
                if inp(k+1) + umax > vmax + lam
                    for i = k0:kp
                        out(i) = vmax;
                    end
                    k = kp + 1; k0 = kp + 1; km = kp + 1; kp = kp + 1;
                    vmin = inp(k) - 2*lam;
                    vmax = inp(k);
                    umin = lam;
                    umax = - lam;
                else
          % #5
                    k = k + 1;
                    umin = umin + inp(k) - vmin;
                    umax = umax + inp(k) - vmax;
          % #6
                    if umin >= lam
                        vmin = vmin + (umin - lam)/(k - k0 + 1);
                        umin = lam;
                        km = k;
                    end
                    
                    if umax <= -lam
                        vmax = vmax + (umax + lam)/(k - k0 + 1);
                        umax = -lam;
                        kp = k;
                    end
                end
            end
          % #7
            if k >= N
          % #8
                if umin < 0
                    for i = k0:km
                        out(i) = vmin;
                    end
                    k = km + 1; k0 = km +1; km = km + 1;
                    vmin = inp(k);
                    umin = lam;
                    umax = inp(k) + lam - vmax;
                else
          % #9
                    if umax > 0
                        for i = k0:kp
                            out(i) = vmax;
                        end
                        k = kp + 1;   k0 = kp + 1; kp = kp + 1;
                        vmax = inp(k);
                        umax = - lam;
                        umin = inp(k) - lam - vmin;
                    else
                        for i = k0:N
                            out(i) = vmin + umin/(k - k0 + 1);
                        end
                        break;
                    end
                end
            end
        end






end