function [u] = GPFM(x,H,Ht,maxeig,lam,Nit,Nh,condition)


% Fast Mapping with Lasso for spike-like activity and fused Lasso 
% for block-like like activity. 

%   u = argmin_{u} 0.5*||Hu- x||_2^2 + lam * ||L{u}||_1

%Inputs
%   x is the input fMRI matrix (Time * voxels). 
%   H and Ht, the toeplitz HRF matrix and its transpose.
%   Condition : spikes (L = Id) or blocks (L = D). D = gradient

% For the spikes case this is 
% Gaudes, C?sar Caballero, et al. "Structured sparse deconvolution for paradigm 
% free mapping of functional MRI data." IEEE ISBI, 2012

% For the blocks case, this is fused Lasso 
% R. Tibshirani et al. J. R. Statist. Soc. B (2005) 67, Part 1, pp. 91?108.



% The solution is found via FISTA iterations. In the fused lasso case, the
% internal iteration use Condat's single-shot algorithm
%`Condat, L. (2013). A direct algorithm for 1-D total variation denoising. IEEE Signal Processing Letters, 20(11), 1054?1057. doi:10.1109/LSP.2013.2278339 <http://dx.doi.org/10.1109/LSP.2013.2278339>`_


% Needs SPM !


% Younes Farouj, @ MIPLAB-EPFL,
% May 30, 2018




y=[x; zeros(Nh-1,1)];


k = 1;
t = 1;
s = zeros(length(x),1);
u = zeros(length(x),1);
    

    if condition == 'blocks'

        

            while (k <= Nit)

              u_l = u;
              z = Ht(y)/maxeig + s - Ht(H(s))/maxeig;
              u = fast_tvprox(z,lam/maxeig);
              t_l = t;
              t = (1+sqrt(1+4*(t^2)))/2;
              s = u + (t_l - 1)/t*(u - u_l);  

              k = k+1;
            end
    end


    if condition == 'spikes'
       tau=lam/maxeig;

            while (k <= Nit)

              u_l = u;
              z = Ht(y)/maxeig + s - Ht(H(s))/maxeig;
              u = max(0,1-tau./max(abs(z),1e-10)).*z; u(u<0)=0;
              t_l = t;
              t = (1+sqrt(1+4*(t^2)))/2;
              s = u + (t_l - 1)/t*(u - u_l);  

              k = k+1;
            end
    end

    
    
    
end
    


