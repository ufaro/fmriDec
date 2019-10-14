clear all
close all
clc

% This is an example of use for the FastMapping_Tool
% Please verify that SPM (http://www.fil.ion.ucl.ac.uk/spm/) 
% is installed and added to the Matlab path. 

% Inputs and outputs are .MAT files. Please be sure that the input TCN.mat is a
% 2-dimesional matrix (T*N) where the first dimension is Time. 

% The output TC.mat is saved to the current folder

% If the number of voxels id low, think about replacing the "parfor" loop by
% a simple "for" loo.

%% Call Inputs and Setting file

InputSetting

%% Initialize parameters 

T = param.Tpoints;
NbrVoxels = param.NbrVoxels;

TR = param.TR;
niter = param.niter;              
condition = param.condition;     
lam0 = param.lam0;                




%% Make HRF operator

p = [6 16 1 1 6 0 30]; 
h=spm_hrf(TR,p);

hrf = [h' zeros(1,T-1)];
c = [hrf(1) zeros(1,T-1)];
xConv = toeplitz(c,hrf);
HRF = xConv';

H =@(u) HRF*u;
Ht =@(u) HRF'*u;

HH=HRF'*HRF;
maxeig=norm(HH,2);

Nh=length(h);

%% Main FISTA loop

TC = zeros(T,NbrVoxels);
nor = 1/.6745;

tic;

parfor v=1:NbrVoxels

    
y = TCN(:,v);    
    
% noise estimation

[coef,len] = wavedec(y,1,'db6');
coef(1:len(1)) = [];                                              
lam = mad(coef,1)*nor;


u = GPFM(y,H,Ht,maxeig,lam0*lam,niter,Nh,condition);

TC(:,v)=u;


end
ftime = toc;

fprintf('It took %f seconds to deconvolve  %d voxels ... \n',ftime,NbrVoxels);

output_dir=fullfile('OUTPUT/',func_name);

if (~exist(output_dir,'dir'))
fprintf('Results path does not exist %s \n creating new directory \n', output_dir);
mkdir(output_dir)
end


save([output_dir,'/TC_lambda_',...
            strrep(num2str(param.lam0,'%d'),'.',''),'_',date],'TC')

save([output_dir,'/TCN'],'TCN')
        
% 
printme = @(filename) print(fullfile(output_dir,filename),'-dpdf');

figure(1)
clf
imagesc(TC')
colormap(jet)
caxis([-2,2])
xlabel('timepoints [s/TR]');
ylabel('voxels')


printme(['TC_lambda_',strrep(num2str(param.lam0,'%d'),'.',''),'_',date])


figure(2)
clf
imagesc(TCN')
colormap(jet)
caxis([-2,2])
xlabel('timepoints [s/TR]');
ylabel('voxels')

printme(['TCN'])


