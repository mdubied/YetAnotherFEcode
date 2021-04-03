%% This function does the normalisation of the output of the fft function
% Input:
%   - data: input data matrix with size [r,c], with r samples and c signals
%   - fsamp: sampling frequency
% Output:
%   - norm_sp: normalized spectrum (positive frequencies)
%   - freq_vec: frequency vector

function [norm_sp, freq_vec]=fft_n(data,fsamp)

dim=size(data);

if dim(2)>dim(1)
    data=data';
end

N=length(data);
df=fsamp/N;

if (N/2)==(floor(N/2))
    
    freq_vec=[0:df:(N/2*df)]';
    NF=length(freq_vec);
    sp=fft(data,[],1);
    norm_sp(1,:)=sp(1,:)/N;
    norm_sp(2:N/2,:)=sp(2:N/2,:)/(N/2);
    norm_sp(N/2+1,:)=sp(N/2+1,:)/N;
    
else
    
    freq_vec=[0:df:((N-1)/2)*df]';
    NF=length(freq_vec);
    sp=fft(data,[],1);
    norm_sp(1,:)=sp(1,:)/N;
    norm_sp(2:(N+1)/2,:)=sp(2:(N+1)/2,:)/(N/2);
    
end




