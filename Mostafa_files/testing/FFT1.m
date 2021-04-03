%  close all
%  clear all

% Ts = h;     
% t=TI_NL.Solution.time(:);
% x=TI_NL.Solution.u(dof,:)';


t=t1;
x=x1;
Ts=5.93568E-05;
figure(1000)
hold on
plot(t,x)
xlabel('Time (seconds)')
ylabel('Amplitude')
y = fft(x);   
fs = 1/Ts;
f = (0:length(y)-1)*fs/length(y);
figure(2000)
hold on
plot(f,abs(y))
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Magnitude')
figure(3000)
hold on
loglog(f,abs(y))
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Magnitude')

[y1,f1]=fft_n([t,x],fs);
figure(4000)
hold on
loglog(f1,abs(y1(:,2)))
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Magnitude')
figure(5000)
hold on
plot(f1,abs(y1(:,2)))
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Magnitude')