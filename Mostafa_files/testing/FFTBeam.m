%  close all
%  clear all

%   This code perform fourier transform on the results of 2D Beam
%   the 4 points are to be monitored A,B and C midpoint of each beam 
%   and the point where the force is being actuated 
%   Use TI_NL if the variables are still in the workspace or Load the Struct
%   file from the directory saved.
%   The assembly should be in the workspace to run the first section or 
%   save offline the dof of each point in the matrix dofM.
%%
dofM=[];
nfsense = find_node(Lx/2,Ly/2,[],nodes); % node where to see the results
DOF = get_index(nfsense, myMesh.nDOFPerNode );
dofM(1)=DOF(2);


nf = find_node(Lx/2,Ly,[],nodes); % node where to put the force
node_force_dofs = get_index(nf, myMesh.nDOFPerNode );
dofM(2)=node_force_dofs(2);

%%
Ts = h;    

for i=1:2
    dof=dofM(i)
t=TI_NL.Solution.time(:);
x=TI_NL.Solution.u(dof,:)';

% t=time(:);
% x=u(dof,:)';


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

% figure(3000)
% hold on
% loglog(f,abs(y))
% xlabel('Frequency (Hz)')
% ylabel('Magnitude')
% title('Magnitude')

[y1,f1]=fft_n([t,x],fs);
figure(4000)
hold on
loglog(f1,abs(y1(:,2)))
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Magnitude')

% figure(5000)
% hold on
% plot(f1,abs(y1(:,2)))
% xlabel('Frequency (Hz)')
% ylabel('Magnitude')
% title('Magnitude')
end