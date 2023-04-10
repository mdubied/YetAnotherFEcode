function AnimateFieldonDeformedMeshActuation2Muscles(Nodes,Elements,ActuationElements,ActuationValues,ActuationElements2,ActuationValues2,S,varargin)
%% function to animate the displacement snapshots in S (solution cell or matrix) 
% More than one solution signals can be simultaneously animated,
% if S is a cell, then each cell component corresponds to a different
% solution matrices and if S is a matrix then it is a matrix whose columns contain
% the displacement snapshots of a single solution

% Additional Name Value Parameters:
% index: numerical array specifying which DOFs on the node correspond to translational
% displacements in X, Y, Z directions for 3D meshes, X,Y direction for 2D
% meshes and simply a X direction for 1 D meshes
% factor: the factor with which displacements should be scaled
% cameraPos: camera placement using view() for 3D plot
% upVec: set vectors that should be oriented vertically
% filename: for storing animation files (with path)
% framerate: numerical rate of frames to be played per second in the video
%
% Last modified: 22/03/2023, Mathieu Dubied, ETH Zurich
figure
[scalefactor,index,cameraPos,upVec,filename,framerate] = parse_inputs(varargin{:});

%% video object
nnodes = size(Nodes,1);
myVideo = VideoWriter(filename,'MPEG-4');
myVideo.FrameRate = framerate;

if iscell(S)
    ns = length(S);
    s = zeros(ns,1);
    for j = 1:ns
        s(j) = size(S{j},2);
    end
    nt = min(s); % collect minimal number of snapshots across solutions
    
    for k = 1:ns     % check if all solutions have same number of snapshots and prune if necessary   
        if size(S{k},2) ~= nt
            warning(['Solution cell components do not have same size: truncating at first' num2str(nt) 'snapshots'] )
            S{k} = S{k}(:,1:nt);
        end
    end
else
    nt = size(S,2);
    S = {S};
    ns = 1;
end


nDOFperNode = size(S{1},1)/nnodes;

M(size(S{1},2)) = struct('cdata',[],'colormap',[]);
color = get(groot, 'defaultAxesColorOrder');
color = [0 0 0; color];
a1 = [];
a2 = [];
a3 = [];
a4 = [];
a5 = [];
normalizedActuationValues = normalize(ActuationValues,'range',[-1 1]);
normalizedActuationValues2 = normalize(ActuationValues2,'range',[-1 1]);

for j = 1:nt
    hold on
    for k = 1:ns
        Solution = S{k};
        meshcolor = color(k,:);
        U = reshape(Solution(:,j),nDOFperNode,[]).';
        disp = U(:,index);
        
        PlotFieldonDeformedMeshActuation2Muscles(Nodes,Elements,ActuationElements,normalizedActuationValues(j),ActuationElements2,normalizedActuationValues2(j),disp,'factor',scalefactor,'cameraPos',cameraPos,'upVec',upVec,'color', meshcolor) ;
        delete(a1);
        delete(a2);
        delete(a3);
        delete(a4);
        delete(a5);
        % muscle 1
        a1 = annotation('textbox', [0.2, 0.2, 0.25, 0.06], 'String', "muscle 1, a=" + ActuationValues(j));
        if normalizedActuationValues(j) >=0
            a2 = annotation('rectangle',[0.155, 0.2, 0.045, 0.06],'FaceColor','red','FaceAlpha',abs(normalizedActuationValues(j)));
        else
            a2 = annotation('rectangle',[0.155, 0.2, 0.045, 0.06],'FaceColor','blue','FaceAlpha',abs(normalizedActuationValues(j)));
        end
        % muscle 2
        a3 = annotation('textbox', [0.2, 0.1, 0.25, 0.06], 'String', "muscle 2, a=" + ActuationValues2(j));
        if normalizedActuationValues2(j) >=0
            a4 = annotation('rectangle',[0.155, 0.1, 0.045, 0.06],'FaceColor','red','FaceAlpha',abs(normalizedActuationValues2(j)));
        else
            a4 = annotation('rectangle',[0.155, 0.1, 0.045, 0.06],'FaceColor','blue','FaceAlpha',abs(normalizedActuationValues2(j)));
        end
        a5 = annotation('textbox', [0.5, 0.2, 0.25, 0.06], 'String', "time: " + num2str(j/framerate,'%4.2f') + "s");
        if j == 1
            brighten(0.6)
        end
        if j*k == 1
            %             set(gca, 'nextplot', 'replacechildren')
            axis manual
            xlim = get(gca,'xlim');
            ylim = get(gca,'ylim');
        end
        set(gca,'xlim',[-xlim(2)*0.2, xlim(2)*1.2],'ylim',ylim*1.5)
    end
    % gif movie
    frame = getframe(gcf);
    M(j) = frame;
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Write to the GIF File
    if j == 1
        %         set(gca, 'nextplot', 'replacechildren');
        axis manual
        imwrite(imind,cm,[filename '.gif'],'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append');
    end
    cla
    hold off
end
open(myVideo)
writeVideo(myVideo,M);
close(myVideo)
close(gcf)
end

function [scalefactor,index,cameraPos,upVec,filename,framerate] = parse_inputs(varargin)
%% parsing inputs
defaultindex = 1;
defaultFactor = 1;
defaultCameraPos = 3;
defaultUpVec = [0;1;0];
defaultfilename = 'test'; % plot norm of displacement\
defaultframerate = 100;

p = inputParser;
addParameter(p,'index',defaultindex, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty','integer','positive'}) );
addParameter(p,'factor',defaultFactor,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty','positive'}) );
addParameter(p,'cameraPos',defaultCameraPos,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );
addParameter(p,'upVec',defaultUpVec,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );
addParameter(p,'filename',defaultfilename,@(x)validateattributes(x, ...
                {'char'},{'nonempty'}))
addParameter(p,'framerate',defaultframerate,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty','integer','positive'}) );

parse(p,varargin{:});

scalefactor = p.Results.factor;
index = p.Results.index;
cameraPos= p.Results.cameraPos;
upVec = p.Results.upVec;
filename = p.Results.filename;
framerate = p.Results.framerate;
end