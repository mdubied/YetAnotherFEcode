function SetColorbar
cbar = colorbar;
% Dimensions of the colorbar     
cpos = get(cbar,'position'); 
cpos(1) = cpos(1)+.15;
cpos(3) = cpos(3) ;   % Reduce the width of colorbar by half
cpos(2) = cpos(2)+.25 ;
cpos(4) = cpos(4)-.4 ;
set(cbar,'Position',cpos) ;
%brighten(0.1); 
     
% Title of the colorbar
% set(get(cbar,'title'),'string','disp [cm]');
%locate = get(cbar,'title');
%tpos = get(locate,'position');
%tpos(3) = tpos(3)+5. ;
%set(locate,'pos',tpos);

% Setting the values on colorbar
%
% get the color limits
clim = caxis;
ylim(cbar,[clim(1) clim(2)]);
numpts = 0;%24 ;    % Number of points to be displayed on colorbar
cbar.Title.Interpreter = 'latex';
kssv = linspace(clim(1),clim(2),numpts);
set(cbar,'YtickMode','manual','YTick',kssv); % Set the tickmode to manual
for i = 1:numpts
    %imep = num2str(kssv(i),'%+3.2E');
    imep = num2str(kssv(i)*100);
    vasu(i) = {imep} ;
end
set(cbar);
% set(cbar,'YTickLabel',vasu(1:numpts),'fontsize',9,'TickLabelInterpreter','latex');
% a =  cbar.Position; %gets the positon and size of the color bar
% set(cbar,'Position',[a(1)+0.1 a(2) 0.02 0.6])% To change size