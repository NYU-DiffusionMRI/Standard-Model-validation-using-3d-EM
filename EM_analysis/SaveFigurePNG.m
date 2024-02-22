function SaveFigurePNG(figname,w,h)
% figname: path/name.png for the Figure
%     w,h: width and height of the Figure (centimeters) 
if(nargin <2)
w = 20;
h = 10;
end
%%%%%%%%%%%%%%%%%%%%
fig = gcf;
fig.PaperPositionMode = 'auto';
%fig.PaperUnits = 'inches';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 w h];
print(figname,'-dpng','-r0')
%%%%%%%%%%%%%%%%%%%%%%
end