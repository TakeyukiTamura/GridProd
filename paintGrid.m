function [] = paintGrid(resultFile)
% 'paintGrid' outputs a heatmap that shows how much target metabolite
% production is achieved for each grid.
% 
%INPUTS 
% resultFile   The file output by 'analyzeResult'. The name of the file
%              should be 'analyzeResult_%d.mat' where %d is ID in the
%              model.
%OUTPUTS
% The heatmap figure is shown.
%
% Jul. 20, 2017   Takeyuki TAMURA
%
load(resultFile);
image(table2,'CDataMapping','scaled')
colormap(jet);
colorbar
grid on
grid minor
ax=gca;
ax.YDir='normal';
xlabel(ax,'Trange');
ylabel(ax,'Brange');
%%set(gca,'XLim',[0 380]);
end

