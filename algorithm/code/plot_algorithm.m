function plot_algorithm(number)
package=project_root();catalog=table2struct(readtable(catalog_file()));out=fullfile(package,'algorithm','figure');if ~isfolder(out),mkdir(out);end
panels={
 'algorithm_BASE_B00','algorithm_BASE_B01',{'alpha','u'},[40 40.5;40.92 41.45],{'B00 aligned','B01'};
 'algorithm_S2_B01','algorithm_S2_B02',{'alpha','u'},[40 44;40 44],{'B01 Euler','B02 SSPRK3'};
 'algorithm_S1_B02','algorithm_S1_B03',{'u'},[40 44],{'B02 on matched B03 clock','B03 IMEX'};
 'algorithm_S07_B03','algorithm_S07_B04',{'u','u'},[40 44;40.28 41],{'B03','B04'}};
entry=panels(number,:);A=load(fullfile(package,catalog(find(strcmp({catalog.key},entry{1}),1)).folder,'state_1_ns.mat'));B=load(fullfile(package,catalog(find(strcmp({catalog.key},entry{2}),1)).folder,'state_1_ns.mat'));
r=A.grid.centers*1e6;C=A.C;
heights=round(1300*[985/2297 956/2165 1100/1634 1123/1828]);
f=figure('Visible','off','Color','w','Position',[20 20 1300 heights(number)]);fields=entry{3};windows=entry{4};
annotation(f,'textbox',[.04 .92 .92 .06],'String',sprintf('%s to %s | t = 1 ns | %g kV | N%d',strrep(entry{1},'algorithm_',''),strrep(entry{2},'algorithm_',''),C.Umax/1000,C.Nr),'Interpreter','none','EdgeColor','none','HorizontalAlignment','center','FontSize',14);
annotation(f,'textbox',[.04 .87 .92 .05],'String',sprintf('MG | C_c=%g | C_v=%g | d_{Nuc}=%g nm | n_0=%.3g m^{-3} | u(a)=0, du/dr(b)=0',C.vofCc,C.vofCv,C.vofDNuc*1e9,C.vofN0),'EdgeColor','none','HorizontalAlignment','center','FontSize',12);
for k=1:numel(fields)
    if isscalar(fields),pos=[.08 .17 .86 .64];else,pos=[.075+.49*(k-1) .17 .40 .64];end
    ax=axes(f,'Position',pos);mask=r>=windows(k,1)&r<=windows(k,2);key=fields{k};a=A.state.(key)(mask);b=B.state.(key)(mask);
    plot(ax,r(mask),a,'--o','Color',[.35 .35 .35],'LineWidth',1.5,'MarkerSize',4);hold(ax,'on');plot(ax,r(mask),b,'-s','Color',[0 .22 .85],'LineWidth',1.6,'MarkerSize',4);
    
    if strcmp(key,'alpha'),label='Liquid fraction \alpha';else,label='Radial velocity u (m/s)';end
    ylabel(ax,label);xlabel(ax,'Radius r (\mum)');xlim(ax,windows(k,:));grid(ax,'on');box(ax,'on');set(ax,'FontName','Arial','FontSize',12);legend(ax,entry{5},'Location','best');
end
set(f,'PaperPositionMode','auto','PaperUnits','inches');p=get(f,'PaperPosition');set(f,'PaperSize',p(3:4));drawnow;
name=sprintf('step_%02d',number);drawnow;set(f,'Visible','on');savefig(f,fullfile(out,[name '.fig']));close(f);f=openfig(fullfile(out,[name '.fig']),'invisible');canvas=f.PaperSize;set(f,'Units','inches','Position',[.2 .2 canvas(1) canvas(2)],'PaperPositionMode','manual','PaperPosition',[0 0 canvas]);drawnow;print(f,fullfile(out,[name '.png']),'-dpng','-r250','-painters');close(f);
end
