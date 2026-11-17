function plot_groups(group_number)
% Native one-factor curves; BASE is the black solid curve.
scope=setup_paths(); %#ok<NASGU>
root=fileparts(fileparts(mfilename('fullpath')));G=parameter_groups(root);g=G(group_number);
colors=[40 104 168;90 167 208;0 0 0;62 143 118;140 87 164]/255;styles={'--',':','-','-.','--'};markers={'o','s','^','d','v'};
symbols={'C_c','C_v','d_{\rm Nuc}','n_0','p_{\rm th}'};symbol=symbols{group_number};
fields={'u','ptotal_MPa','alpha'};labels={'u (m/s)','p_{total} (MPa)','\alpha'};q=cell(size(g.ids));names=cell(size(g.ids));
for j=1:numel(g.ids)
 q{j}=readtable(fullfile(root,'results/profiles',[g.ids{j} '_1ns.csv']));v=g.levels(j);
 if strcmp(g.name,'n0'),names{j}=sprintf('%s = 10^{%g} m^{-3}',symbol,round(log10(v)));elseif strcmp(g.name,'d'),names{j}=sprintf('%s = %g nm',symbol,v);elseif strcmp(g.name,'pth'),names{j}=sprintf('%s = %g MPa',symbol,v);else,names{j}=sprintf('%s = %g',symbol,v);end
end
f=figure('Visible','off','Color','w','Units','inches','Position',[.2 .2 15 4.7],'PaperPositionMode','auto');base=find(strcmp(g.ids,'BASE'));order=[setdiff(1:numel(g.ids),base) base];
for col=1:3
 ax=axes(f,'Position',[.065+.325*(col-1) .19 .265 .73]);hold(ax,'on');h=gobjects(size(g.ids));
 for j=order
  ci=g.color_indices(j);width=1.2;if ci==3,width=2.5;end
  h(j)=plot(ax,q{j}.r_um,q{j}.(fields{col}),'Color',colors(ci,:),'LineStyle',styles{ci},'LineWidth',width,'DisplayName',names{j},'Tag','native','UserData',struct('source',['results/profiles/' g.ids{j} '_1ns.csv'],'field',fields{col}));
 end
 for j=1:numel(g.ids)
  ci=g.color_indices(j);targets=40.5+(ci-1)*.65:4:59.8;idx=zeros(size(targets));for z=1:numel(targets),[~,idx(z)]=min(abs(q{j}.r_um-targets(z)));end
  plot(ax,q{j}.r_um(idx),q{j}.(fields{col})(idx),'LineStyle','none','Marker',markers{ci},'MarkerSize',4.8,'MarkerFaceColor','w','Color',colors(ci,:),'LineWidth',1,'Tag','native_markers','HandleVisibility','off','UserData',struct('source',['results/profiles/' g.ids{j} '_1ns.csv'],'field',fields{col},'indices',idx));
 end
 xlim(ax,[40 60]);xticks(ax,40:5:60);grid(ax,'on');box(ax,'on');xlabel(ax,'r (\mum)');ylabel(ax,labels{col});set(ax,'FontSize',12,'FontName','Arial');
 if col==3,legend(ax,h,names,'Location','best','FontSize',10,'Interpreter','tex');end
end
set(f,'PaperUnits','inches','PaperSize',[15 4.7],'PaperPosition',[0 0 15 4.7]);drawnow;name=['group_' g.name];out=fullfile(root,'figure/groups');if ~isfolder(out),mkdir(out);end;set(f,'Visible','on');savefig(f,fullfile(out,[name '.fig']));exportgraphics(f,fullfile(out,[name '.png']),'Resolution',250,'Width',15,'Height',4.7,'Units','inches','Padding','figure');close(f);
end
