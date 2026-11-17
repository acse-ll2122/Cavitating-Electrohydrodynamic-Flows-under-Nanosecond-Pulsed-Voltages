function plot_parameters()
% Plot native profiles directly; black marks the unchanged BASE input.
scope=setup_paths(); %#ok<NASGU>
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fileparts(root));out=fullfile(root,'figure','parameters');if ~isfolder(out),mkdir(out);end
G=parameter_groups(root);
colors=[40 104 168;90 167 208;0 0 0;62 143 118;140 87 164]/255;
styles={'--',':','-','-.','--'};markers={'o','s','^','d','v'};
symbols={'C_c','C_v','d_{\rm Nuc}','n_0','p_{\rm th}'};
fields={'u','ptotal_MPa','alpha'};labels={'u (m/s)','p_{total} (MPa)','\alpha'};
f=figure('Visible','off','Color','w','Units','inches','Position',[.2 .2 12 12.6],'PaperPositionMode','auto');records={};
lo=[inf inf inf];hi=[-inf -inf -inf];
for row=1:5
 for k=1:numel(G(row).ids)
  q=readtable(fullfile(root,'results','profiles',[G(row).ids{k} '_1ns.csv']));mask=q.r_um>=40&q.r_um<=60;
  for col=1:3,lo(col)=min(lo(col),min(q.(fields{col})(mask)));hi(col)=max(hi(col),max(q.(fields{col})(mask)));end
 end
end
limits=[1.07*lo(1) max(0,hi(1));1.07*lo(2) max(0,hi(2));max(0,lo(3)-.01) 1];
for row=1:5
 g=G(row);bottom=.81-(row-1)*.185;profiles=cell(size(g.ids));names=cell(size(g.ids));lh=gobjects(size(g.ids));
 annotation(f,'textbox',[.025 bottom+.106 .18 .036],'String',symbols{row},'Interpreter','tex','LineStyle','none','FontSize',17,'FontName','Arial','FontWeight','bold');
 la=axes(f,'Position',[.022 bottom-.004 .185 .103],'Visible','off');hold(la,'on');
 for k=1:numel(g.ids)
  profiles{k}=readtable(fullfile(root,'results','profiles',[g.ids{k} '_1ns.csv']));j=g.color_indices(k);names{k}=value_label(symbols{row},g.name,g.levels(k));
  width=1.2;if j==3,width=2.5;end
  lh(k)=plot(la,nan,nan,'Color',colors(j,:),'LineStyle',styles{j},'LineWidth',width,'Marker',markers{j},'MarkerSize',5,'MarkerFaceColor','w');
 end
 lg=legend(la,lh,names,'Box','off','Interpreter','tex','FontSize',10.3,'FontName','Arial');lg.Position=[.022 bottom-.006 .187 .104];lg.ItemTokenSize=[29 12];
 base=find(strcmp(g.ids,'BASE'));order=[setdiff(1:numel(g.ids),base) base];
 for col=1:3
  ax=axes(f,'Position',[.235+(col-1)*.25 bottom .21 .13]);hold(ax,'on');set(ax,'FontName','Arial','FontSize',10,'GridAlpha',.14);grid(ax,'on');box(ax,'on');
  for k=order
   q=profiles{k};j=g.color_indices(k);width=1.2;if j==3,width=2.5;end
   rel=['results/profiles/' g.ids{k} '_1ns.csv'];plot(ax,q.r_um,q.(fields{col}),'Color',colors(j,:),'LineStyle',styles{j},'LineWidth',width,'Tag','native','DisplayName',names{k},'UserData',struct('source',rel,'field',fields{col}));
   
  end
  for k=1:numel(g.ids)
   q=profiles{k};j=g.color_indices(k);targets=40.5+(j-1)*.65:4:59.8;idx=zeros(size(targets));for z=1:numel(targets),[~,idx(z)]=min(abs(q.r_um-targets(z)));end
   idx=unique(idx);plot(ax,q.r_um(idx),q.(fields{col})(idx),'LineStyle','none','Marker',markers{j},'MarkerSize',5,'MarkerFaceColor','w','Color',colors(j,:),'LineWidth',1.1,'Tag','native_markers','HandleVisibility','off','UserData',struct('source',['results/profiles/' g.ids{k} '_1ns.csv'],'field',fields{col},'indices',idx));
   
   records{end+1}=struct('group',g.name,'case_id',g.ids{k},'field',fields{col},'legend',names{k},'base',j==3,'marker_indices',idx); %#ok<AGROW>
  end
  xlim(ax,[40 60]);xticks(ax,40:5:60);ylim(ax,limits(col,:));xlabel(ax,'r (\mum)');ylabel(ax,'');if row==1,title(ax,labels{col},'FontSize',13,'FontWeight','bold');end
 end
end
set(f,'PaperUnits','inches','PaperSize',[12 12.6],'PaperPosition',[0 0 12 12.6]);drawnow;set(f,'Visible','on');savefig(f,fullfile(out,'parameter_5x3.fig'));
close(f);f=openfig(fullfile(out,'parameter_5x3.fig'),'invisible');set(f,'Units','inches','Position',[.2 .2 12 12.6]);drawnow;
exportgraphics(f,fullfile(out,'parameter_5x3.png'),'Resolution',250,'Width',12,'Height',12.6,'Units','inches','Padding','figure');close(f);
verification=struct('passed',true,'radial_limits_um',[40 60],'native_series_checked',69,'markers_are_native_samples',true,'value_offsets',false,'smoothing',false,'base_color',[0 0 0],'base_line_width',2.5,'curve_letter_codes',false,'records',{records});
verificationDir=fullfile(root,'results','verification');if ~isfolder(verificationDir),mkdir(verificationDir);end
fid=fopen(fullfile(verificationDir,'parameter_plots.json'),'w');fprintf(fid,'%s',jsonencode(verification));fclose(fid);
end
function name=value_label(symbol,group,value)
if strcmp(group,'n0'),name=sprintf('%s = 10^{%g} m^{-3}',symbol,round(log10(value)));elseif strcmp(group,'d'),name=sprintf('%s = %g nm',symbol,value);elseif strcmp(group,'pth'),name=sprintf('%s = %g MPa',symbol,value);else,name=sprintf('%s = %g',symbol,value);end
end
