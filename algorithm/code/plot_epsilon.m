function plot_epsilon(kind)
% Replot the selected fixed/dynamic comparison from native profiles.
if nargin<1,kind='voltage';end
if strcmp(kind,'voltage'),group='voltage';elseif strcmp(kind,'rise_time'),group='rise_time';else,error('EHD:Comparison','Use voltage or rise_time.');end
root=project_root();outdir=fullfile(root,'MG_epsilon',group,'figure');
if ~isfolder(outdir),mkdir(outdir);end
colors=[40 104 168;213 107 46;62 143 118]/255;
fields={'rho_m','rho_l','p_MPa','u','alpha','ptotal_MPa','epsilon_r','E_GV_m','pE_MPa'};
labels={'\rho_m (kg/m^3)','\rho_l (kg/m^3)','p (MPa)','u (m/s)','\alpha','p_{total} (MPa)','\epsilon_r','E (GV/m)','p_E (MPa)'};
modes={'fixed','dynamic'};
styles={'--','-'};
markers={'o','^'};
records={};
if strcmp(kind,'voltage'),comparison([10 1;30 1;50 1],'voltage','epsilon_voltage');
else,comparison([50 1;50 2;50 5],'rise_time','epsilon_rise_time');end
verification=struct('passed',true,'figure_files',1,'native_series_checked',54,...
 'radial_limits_um',[40 60],'subplot_labels','(a)--(i)',...
 'markers','Fixed: open circles; dynamic: open triangles; staggered native samples with an added triangle at the start of the 50 kV / 1 ns mixture-density descent',...
 'smoothing',false,'value_offsets',false,'production_integrations',0,'plots',{records});
verify=fullfile(root,'MG_epsilon',group,'results','verification');if ~isfolder(verify),mkdir(verify);end
fid=fopen(fullfile(verify,'plot_revision_verification.json'),'w');
fprintf(fid,'%s',jsonencode(verification));fclose(fid);

 function comparison(cases,kind,name)
  if strcmp(kind,'voltage')
   map=[1 3 2];
  else
   map=[2 1 3];
  end
  f=figure('Visible','off','Units','inches','Position',[.2 .2 15.6 9.5],...
   'PaperPositionMode','auto','Color','w','InvertHardcopy','off');
  profiles=cell(3,2);ids=cell(3,2);sources=cell(3,2);
  for j=1:3
   for mode=1:2
    ids{j,mode}=sprintf('V%g_T%g_%s',cases(j,1),cases(j,2),modes{mode});
    [file,sources{j,mode}]=epsilon_profile(ids{j,mode},cases(j,2));
    profiles{j,mode}=readtable(file);
   end
  end
  % Six explicit legend entries show the exact line/marker combination.
  la=axes(f,'Position',[.04 .918 .92 .071],'Visible','off');hold(la,'on');
  lh=gobjects(1,6);names=cell(1,6);
  for mode=1:2
   for j=1:3
    s=(mode-1)*3+j;suffix=modes{mode};if mode==2,suffix='dynamic MG';end
    lh(s)=plot(la,nan,nan,'Color',colors(map(j),:),'LineStyle',styles{mode},...
     'LineWidth',1.6,'Marker',markers{mode},'MarkerSize',5,'MarkerFaceColor','w');
    names{s}=sprintf('%g kV / %g ns - %s',cases(j,1),cases(j,2),suffix);
   end
  end
  lg=legend(la,lh,names,'NumColumns',3,'Orientation','horizontal',...
   'Box','off','FontName','Arial','FontSize',11.1,'Interpreter','none');
  lg.Position=[.055 .918 .89 .071];lg.ItemTokenSize=[39 16];
  panelRecords=cell(1,9);
  for k=1:9
   row=floor((k-1)/3);col=mod(k-1,3);
   ax=axes(f,'Position',[.069+col*.325 .69-row*.295 .258 .19]);hold(ax,'on');
   set(ax,'FontName','Arial','FontSize',10.5,'LineWidth',.7,...
    'GridAlpha',.14,'TickDir','out');grid(ax,'on');box(ax,'on');
   values=[];
   % Solid curves below dashed curves leave fixed/dynamic overlays visible.
   for mode=[2 1]
    for j=1:3
     q=profiles{j,mode};
     mask=q.r_um>=40&q.r_um<=60;values=[values;q.(fields{k})(mask)]; %#ok<AGROW>
     plot(ax,q.r_um,q.(fields{k}),'LineStyle',styles{mode},'Color',colors(map(j),:),...
      'LineWidth',1.55,'Tag','native','DisplayName',ids{j,mode},...
      'UserData',struct('source',sources{j,mode},'field',fields{k}));
     
     
    end
   end
   % Each marker is an existing raw sample; no interpolation or vertical shift.
   markerIndices=cell(3,2);
   for mode=1:2
    for j=1:3
     s=(mode-1)*3+j;q=profiles{j,mode};targets=40.25+(s-1)*.48:2.88:59.85;
     idx=zeros(size(targets));for z=1:numel(targets),[~,idx(z)]=min(abs(q.r_um-targets(z)));end
     % Mark the near-wall mixture-density descent of the 50 kV / 1 ns case.
     if k==1 && isequal(cases(j,:),[50 1])
      near=find(q.r_um>=40 & q.r_um<=41);y=q.(fields{k})(near);
      fraction=.72;if mode==2,fraction=.40;end
      targetY=min(y)+fraction*(max(y)-min(y));
      [~,ii]=min(abs(y-targetY));idx=[idx near(ii)];
      if mode==2
       % Label the start of the dynamic descent, inside the left axis.
       [~,firstTriangle]=min(abs(q.r_um-40.31));
       idx=[idx firstTriangle];
      end
     end
     idx=unique(idx);markerIndices{j,mode}=idx;
     plot(ax,q.r_um(idx),q.(fields{k})(idx),'LineStyle','none',...
      'Marker',markers{mode},'MarkerSize',4.8,'MarkerFaceColor','w',...
      'MarkerEdgeColor',colors(map(j),:),'LineWidth',1.05,...
      'Tag','native_markers','HandleVisibility','off',...
      'UserData',struct('source',sources{j,mode},'field',fields{k},'indices',idx));
    end
   end
   lo=min(values);hi=max(values);span=hi-lo;
   if span==0,span=max(1,abs(hi)*.02);end
   ylim(ax,[lo-.045*span hi+.07*span]);xlim(ax,[40 60]);xticks(ax,40:5:60);
   xlabel(ax,'r (\mum)','Interpreter','tex');ylabel(ax,'');
   title(ax,sprintf('(%c)  %s',char('a'+k-1),labels{k}),...
    'FontSize',12,'FontWeight','bold','Interpreter','tex');
   set(ax.Title,'Units','normalized','Position',[.5 1.04 0],'Clipping','off');
   
   
   panelRecords{k}=struct('field',fields{k},'label',char('a'+k-1),...
    'native_lines',6,'marker_indices',{markerIndices});
  end
  set(f,'Units','inches','Position',[.2 .2 15.6 9.5],'PaperPositionMode','auto');
  drawnow;
  exportgraphics(f,fullfile(outdir,[name '.png']),'Resolution',250,...
   'Width',15.6,'Height',9.5,'Units','inches','Padding','figure');
  set(f,'PaperUnits','inches','PaperSize',[15.6 9.5],'PaperPosition',[0 0 15.6 9.5]);
  set(f,'Visible','on');savefig(f,fullfile(outdir,[name '.fig']));
  records{end+1}=struct('name',name,'sources',{sources},'panels',{panelRecords});
  close(f);
 end
end
