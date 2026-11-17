function plot_force()
% Plot only accepted native force-density components, retaining signed values.
code=fileparts(mfilename('fullpath'));package=fileparts(code);out=fullfile(package,'figure','force');if ~isfolder(out),mkdir(out);end
fixed=readtable(fullfile(package,'results/force/fixed_F_split.csv'));dynamic=readtable(fullfile(package,'results/force/dynamic_F_split.csv'));
f=figure('Visible','off','Color','w','Units','inches','Position',[.3 .3 12.2 4.2],'PaperPositionMode','auto');
for k=1:2
    ax=subplot(1,2,k);if k==1,q=fixed;heading='Fixed epsilon_r = 81';else,q=dynamic;heading='Dynamic Maxwell-Garnett';end
    show(ax,q);xlim(ax,[40 200]);title(ax,heading,'FontWeight','normal');
end
save_plot(f,out,'full_domain');
f=figure('Visible','off','Color','w','Units','inches','Position',[.3 .3 12.2 4.2],'PaperPositionMode','auto');
ax=subplot(1,2,1);show(ax,dynamic);xlim(ax,[40 60]);title(ax,'Dynamic MG two-term force density','FontWeight','normal');
ax=subplot(1,2,2);plot(ax,dynamic.r_um,dynamic.F1_TN_m3,'-','Color',[.08 .35 .70],'LineWidth',1.6);
grid(ax,'on');box(ax,'on');xlim(ax,[40 60]);xlabel(ax,'Radius r (\mum)');ylabel(ax,'F_1 (10^{12} N/m^3)');title(ax,'Dielectric-gradient term enlarged','FontWeight','normal');set(ax,'FontName','Arial','FontSize',11);
save_plot(f,out,'near_wall');
end
function show(ax,q)
h1=plot(ax,q.r_um,q.F1_TN_m3,'-','Color',[.08 .35 .70],'LineWidth',1.5);hold(ax,'on');h2=plot(ax,q.r_um,q.F2_TN_m3,'--','Color',[.88 .37 .12],'LineWidth',1.5);h3=plot(ax,q.r_um,q.F_TN_m3,':','Color',[.05 .05 .05],'LineWidth',1.5);

xlabel(ax,'Radius r (\mum)');ylabel(ax,'Force density (10^{12} N/m^3)');grid(ax,'on');box(ax,'on');set(ax,'FontName','Arial','FontSize',11);legend(ax,[h1 h2 h3],{'F_1 dielectric gradient','F_2 electrostrictive','F total'},'Location','best','FontSize',9);
end
function save_plot(f,out,name)
drawnow;set(f,'Units','inches','Position',[.3 .3 12.2 4.2],'PaperPositionMode','auto');drawnow;
print(f,fullfile(out,[name '.png']),'-dpng','-r250','-painters');set(f,'Visible','on');savefig(f,fullfile(out,[name '.fig']));close(f);
end
