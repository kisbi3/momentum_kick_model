%figure format

fig_width = 8; % inch
fig_height = 3; % inch
margin_vert_u = 0.3; % inch
margin_vert_d = 0.45; % inch
margin_horz_l = 0.25; % inch
margin_horz_r = 0.25; % inch
fontname = 'times new roman';
set(0,'defaultaxesfontname',fontname);
set(0,'defaulttextfontname',fontname);
fontsize = 9; % pt
set(0,'defaultaxesfontsize',fontsize);
set(0,'defaulttextfontsize',fontsize);
set(0,'fixedwidthfontname','times');
hFig = figure;
set(hFig,'renderer','painters');
set(hFig,'units','inches');
set(hFig,'position',[3 6 fig_width fig_height]);
set(hFig,'PaperUnits','inches');
set(hFig,'PaperSize', [fig_width fig_height]);
set(hFig,'PaperPositionMode', 'manual');
set(hFig,'PaperPosition',[0 0 fig_width fig_height]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pt = linspace(0,10,20);
  
Tjet_STAR = 0.55 + 0.*pt;
Tjet_PHENIX = 0.19 + 0.06.*pt;
Tjet_CMS = 0.228 + 0.072.*pt;
Tjet_pp = 0.266 + 0.084.*pt;
subplot(1,3,1);
hold on
ax = plot(pt, Tjet_STAR, 'r --', pt, Tjet_PHENIX, '. b :', pt, Tjet_pp, 'g -.', pt, Tjet_CMS, 'k -');
set(gca,'XTick',(0:2:12));
        title('(a)');
        xlabel('p_T^{trig}','FontSize',8,'Interpreter','tex');
        ylabel('$T_{jet}$','FontSize',12,'Interpreter','latex');
        legend(ax,'STAR','PHENIX','CMS (pp)','CMS (PbPb)','Location','best')

Njet_STAR = 0.75 + 0.*pt;
Njet_PHENIX = 0.15 + 0.10.*pt;
Njet_CMS = 0.602 + 0.17.*pt;
subplot(1,3,2);
hold on
plot(pt, Njet_STAR, 'r --', pt, Njet_PHENIX, '. b :', pt, Njet_CMS, 'k -');
set(gca,'XTick',(0:2:12));
        title('(b)');
        xlabel('p_T^{trig}','FontSize',8,'Interpreter','tex');
        ylabel('$N_{jet}$','FontSize',12,'Interpreter','latex');

fRNk_STAR = 3.8 + 0.*pt;
fRNk_PHENIX = 3.0 + 0.*pt;
fRNk_CMS = exp(-1.137./pt) .* 10.69652482 .* exp(-0.1692 .* pt);
fRNk_pp = 1.5 + 0.*pt;
subplot(1,3,3);
hold on
plot(pt, fRNk_STAR, 'r --', pt, fRNk_PHENIX, '. b :', pt, fRNk_pp, 'g -.', pt, fRNk_CMS, 'k -');
set(gca,'XTick',(0:2:12));
        title('(c)');
        xlabel('p_T^{trig}','FontSize',8,'Interpreter','tex');
        ylabel('$f_R\langle N_{k} \rangle$','FontSize',12,'Interpreter','latex');

saveas(hFig,'Test.pdf')