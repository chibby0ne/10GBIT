clear all;
legendval = {'' };
graphtitle = 'WPAN';
figfilename = 'WPAN';
xtitle = 'E_b/N_0 / dB';
ytitle = 'FER';
%Farbvektor initialisieren, entsprechende Kurven haben gleiche Farben;
		farb=zeros(1,3);
		%Moegliche Marker
		mark=['+','o','*','s','d','h'];
		%Graph erzeugen
		figf=figure('Visible','off');
		for k=1:length(legendval)
		    %Graph zeichnen, Farbe zufaellig, marker aus Vektor mark
		    plot(xval{k}(:), yval{k}(:), 'linestyle','--', 'linewidth',1.2, 'color',rand(1,3), 'marker',mark(mod(k,6)+1), 'DisplayName',legendval{k});
		    hold on;
		end
		%Graph formatieren
		xlabel(xtitle);
		ylabel(ytitle);
		title(graphtitle);
		axis tight;
		grid on;
set(gca,'YScale','log');
set(gca,'YLimMode','auto');
		set(gca,'Position',[0.05511 0.07292 0.9389 0.8847]);  %fits for create_pdf.m, paper size 24x18
		box on;
		hold off;
		leg = legend('toggle');
		legend('boxon');
		set(leg,'Location','southwest','FontSize',10);
		% Grafiken anzeigen
		set(figf,'Visible','on');
		% Grafik speichern
		%saveas(figf,[figfilename, '.fig'],'fig');
		%h = gcf();
		%set(h,'PaperSize',[24 18]);
		%set(h,'PaperUnits','centimeters');
		%set(h,'PaperOrientation','portrait');
		%set(h,'PaperPosition',[0 0 24 18]);
		%print('-dpdf',[figfilename, '.pdf']);
