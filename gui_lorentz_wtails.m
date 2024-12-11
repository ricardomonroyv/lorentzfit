function gui_lorentz_wtails
% SIMPLE_GUI2 Select a data set from the pop-up menu, then
% click one of the plot-type push buttons. Clicking the button
% plots the selected data in the axes.
% Create and then hide the UI as it is being constructed.
f = figure('Visible','off','Position',[360,500,450,285]);

% Construct the components.
%hsurf = uicontrol('Style','pushbutton',...
%        'String','Surf','Position',[315,220,70,25],...
%        'Callback',@surfbutton_Callback);
% hmesh = uicontrol('Style','pushbutton',...
%         'String','Mesh','Position',[315,180,70,25],...
%         'Callback',@meshbutton_Callback);
% hcontour = uicontrol('Style','pushbutton',...
%         'String','Contour','Position',[315,135,70,25],...
%         'Callback',@contourbutton_Callback);
hqt = uicontrol('Style','pushbutton',...
        'String','Calculate QL','Position',[350,100,70,25],...
        'Callback',@qfbutton_Callback);
htext = uicontrol('Style','text','String','Select a folder',...
        'Position',[325,60,60,15]);
    
Fldrfiles=dir('*spectra_*');
Fldrnms={};
for i=1:length(Fldrfiles)
    Fldrnms{i}=Fldrfiles(i).name;
end

hpopup = uicontrol('Style','popupmenu',...
        'String',Fldrnms,...
        'Position',[300,30,100,25],...
        'Callback',@popup_menu_Callback);
ha = axes('Units','pixels','Position',[50,60,200,185]);
align([hqt,htext,hpopup],'Center','None');

% Initialize the UI.
% Change units to normalized so components resize automatically.
f.Units = 'normalized';
ha.Units = 'normalized';
%hsurf.Units = 'normalized';
%hmesh.Units = 'normalized';
%hcontour.Units = 'normalized';
hqt.Units = 'normalized';
htext.Units = 'normalized';
hpopup.Units = 'normalized';

% Generate the data to plot.
%Create a Simple UI Programmatically
%3-11




% Assign the a name to appear in the window title.
f.Name = 'Simple GUI';

% Move the window to the center of the screen.
movegui(f,'center')

% Make the window visible.
f.Visible = 'on';
set(gcf,'units','normalized','outerposition',[0 0 1 1])
%f.WindowState = 'Maximized';
% Pop-up menu callback. Read the pop-up menu Value property to

function popup_menu_Callback(source,eventdata)
% Determine the selected data set.
str = get(source, 'String');
val = get(source,'Value');
global Adat
foldername=str{val};
tracesave=sprintf('%s/%s/',pwd,foldername);
filenameS=sprintf('%s*AmpvsFreq-spectra-*',tracesave);
Files=dir(filenameS);
pos=1; 
IDfilenameS=sprintf('%s%s',tracesave,Files(pos).name);
fid=fopen(IDfilenameS,'r');
datap(1,1)=textscan(fid,'%f %f','HeaderLines',1,'Delimiter',' ','CollectOutput',1);
fclose(fid);
Adat=datap{1,1};


freqfull=Adat(:,1);
ampfull=Adat(:,2);
%stdfreqfit=Ap(:,3);
plot(freqfull*10^-9,ampfull);
%set(gca,'YScale','log','FontSize',30);
set(gca,'FontSize',25);
grid on;
xlim([0.4 1.5])
ylim([-150 0])
xlabel('Frequency [GHz]')
ylabel('Amplitude [dB]')
% Set current data to the selected data set.
% switch str{val}
%     case 'Peaks' % User selects Peaks.
%     current_data = peaks_data;
%     case 'Membrane' % User selects Membrane.
%     current_data = membrane_data;
%     case 'Sinc' % User selects Sinc.
%     current_data = sinc_data;
% end
end


function qfbutton_Callback(source,eventdata)
% Display mesh plot of the currently selected data.
global Adat
figstruct=get(f);
figaxes=figstruct.CurrentAxes;
xlimfig=figaxes.XLim;
freqfull=Adat(:,1);
ampfull=Adat(:,2);

[~,finit] = min(abs(freqfull*10^-9-xlimfig(1)));
[~,ffins] = min(abs(freqfull*10^-9-xlimfig(2)));
new_freq=freqfull(finit:ffins);
new_amp=ampfull(finit:ffins);

%
[~,amp0] =max(new_amp);
[~,fmnus3dB] = min(abs(new_amp(1:amp0)-(new_amp(amp0)-3)));
[~,fplus3dB_b] = min(abs(new_amp(amp0:end)-(new_amp(amp0)-3)));
fplus3dB=amp0+fplus3dB_b-1;
bndwidth=abs(new_freq(fplus3dB)-new_freq(fmnus3dB));
quafL=new_freq(amp0)/abs(new_freq(fplus3dB)-new_freq(fmnus3dB));

[~,fmnus1dB] = min(abs(new_amp(1:amp0)-(new_amp(amp0)-1)));
[~,fplus1dB_b] = min(abs(new_amp(amp0:end)-(new_amp(amp0)-1)));
fplus1dB=amp0+fplus1dB_b-1;
bndwidth1=abs(new_freq(fplus1dB)-new_freq(fmnus1dB));
quafL1=new_freq(amp0)/abs(new_freq(fplus1dB)-new_freq(fmnus1dB))/2;
quafL1_L=new_freq(amp0)/(2*abs(new_freq(amp0)-new_freq(fmnus1dB)))/2;
quafL1_R=new_freq(amp0)/(2*abs(new_freq(fplus1dB)-new_freq(amp0)))/2;

freqdata=new_freq*10^-9; 
ampdata=10.^(new_amp/10);
ampdatanorm=ampdata/max(ampdata);
%options = optimset('PlotFcns',@optimplotfval);
arrayfreqdata=freqdata(1):3*10^-9:freqdata(end);
%
t = @(aa,bd,f0,x) abs(	aa./(1+1i*(1/bd)*(x-f0))).^2;
%g = @(a,b,c,d,x0,xs0,x)((b*a^2./((x-x0).^2+a^2)+c*d^2./((x-xs0).^2+d^2)).^2+(a*b*(x-x0)./((x-x0).^2+a^2)+c*d*(x-xs0)./((x-xs0).^2+d^2)).^2);
hft_ = fittype(t,...
     'dependent',{'y'},'independent',{'x'},...
     'coefficients',{'aa','bd','f0'});
 hst_ = [1,0.5*bndwidth*10^-9,1.0*new_freq(amp0)*10^-9];
 hilo_=[0.01,0.05*0.5*bndwidth*10^-9,0.9500*new_freq(amp0)*10^-9];    
 hiup_=[1,2.1*0.5*bndwidth*10^-9,1.0500*new_freq(amp0)*10^-9];
 hfo_ = fitoptions('method','NonlinearLeastSquares','Lower',hilo_,'Upper',hiup_);
 set(hfo_,'Startpoint',hst_);
[hcf_,gof_,output_] = fit(freqdata, ampdatanorm,hft_,hfo_);
yt=t(hcf_.aa,hcf_.bd,hcf_.f0,freqdata);
yyyyyt=t(hcf_.aa,hcf_.bd,hcf_.f0,arrayfreqdata);

DnoZ_ampdatanorm=10*log10(ampdatanorm);
%DnoZ_ampdatanorm(DnoZ_ampdatanorm==0) = [];

FnoZ_ampdatanorm=10*log10(yt);
%FnoZ_ampdatanorm(FnoZ_ampdatanorm==0) = [];

chi2t =abs(nansum(((DnoZ_ampdatanorm-FnoZ_ampdatanorm).^2)./(FnoZ_ampdatanorm)));
sprintf('%f',chi2t);
chi2redt=chi2t/(length(DnoZ_ampdatanorm)-3);
sprintf('%f',chi2redt);
f0sing=hcf_.f0;
quafL_fit_s=hcf_.f0/(2*hcf_.bd);
ci = confint(hcf_, 0.95);
tt_ = tinv((1+0.95)/2, 3); 
se = (ci(2,:)-ci(1,:)) ./ (2*tt_);
as_err=se(1);
bds_err=se(2);
fs_err=se(3);
stdQL_s=sqrt(((hcf_.f0/(2*hcf_.bd))^2)*((se(3)/(hcf_.f0))^2+(se(2)/(hcf_.bd))^2));
hcf_.f0/(2*hcf_.bd)
%

g = @(ac,bdc,f0c,aL,bdL,f0L,aR,bdR,f0R,x) abs(	ac./(1+1i*(1/bdc)*(x-f0c))+ ... 
                                                aL./(1i*(1/bdL)*(x-f0L))+ ...
                                                aR./(1i*(1/bdR)*(x-f0R)) ).^2;
%g = @(a,b,c,d,x0,xs0,x)((b*a^2./((x-x0).^2+a^2)+c*d^2./((x-xs0).^2+d^2)).^2+(a*b*(x-x0)./((x-x0).^2+a^2)+c*d*(x-xs0)./((x-xs0).^2+d^2)).^2);
hft_ = fittype(g,...
     'dependent',{'y'},'independent',{'x'},...
     'coefficients',{'ac','bdc','f0c','aL','bdL','f0L','aR','bdR','f0R'});
 hst_ = [1,0.5*bndwidth*10^-9,1.0*new_freq(amp0)*10^-9,...
         1,0.5*bndwidth*10^-9,0.9*new_freq(amp0)*10^-9,...
         1,0.5*bndwidth*10^-9,1.1*new_freq(amp0)*10^-9];
 hilo_=[0.01,0.05*0.5*bndwidth*10^-9,0.9500*new_freq(amp0)*10^-9,...
        0.01,0.05*0.5*bndwidth*10^-9,0.8250*new_freq(amp0)*10^-9,...
        0.01,0.05*0.5*bndwidth*10^-9,1.0000*new_freq(amp0)*10^-9];    
 hiup_=[1,2.1*0.5*bndwidth*10^-9,1.0500*new_freq(amp0)*10^-9,...
        1,2.1*0.5*bndwidth*10^-9,1.0000*new_freq(amp0)*10^-9,...
        1,2.1*0.5*bndwidth*10^-9,1.1750*new_freq(amp0)*10^-9];
 hfo_ = fitoptions('method','NonlinearLeastSquares','Lower',hilo_,'Upper',hiup_);
 set(hfo_,'Startpoint',hst_);
[hcf_,gof_,output_] = fit(freqdata, ampdatanorm,hft_,hfo_);


fwhm=2*hcf_.bdc;
quafL_fit=hcf_.f0c/(2*hcf_.bdc);

fwhmL=2*hcf_.bdL;
quafL_fitL=hcf_.f0L/(2*hcf_.bdL);

fwhmR=2*hcf_.bdR;
quafL_fitR=hcf_.f0R/(2*hcf_.bdR);

ci = confint(hcf_, 0.95);
tt_ = tinv((1+0.95)/2, 9); 
se = (ci(2,:)-ci(1,:)) ./ (2*tt_);

%sigma= sqrt(diag(inv((output_.Jacobian).'*output_.Jacobian))); %standard errors

stdQL=sqrt(((hcf_.f0c/(2*hcf_.bdc))^2)*((se(3)/(hcf_.f0c))^2+(se(2)/(hcf_.bdc))^2));


y=g(hcf_.ac,hcf_.bdc,hcf_.f0c,hcf_.aL,hcf_.bdL,hcf_.f0L,hcf_.aR,hcf_.bdR,hcf_.f0R,freqdata);
%y=hcf_.b*(hcf_.a^2./((freqdata-hcf_.x0).^2+hcf_.a^2));
%yy=hcf_.b*(hcf_.a^2./((arrayfreqdata-hcf_.x0).^2+hcf_.a^2));
%yyy=(hcf_.b^2)*(hcf_.a^2./((arrayfreqdata-hcf_.x0).^2+hcf_.a^2)).^2+(hcf_.a*hcf_.b*(arrayfreqdata-hcf_.x0)./((arrayfreqdata-hcf_.x0).^2+hcf_.a^2)+(hcf_.c*hcf_.d./(arrayfreqdata-hcf_.xs0))).^2;
%yyyy=(hcf_.b*hcf_.a^2./((arrayfreqdata-hcf_.x0).^2+hcf_.a^2)+hcf_.c*hcf_.d^2./((arrayfreqdata-hcf_.xs0).^2+hcf_.d^2)).^2+(hcf_.a*hcf_.b*(arrayfreqdata-hcf_.x0)./((arrayfreqdata-hcf_.x0).^2+hcf_.a^2)+hcf_.c*hcf_.d*(arrayfreqdata-hcf_.xs0)./((arrayfreqdata-hcf_.xs0).^2+hcf_.d^2)).^2;
yyyyy=g(hcf_.ac,hcf_.bdc,hcf_.f0c,hcf_.aL,hcf_.bdL,hcf_.f0L,hcf_.aR,hcf_.bdR,hcf_.f0R,arrayfreqdata);

figure(10);
plot(freqdata,10*log10(ampdatanorm),'LineWidth',1.2,'Marker','o','LineStyle','none','MarkerSize',4)
hold on;
plot(arrayfreqdata,10*log10(yyyyyt),'color',[0.9290 0.6940 0.1250],'LineStyle','-','LineWidth',1.2)
plot(arrayfreqdata,10*log10(yyyyy),'r','LineStyle','-','LineWidth',1.2)
scatter(new_freq(amp0)*10^-9,10*log10(ampdatanorm(amp0)),'MarkerEdgeColor',[242, 142, 0]/255, 'MarkerFaceColor',[242, 142, 0]/255,'LineWidth',1.5);
scatter(new_freq(fmnus3dB)*10^-9,10*log10(ampdatanorm(fmnus3dB)),'MarkerEdgeColor',[0, 166, 235]/255, 'MarkerFaceColor',[0, 166, 235]/255,'LineWidth',1.5);
scatter(new_freq(fplus3dB)*10^-9,10*log10(ampdatanorm(fplus3dB)),'MarkerEdgeColor',[217, 123, 244]/255, 'MarkerFaceColor',[217, 123, 244]/255,'LineWidth',1.5);
hold off;
set(gca,'FontSize',25);
grid on;
%xlim([xlimfig(1) xlimfig(2)])
%ylim([-150 0])
xlabel('Frequency / GHz')
%ylabel('Norm. Amplitude [A. U.]')
ylabel('Normalized |S_{21}|^2 / A. U.')
legend_2=legend('Exp. Data','Lorentz single', 'Lorentz w/ 2 tails');
set(legend_2, 'FontSize', 24);

DnoZ_ampdatanorm=10*log10(ampdatanorm);
%DnoZ_ampdatanorm(DnoZ_ampdatanorm==0) = [];

FnoZ_ampdatanorm=10*log10(y);
%FnoZ_ampdatanorm(FnoZ_ampdatanorm==0) = [];

chi2 =abs(nansum(((DnoZ_ampdatanorm-FnoZ_ampdatanorm).^2)./(FnoZ_ampdatanorm)));
sprintf('%f',chi2);
chi2red=chi2/(length(DnoZ_ampdatanorm)-9);
sprintf('%f',chi2red);

sprintf('Exp. Data:    \n[Bandwidth[Hz]      = %s,\n Q_L                = %s,\n f_max[Hz]          = %s],\n A_max [dB]         = %s',...
    num2str(bndwidth) , num2str(quafL),  num2str(new_freq(amp0)),num2str(new_amp(amp0)))

rescent=sprintf('Lorentz Curve:\n[Bandwidth[Hz]      = %s     +/- %s,\n Q_L                = %s       +/- %s,\n f_max[Hz]          = %s +/- %s, \n',...
    num2str(fwhm*10^9),num2str(2*se(2)*10^9),num2str(quafL_fit),num2str(stdQL),num2str(hcf_.f0c*10^9),num2str(se(3)*10^9));
resleft=sprintf('             \n Bandwidth[Hz]      = %s     +/- %s,\n Q_L                = %s       +/- %s,\n f_max[Hz]          = %s +/- %s, \n',...
    num2str(fwhmL*10^9),num2str(2*se(5)*10^9),num2str(quafL_fitL),num2str(stdQL),num2str(hcf_.f0L*10^9),num2str(se(6)*10^9));
resrigh=sprintf('             \n Bandwidth[Hz]      = %s     +/- %s,\n Q_L                = %s       +/- %s,\n f_max[Hz]          = %s +/- %s,\n ',...
    num2str(fwhmR*10^9),num2str(2*se(8)*10^9),num2str(quafL_fitR),num2str(stdQL),num2str(hcf_.f0R*10^9),num2str(se(9)*10^9));
resres=sprintf('             \n A_max[dB]          = %s\n Chi^2_red[dB]      = %s]',...
    num2str(max(abs(10*log10(yyyyy)))), num2str(chi2red));

sprintf('%s %s %s %s',rescent,resleft,resrigh,resres)
hcf_
stringQ3  =sprintf('Q_L(3dB-Sym)=%f',quafL);
stringQ1  =sprintf('Q_L(1dB-Sym)=%f',quafL1);
stringQ1_L=sprintf('Q_L(1dB-Lef)=%f',quafL1_L);
stringQ1_R=sprintf('Q_L(1dB-Rig)=%f',quafL1_R);
stringQLcs =sprintf('Q_L(Lorentz s)=%f +/- %f',quafL_fit_s,stdQL_s);
stringQLc =sprintf('Q_L(Lorentz t)=%f +/- %f',quafL_fit,stdQL);
%msgbox([stringQ3;stringQ1;stringQ1_L;stringQ1_R])
sprintf('%s\n%s\n%s\n%s\n%s\n%s',stringQ3,stringQ1,stringQ1_L,stringQ1_R,stringQLcs,stringQLc)
sprintf(' Raw    f_0(MHz)=%s\n Single f_0(MHz)=%s +/- %s\n Tails  f_0(MHz)=%s +/- %s',num2str(new_freq(amp0)*10^-6,'%0.6f'),num2str(f0sing*10^9*10^-6,'%0.6f'),num2str(fs_err*10^9*10^-6),num2str(hcf_.f0c*10^9*10^-6,'%0.6f'),num2str(se(3)*10^9*10^-6))
sprintf(' Single Chi^2(dB)=%s,\tChi^2_red(dB)=%s\n Tails  Chi^2(dB)=%s,\tChi^2_red(dB)=%s',num2str(chi2t),num2str(chi2redt),num2str(chi2),num2str(chi2red))
end
end