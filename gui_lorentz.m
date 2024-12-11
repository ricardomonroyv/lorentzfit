function gui_lorentz
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
new_freq(amp0)
freqdata=new_freq*10^-9; 
ampdata=10.^(new_amp/10);
ampdatanorm=ampdata/max(ampdata);
%options = optimset('PlotFcns',@optimplotfval);

g = @(a,b,x0,x)(b*(a^2./((x-x0).^2+a^2)));
hft_ = fittype(g,...
     'dependent',{'y'},'independent',{'x'},...
     'coefficients',{'a','b','x0'});
 hst_ = [0.5*bndwidth*10^-9,1,new_freq(amp0)*10^-9];
 hfo_ = fitoptions('method','NonlinearLeastSquares');
 set(hfo_,'Startpoint',hst_);
[hcf_,gof_,output_] = fit(freqdata, ampdatanorm,hft_,hfo_);

fwhm=2*hcf_.a;
quafL_fit=hcf_.x0/(2*hcf_.a);

ci = confint(hcf_, 0.95);
tt_ = tinv((1+0.95)/2, 3); 
se = (ci(2,:)-ci(1,:)) ./ (2*tt_);

%sigma= sqrt(diag(inv((output_.Jacobian).'*output_.Jacobian))); %standard errors

stdQL=sqrt(((hcf_.x0/(2*hcf_.a))^2)*((se(3)/(hcf_.x0))^2+(se(1)/(hcf_.a))^2));

arrayfreqdata=freqdata(1):3*10^-9:freqdata(end);
y=hcf_.b*(hcf_.a^2./((freqdata-hcf_.x0).^2+hcf_.a^2));
yy=hcf_.b*(hcf_.a^2./((arrayfreqdata-hcf_.x0).^2+hcf_.a^2));



figure(10);
plot(freqdata,10*log10(ampdatanorm),'LineWidth',1.2)
hold on;
plot(arrayfreqdata,10*log10(yy),'r','LineStyle','--','LineWidth',1.2)
scatter(new_freq(amp0)*10^-9,10*log10(ampdatanorm(amp0)),'MarkerEdgeColor',[242, 142, 0]/255, 'MarkerFaceColor',[242, 142, 0]/255,'LineWidth',1.5);
scatter(new_freq(fmnus3dB)*10^-9,10*log10(ampdatanorm(fmnus3dB)),'MarkerEdgeColor',[0, 166, 235]/255, 'MarkerFaceColor',[0, 166, 235]/255,'LineWidth',1.5);
scatter(new_freq(fplus3dB)*10^-9,10*log10(ampdatanorm(fplus3dB)),'MarkerEdgeColor',[217, 123, 244]/255, 'MarkerFaceColor',[217, 123, 244]/255,'LineWidth',1.5);
hold off;
set(gca,'FontSize',25);
grid on;
%xlim([xlimfig(1) xlimfig(2)])
%ylim([-150 0])
xlabel('Frequency [GHz]')
ylabel('Norm. Amplitude [A. U.]')
legend_2=legend('Exp. Data','Lorentz curve');
set(legend_2, 'FontSize', 24);

DnoZ_ampdatanorm=10*log10(ampdatanorm);
%DnoZ_ampdatanorm(DnoZ_ampdatanorm==0) = [];

FnoZ_ampdatanorm=10*log10(y);
%FnoZ_ampdatanorm(FnoZ_ampdatanorm==0) = [];

chi2 =abs(nansum(((DnoZ_ampdatanorm-FnoZ_ampdatanorm).^2)./(FnoZ_ampdatanorm)));
sprintf('%f',chi2);
chi2red=chi2/(length(DnoZ_ampdatanorm)-3);
sprintf('%f',chi2red);

sprintf('Exp. Data:    \n[Bandwidth[Hz]      = %s,\n Q_L                = %s,\n f_max[Hz]          = %s]',...
    num2str(bndwidth) , num2str(quafL),  num2str(new_freq(amp0)))

sprintf('Lorentz Curve:\n[Bandwidth[Hz]      = %s     +/- %s,\n Q_L                = %s       +/- %s,\n f_max[Hz]          = %s +/- %s\n Chi^2_red[dB]      = %s]',...
    num2str(fwhm*10^9),num2str(2*se(1)*10^9),num2str(quafL_fit),num2str(stdQL),num2str(hcf_.x0*10^9),num2str(se(3)*10^9), num2str(chi2red))

end
end