list=getfn('D:\Promotion\Daten\SpinningDisk\Eventanalysis','judge.mat');
list=list';
%%
for k=1:size(list,1)
load(list{k});
dmy=list{k};
dmy=dmy(1:end-9)
    for i=1:size(jstr,1)
        if strcmp(jstr{i},'upice')
            load([dmy,'results.mat']);
            figure(1)
            dmyl=results.Trace_l(:,i)-results.Trace_BG_l-(0.146*(results.Trace_r(:,i)-results.Trace_BG_r));
            plot(dmyl)
            hold on
            plot(movmean(dmyl,50))
            hold off
            figure(2)
            plot(results.Trace_r(:,i)-results.Trace_BG_r)
            hold on
            plot(movmean(results.Trace_r(:,i)-results.Trace_BG_r,50))
            hold off
            drawnow            
            str=input('(n)ext, (e)mpty [upice], (c)ontent remained [upicc]','s');
            switch str
                case 'n'                    
                case 'e'
                    jstr{i}='upice';
                case 'c'
                    jstr{i}='upicc';
            end            
            %Eingabe von String addition mit Rückspeicherung
        end        
    end
end

