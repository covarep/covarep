%SINVEST1 shows the parameters of a time series calculated by INVEST1
% only called by INVEST1

%       $Id: sinvest1.m 5090 2008-06-05 08:12:04Z schloegl $
%	Copyright (C) 1998-2002,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Fs=flag_implicit_samplerate;
M=size(AutoCorr,2);
oo=M-1;
while 1,
        K=menu('Select','Autocovariance ACOVF','Autocorrelation ACF', ...
                'Partial ACF PACF','Coeff. of Determination R²', ...
                'Error curve',...
                'Autoregressive Parameters',...
                'Information Criteria', ...
                'Matched Filter', ... 
                'Log PSD and Phase', ... 
                'Poles', ... 
                'Inverse Filtering', ... 
                'Spectra H(f,p)', ... 
                'Entropy H=ln(det(R))', ... 
                'Histogram of MOPS', ...
                'end');
        subplot(111);	  
        if K==1
                plot(0:M,AutoCov);
                title('Autocovariance function ACOVF(k)');
                xlabel('Lag k');
        elseif K==2
                if exist('OCTAVE_VERSION')==5             %%%%% Fuer OCTAVE 
                        plot(1:M,AutoCorr);
                elseif strcmp(version,'MIDEVA')           %%%%% fuer MatCom
                        plot(1:M,AutoCorr);
                else                                      %%%%% fuer Matlab
                        %               size(AutoCorr),size(ACFsd)]
                        %                errorbar(ones(nr,1)*(1:M),AutoCorr,ACFsd,ACFsd);
                        %               errorbar(1:M,AutoCorr,AutoCorr-ACFsd,AutoCorr+ACFsd);
                        plot(1:M,AutoCorr,'b',[1,M],[-1,1;-1,1]/sqrt(min(NC)),'b:');
                        if exist('OCTAVE_VERSION')<5
                                legend({'ACF','1/sqrt(N)'})
                        end;
                end;
                title('Autocorrelation function ACF(k)');
                xlabel('Lag k');
                
        elseif K==3
                %rc=ARPMX(:,(1:M).*(2:M+1)/2);
                %plot(1:M,PartACF);
                if size(ARPMX,2)==2*Pmax,
                        RC=ARPMX(:,Pmax+1:2*Pmax);
                else
                        RC=ARPMX(:,(1:M).*(2:M+1)/2);
                end;
                % according to http://www.itl.nist.gov/div898/handbook/pmc/section4/pmc4463.htm
                % is the 95% confidence interval 2/sqrt(N) 
                %plot(1:M,RC,'b',[1,M],[3,3;2,2;1,1;-1,-1;-2,-2;-3,-3]*(1/sqrt(min(NC))),'b:');
                %legend({'Part. ACF','1/sqrt(N)','2/sqrt(N) = 95% confidence interval','3/sqrt(N)'})
                plot(1:M,RC,'b',[1,M],[2,2;-2,-2]'*(1/sqrt(min(NC))),'b:');
                if exist('OCTAVE_VERSION')<5
                        legend({'Part. ACF','2/sqrt(N) = 95% confidence interval'})
                end;	
                title('Partial Autocorrelation function PACF(k)');
                xlabel('Lag k');
        elseif K==4
                plot(0:M,(E(1)-E)/E(1));
                title('Determination of Regression R²=1-var{E}/var{Y}');
                xlabel('Model order p');
        elseif K==5
                plot(0:M,E,'r');
                if exist('OCTAVE_VERSION')<5
                        v=axis; v(3)=min([v(3); 0]); axis(v);
                end
                title('Mean Square (prediction) Error decrease with increasing model order');
                xlabel('Model order p');
        elseif K==6 
                plot(1:oo,ARPMX(:,oo/2*(oo-1)+(1:oo))','o-');
                title( ['AutoRegressive Parameters with model order' int2str(oo)])
        elseif K==7 
                while 1
                        K=menu('Select Criterion',...
                                'Final Predection Error FPE', ...
                                'Akaike Information Criterion AIC', ... 
                                'Bayesian Akaike Information Criterion BIC', ... 
                                'Schwartz''s Bayesian Crit. SBC (=MDL)', ... 
                                'Parzen''s CAT Criterion ', ... 
                                'PHI Criterion', ... 
                                'end');
                        
                        if K==1
                                tmp=min(2*optFPE+2,M);                           
                                oo=optFPE;
                                plot(0:tmp-1,FPE(1:tmp),'r',optFPE,FPE(optFPE+1),'ro');
                                if exist('OCTAVE_VERSION')<5
                                        text(optFPE,FPE(optFPE+1),sprintf('%i',optFPE));
                                        v=axis; v(3)=min([v(3); 0]); axis(v);
                                end;
                                title('Final Prediction Error FPE criterion');
                        elseif K==2 
                                tmp=min(2*optAIC+2,M);
                                oo=optAIC;
                                plot(0:tmp-1,AIC(1:tmp),'r',optAIC,AIC(optAIC+1),'ro');
                                if exist('OCTAVE_VERSION')<5
                                        text(optAIC,AIC(optAIC+1),sprintf('%i',optAIC));
                                        v=axis; v(3)=min([v(3); 0]); axis(v);
                                end;
                                title('Akaike''s Information Criterion AIC');
                        elseif K==3
                                tmp=min(2*optBIC+2,M);
                                oo=optBIC;
                                plot(0:tmp-1,BIC(1:tmp),'r',optBIC,BIC(optBIC+1),'ro');
                                if exist('OCTAVE_VERSION')<5
                                        text(optBIC,BIC(optBIC+1),sprintf('%i',optBIC));
                                        v=axis; v(3)=min([v(3); 0]); axis(v);
                                end;
                                title('Bayesian Akaike Information Criterion BIC');
                        elseif K==4
                                tmp=min(2*optSBC+2,M);
                                oo=optSBC;
                                plot(0:tmp-1,SBC(1:tmp),'r',optSBC,SBC(optSBC+1),'ro');
                                if exist('OCTAVE_VERSION')<5
                                        text(optSBC,SBC(optSBC+1),sprintf('%i',optSBC));
                                        v=axis; v(3)=min([v(3); 0]); axis(v);
                                end;
                                title('Schwartz''s Bayesian Criterion SBC');
                                %elseif K==5
                                %	tmp=min(2*optMDL,M);
                                %	oo=optMDL;
                                %	plot(0:tmp-1,MDL(1:tmp),'r',optMDL,MDL(optMDL+1),'ro');
                                %	v=axis; v(3)=min([v(3); 0]); axis(v);
                                %	text(optMDL,MDL(optMDL+1),sprintf('%i',optMDL))
                                %	title('Minimal Description length Criterion MDL');
                        elseif K==5
                                tmp=min(2*optCAT+2,M);
                                oo=optCAT;
                                plot(0:tmp-1,CATcrit(1:tmp),'r',optCAT,CATcrit(optCAT+1),'ro');
                                if exist('OCTAVE_VERSION')<5
                                        text(optCAT,CATcrit(optCAT+1),sprintf('%i',optCAT));
                                        v=axis; v(3)=min([v(3); 0]); axis(v);
                                end;
                                title('Parzen''s CAT Criterion ');
                        elseif K==6
                                tmp=min(2*optPHI+2,M);
                                oo=optPHI;
                                plot(0:tmp-1,PHI(1:tmp),'r',optPHI,PHI(optPHI+1),'ro');
                                if exist('OCTAVE_VERSION')<5
                                        text(optPHI,PHI(optPHI+1),sprintf('%i',optPHI));
                                        v=axis; v(3)=min([v(3); 0]); axis(v);
                                end;
                                title('Phi criterion ');   
                        elseif K==7
                                %[ARP,rc,res] =durlev(sum(AutoCov(:,1:(oo+1)),1));
                                ARP=ARPMX(:,oo/2*(oo-1)+(1:oo));
                                break;	
                        end;	%IF
                end;	%WHILE
        elseif K==8
                % if model order p is given then the filter parameters A are
                % A=[1; -earpyw(signal,p)]
                % inverse filtering is    invfiltsignal=filter(A,1,signal);
                
                h=zeros(nr,512)';
                w=zeros(nr,512)';
                for k=1:nr,
                        tmp=freqz(sqrt(E(k,oo+1)),[1 -ARPMX(k,oo/2*(oo-1)+(1:oo))],512);
                        h(:,k)=tmp(:);
                end
                
                plot((0:511)/512/2*Fs,abs(h));
                %plot((0:511)/512/2,abs(freqz(1,[1 -ARP]')));
                title('Matched Filter');
                xlabel('Frequency f');
                ylabel('|H(f)|');
        elseif K==9
                h=zeros(nr,512)';
                w=zeros(nr,512)';
                for k=1:nr,
                        tmp=freqz(sqrt(E(k,oo+1)),[1 -ARPMX(k,oo/2*(oo-1)+(1:oo))],512);
                        h(:,k)=tmp(:);
                end;        
                subplot(211);
                semilogy((0:511)/512/2*Fs,abs(h'))
                %semilogy((0:511)/512/2,abs(freqz(1,[1 -ARP]')))
                title('Logarithmic Spectral Density Fct.');
                subplot(212);
                plot((0:511)/512/2*Fs,angle(h)');
                %plot((0:511)/512/2,angle(freqz(1,[1 -ARP]')));
                ylabel('rad');
        elseif K==10           
                %	clf;
                %r = roots([1 -ARP]);
                r = roots([1 -ARPMX(:,oo/2*(oo-1)+(1:oo))]);
                t = 0:1/70:2*pi;
                plot(cos(t), sin(t), 'b:',real(r), imag(r), 'rx');
                
                %	zplane([],[1 -ARPmx(oo+1,1:oo)]);
                title('Pole Diagram');
                xlabel('real(z)');
                ylabel('imag(z)');
                
                MATLAB_VERSION = version;
                if MATLAB_VERSION(1)=='4'
                        ax = gca;
                        tmp = get(ax,'Aspect');
                        set(ax,'Aspect',[tmp(1),1]);
                elseif MATLAB_VERSION(1)=='5'
                        ax = gca;
                        tmp = get(ax,'DataAspectRatio');
                        set(ax,'PlotBoxAspectRatio',tmp);
                end;
        elseif K==11
                plot([Y(:) filter([1 -ARPMX(:,oo/2*(oo-1)+(1:oo))],1,Y(:))-max(Y)+min(Y)]);
        elseif K==12
                %[tmp,ARPmx,PE]= acf2pacf(AutoCov(2:length(AutoCov))/AutoCov(1),AutoCov(1));
                %[arp,rc,PE,ARPMX] = durlev(AutoCov);
                %[tmp,ARPmx]=arp2pacf(AR);
                N = Pmax-1; %2*oo-1; %2^(ceil(log(oo)/log(2)));
                sdf = zeros(512,N); %length(AR));
                for k = 1:N; %[k,size(sdf),N],%length(AR);
                        % [sdf(:,k),F] = freqz(1,[1 -ARPmx(k,1:k)],N);
                        [tmp,F] = freqz(1,[1 -ARPMX(1,k/2*(k+1)+(1:k))]',512);
                        sdf(:,k)=tmp(:);
                        % sdf(:,k)=sqrt(E(k+1))*sdf(:,k);
                end;
                mesh(F/2/pi*Fs,1:N,log10(abs(sdf)'));
                zlabel('log10 |H(f,p)|');ylabel('model order p'); xlabel('frequency f [2*pi rad/s]');
                if exist('OCTAVE_VERSION')<5
                        view(30,45);
                end;
        elseif K==13
                for k=1:M,
                        xxx=eig(toeplitz(AutoCov(1:k)/E(k)));
                        H(k) = .5*sum(log(xxx));
                        H1(k)= H(k)/(k);
                        %xxx1=eig(toeplitz(AutoCov(1:k)/E(k)*E(1)));
                        %xxx1=eig(toeplitz(AutoCorr(1:k)));
                        %H2(k) =.5*sum(log(xxx1));
                        %H3(k)=.5*sum(log(xxx1))/(k);
                        %if any(xxx<=0) fprintf(1,'COV positive definite for p up to %i\n',k-1); break; end;
                end;
                if 0 
                        subplot(211);
                        plot(0:k-1,H2');
                        title('Entropy H(k) depending on number of coefficients')
                        subplot(212);
                        plot(0:k-1,H3');
                        title('Entropy rate H(k) depending on number of coefficients')
                else
                        subplot(311);
                        plot(0:k-1,H');
                        title('Entropy H(k) depending on number of coefficients')
                        subplot(312);
                        xxx=(-diff(log(E(:))))/2; %[size(xxx)size(H) M k]
                        plot(0:k-1,H1','b',1:k,xxx(1:k),'g');
                        %plot(0:k-1,H1','b',0:k-1,log(E(1:k))-log(E(1)),'b',1:k,xxx(1:k),'g');
                        %plot(0:k-1,H1','b',1:k,xxx,'g',1:k-1,cumprod(1+xxx(1:M-1))./H(2:M)','r');
                        title('Entropy rate H(k) and Entropy difference H(k)-H(k+1)')
                        subplot(313);
                        plot(1:k,log(xxx(1:k)),'g');
                        title('LOG Entropy difference for p->p+1: H(k)-H(k+1)')
                end;
        elseif K==14,
                tmp = histo3(MOPS(1:max(1,size(MOPS,1)-1),:));
                if exist('OCTAVE_VERSION')<5,
                        bar(tmp.X,tmp.H,'stacked');
                else
                        bar(tmp.X,sum(tmp.H,2));
                end;
                xlabel('model order p')
                if exist('OCTAVE_VERSION')<5
                        %legend({'FPE','AIC','BIC','SBC','MDL','CAT','PHI','JEW','HAR'});
                        legend('FPE','AIC','BIC','SBC','MDL','CAT','PHI','JEW','HAR');
                end;        
        elseif K==15
                break; 
        end;
end;
