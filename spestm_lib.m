classdef spestm_lib
    % When using this library,
    %   x: input with size [N,d], where d dimension N element count
    %   p: AR parameter order and noise subspace lower eigenvector limit
    %   q: MA parameter order
    %   g: Capon model order
    %   a: AR parameters
    %   b: MA parameters
    %   M: Noise subspace upper eigenvector limit
    %   sig2: driving noise variance
    
    properties
        x=[]; % input signal
        fs=[]; % input signal sample rate, if not known will be 1/N
        
        p=12;
        q=2;
        g=2;
        M=[]; %if = [], will be automatic
        Nf=[]; % spectrum frequency point count, if = [], will be automatic
        window_type='rectwin'; % look for https://www.mathworks.com/help/signal/ref/window.html for options
        f_range='full'; % 'full' for [0,2pi] range or 'half' for [0,pi] range
        
        % Constructed parameters
        N=[];
        d=[];
        f=[];
        e=[];
    end
    
    methods
        %% Initialization
        function spestm_obj=init(spestm_obj)
            % Parameter extraction
            spestm_obj.N=size(spestm_obj.x,1);
            spestm_obj.d=size(spestm_obj.x,2);
            
            if isempty(spestm_obj.fs)==1
                spestm_obj.fs=1/spestm_obj.N;
            end
            
            if isempty(spestm_obj.M)==1
                spestm_obj.M=spestm_obj.N;
            end
            
            if isempty(spestm_obj.Nf)==1
                spestm_obj.Nf=2^nextpow2(spestm_obj.N);
            end

            if strcmp(spestm_obj.f_range,'full')==1
                spestm_obj.f=linspace(0,1.0*(1-1/spestm_obj.Nf),spestm_obj.Nf);
            else
                spestm_obj.f=linspace(0,0.5*(1-1/spestm_obj.Nf),spestm_obj.Nf);
            end
                        
            % e Matrix
            k = 0:(spestm_obj.M-1);
            spestm_obj.e=zeros(spestm_obj.M,spestm_obj.Nf);
            for n = 1:spestm_obj.M
                spestm_obj.e(n,:) = exp(1i*2*pi*k(n)*spestm_obj.f);
            end
        end
        
        %% PSD Functions
        function [y,F]=psd_periodogram(spestm_obj)
            % Window generation
            eval(['w=window(@', spestm_obj.window_type, ',', num2str(spestm_obj.N) ');']);
            
            % PSD estimation
            xw=spestm_obj.x.*w;
            rxx=cell(1,spestm_obj.d);
            for i=1:spestm_obj.d
                rxx{i}=xcorr(xw(:,i),xw(:,i))/spestm_obj.fs;
                rxx{i}=rxx{i}(spestm_obj.N:end);
            end

            if strcmp(spestm_obj.f_range,'full')==1
                y=zeros(spestm_obj.Nf,spestm_obj.d);
                for i=1:spestm_obj.d
                    y(:,i)=abs(fft(rxx{i},spestm_obj.Nf,1));
                end
            else
                y=zeros(2*spestm_obj.Nf,spestm_obj.d);
                for i=1:spestm_obj.d
                    y(:,i)=abs(fft(rxx{i},2*spestm_obj.Nf,1));
                end
                y=y(1:spestm_obj.Nf,:);
            end
            F=spestm_obj.f;
        end
        
        function [y,F]=psd_blackmantukey(spestm_obj)  
            % Window generation
            eval(['w=window(@', spestm_obj.window_type, ',', num2str(2*spestm_obj.N) ');']);
            w=w((spestm_obj.N+1):end);
            
            % PSD estimation
            rxx=cell(1,spestm_obj.d);
            for i=1:spestm_obj.d
                rxx{i}=xcorr(spestm_obj.x(:,i),spestm_obj.x(:,i))/spestm_obj.fs;
                rxx{i}=w.*rxx{i}(spestm_obj.N:end);
            end

            if strcmp(spestm_obj.f_range,'full')==1
                y=zeros(spestm_obj.Nf,spestm_obj.d);
                for i=1:spestm_obj.d
                    y(:,i)=abs(fft(rxx{i},spestm_obj.Nf,1));
                end
            else
                y=zeros(2*spestm_obj.Nf,spestm_obj.d);
                for i=1:spestm_obj.d
                    y(:,i)=abs(fft(rxx{i},2*spestm_obj.Nf,1));
                end
                y=y(1:spestm_obj.Nf,:);
            end
            F=spestm_obj.f*spestm_obj.fs;
        end
        
        function [y,F] = psd_aryulewalker(spestm_obj)
            % Window generation
            eval(['w=window(@', spestm_obj.window_type, ',', num2str(spestm_obj.N) ');']);
            
            % Parameter estimation
            xw=spestm_obj.x.*w;
            a=zeros(spestm_obj.p+1,spestm_obj.d);
            sig2=zeros(spestm_obj.d,1);
            for i=1:spestm_obj.d
                [a(:,i),sig2(i)]=spestm_obj.arp_yulewalker(xw(:,i));
            end

            % PSD estimation
            y=zeros(spestm_obj.Nf,spestm_obj.d);
            for i=1:spestm_obj.d
                y(:,i)=sig2(i)./(eps+(abs(ctranspose(spestm_obj.e(1:(spestm_obj.p+1),:))*a(:,i)).^2)');
            end
            F=spestm_obj.f;
        end
        
        function [y,F] = psd_capon(spestm_obj)
            % Window generation
            eval(['w=window(@', spestm_obj.window_type, ',', num2str(spestm_obj.N) ');']);

            y=zeros(spestm_obj.Nf,spestm_obj.d);
            for i=1:spestm_obj.d
                % Autocorrelation matrix
                xw=spestm_obj.x(:,i).*w;
                Rxx = spestm_obj.side_autocorrmat(xw,0,spestm_obj.p);
                
                Rxx_inv=pinv(Rxx);

                Rxx_inv_g_1=power(Rxx_inv,spestm_obj.g-1);
                Rxx_inv_g=Rxx_inv*Rxx_inv_g_1;
                
                % PSD Estimation
                for n = 1:spestm_obj.Nf
                    econj=conj(spestm_obj.e(1:(spestm_obj.p+1),n));

                    y(n,i) = (econj'*Rxx_inv_g_1*econj)./(econj'*Rxx_inv_g*econj);
                end
            end
            y=abs(y);
            F=spestm_obj.f;
        end
        
        function [y,F] = psd_armodcov(spestm_obj)
            % Window generation
            eval(['w=window(@', spestm_obj.window_type, ',', num2str(spestm_obj.N) ');']);
            
            % Parameter estimation
            xw=spestm_obj.x.*w;
            a=cell(1,spestm_obj.d);
            sig2=zeros(spestm_obj.d,1);
            for i=1:spestm_obj.d
                [a{i},sig2(i)]=spestm_obj.arp_modcov(xw(:,i));
            end

            % PSD estimation
            y=zeros(spestm_obj.Nf,spestm_obj.d);
            for i=1:spestm_obj.d
                y(:,i)=sig2(i)./(eps+(abs(ctranspose(spestm_obj.e(1:size(a{i},1),:))*a{i}).^2));
            end
            F=spestm_obj.f;
        end
        
        function [y,F] = psd_arburg(spestm_obj)
            % Window generation
            eval(['w=window(@', spestm_obj.window_type, ',', num2str(spestm_obj.N) ');']);
            
            % Parameter estimation
            xw=spestm_obj.x.*w;
            a=zeros(spestm_obj.p+1,spestm_obj.d);
            sig2=zeros(spestm_obj.d,1);
            for i=1:spestm_obj.d
                [a(:,i),sig2(i)]=spestm_obj.arp_burg(xw(:,i));
            end

            % PSD estimation
            y=zeros(spestm_obj.Nf,spestm_obj.d);
            for i=1:spestm_obj.d
                y(:,i)=sig2(i)./(eps+(abs(ctranspose(spestm_obj.e(1:(spestm_obj.p+1),:))*a(:,i)).^2)');
            end
            F=spestm_obj.f;
        end
        
        function [y,F] = psd_arma(spestm_obj)
            % Window generation
            eval(['w=window(@', spestm_obj.window_type, ',', num2str(spestm_obj.N) ');']);
            
            % Parameter estimation
            xw=spestm_obj.x.*w;
            a=zeros(spestm_obj.p+1,spestm_obj.d);
            b=zeros(spestm_obj.q+1,spestm_obj.d);
            sig2=zeros(spestm_obj.d,1);
            for i=1:spestm_obj.d
                [b(:,i),a(:,i),sig2(i)]=spestm_obj.arp_modyulewalker(xw(:,i));
            end

            % PSD estimation
            y=zeros(spestm_obj.Nf,spestm_obj.d);
            for i=1:spestm_obj.d
                y_nom=(eps+(abs(ctranspose(spestm_obj.e(1:(spestm_obj.q+1),:))*b(:,i)).^2)');
                y_den=sig2(i)./(eps+(abs(ctranspose(spestm_obj.e(1:(spestm_obj.p+1),:))*a(:,i)).^2)');
                y(:,i)=y_nom.*y_den;
            end
            F=spestm_obj.f;
        end

        function [y,F] = psd_music(spestm_obj)
            % Window generation
            eval(['w=window(@', spestm_obj.window_type, ',', num2str(spestm_obj.N) ');']);

            y=zeros(spestm_obj.Nf,spestm_obj.d);
            for i=1:spestm_obj.d
                % Autocorrelation matrix
                xw=spestm_obj.x(:,i).*w;
                Rxx=spestm_obj.side_autocorrmat(xw,0,spestm_obj.M-1);

                % Eigendecomposition
                [v,lambda]=eig(Rxx);
                lambda=diag(lambda);
                [~, c] = sort(lambda,'descend');

                % PSD Estimation
                p_noise = 0;
                for j = (spestm_obj.p+1):spestm_obj.M
                    ev = (spestm_obj.e.')*(v(:,c(j)));
                    p_noise = p_noise + abs(ev).^2;
                end
                y(:,i) = 1./(eps+p_noise');
            end
            F=spestm_obj.f*spestm_obj.fs;
        end

        function [y,F]=psd_minnorm(spestm_obj)    
            % Window generation
            eval(['w=window(@', spestm_obj.window_type, ',', num2str(spestm_obj.N) ');']);
            
            % Parameter estimation
            xw=spestm_obj.x.*w;
            a=zeros(spestm_obj.M,spestm_obj.d);
            sig2=zeros(spestm_obj.d,1);
            for i=1:spestm_obj.d
                [a(:,i),sig2(i)]=spestm_obj.arp_minnorm(xw(:,i),spestm_obj.M-1);
            end

            % PSD estimation
            y=zeros(spestm_obj.Nf,spestm_obj.d);
            for i=1:spestm_obj.d
                y(:,i)=sig2(i)./(eps+(abs(ctranspose(spestm_obj.e)*a(:,i)).^2)');
            end
            F=spestm_obj.f*spestm_obj.fs;
        end

        %% AR Parameter Estimation Functions
        function [a,sig2]=arp_yulewalker(spestm_obj,x)
            % Autocorrelation Matrix
            Rxx=spestm_obj.side_autocorrmat(x,0,spestm_obj.p);
            
            % Parameter a Estimation
            A=Rxx(1:spestm_obj.p,1:spestm_obj.p);
            detA=det(A);
            b=Rxx(2:(spestm_obj.p+1),1);
            a=ones(spestm_obj.p+1,1);
            for i=1:spestm_obj.p
                A_k = A;
                A_k(:,i) = -b;
                a(i+1)=det(A_k)./detA;
            end
            
            % Parameter sig2 Estimation
            sig2=sum(a'.*Rxx(1,:));
        end
        
        function [a,sig2]=arp_modcov(spestm_obj,x)
            % Autocorrelation Matrix
            Cxx=spestm_obj.side_modcovmat(x,0,spestm_obj.p);
            
            % Parameter a Estimation
            p_new=rank(Cxx);
            if p_new~=(spestm_obj.p+1)
                Cxx=spestm_obj.side_modcovmat(x,0,spestm_obj.p+1);
                display(['Rank deficient, p is changed to ' num2str(p_new)]);
            end

            A=Cxx(2:(p_new),2:(p_new));
            detA=det(A);
            b=Cxx(2:(p_new),1);
            a=ones(p_new,1);
            for i=1:(p_new-1)
                A_k = A;
                A_k(:,i) = -b;
                a(i+1)=det(A_k)./detA;
            end
            
            % Parameter sig2 Estimation
            sig2=sum(a'.*Cxx(1,1:p_new));
        end

        function [a,sig2]=arp_burg(spestm_obj,x)
            [a,sig2]=arburg(x,spestm_obj.p);
        end
        
        function [b,a,sig2]=arp_modyulewalker(spestm_obj,x)
            % STEP 1
            % Autocorrelation Matrix
            Rxx=spestm_obj.side_autocorrmat(x,spestm_obj.q,spestm_obj.p);
            
            % Parameter a Estimation
            A=Rxx(1:spestm_obj.p,1:spestm_obj.p);
            detA=det(A);
            b=Rxx(2:(spestm_obj.p+1),1);
            a=ones(spestm_obj.p+1,1);
            for i=1:spestm_obj.p
                A_k = A;
                A_k(:,i) = -b;
                a(i+1)=det(A_k)./detA;
            end
                        
            % STEP 2
            % Filter
            x_f=filter(a,1,x);
            % Autocorrelation Matrix
            Rxx=spestm_obj.side_autocorrmat(x_f,0,spestm_obj.p);
            
            % Parameter b_a Estimation
            A=Rxx(1:spestm_obj.p,1:spestm_obj.p);
            detA=det(A);
            b=Rxx(2:(spestm_obj.p+1),1);
            b_a=ones(spestm_obj.p+1,1);
            for i=1:spestm_obj.p
                A_k = A;
                A_k(:,i) = -b;
                b_a(i+1)=det(A_k)./detA;
            end
            
            % Parameter sig2 Estimation
            sig2=sum(b_a'.*Rxx(1,:));
            
            % STEP 3
            % Autocorrelation Matrix
            Raa=spestm_obj.side_autocorrmat(a,0,spestm_obj.q);
            
            % Parameter b_a Estimation
            A=Raa(1:spestm_obj.q,1:spestm_obj.q);
            detA=det(A);
            b=Raa(2:(spestm_obj.q+1),1);
            b_b=ones(spestm_obj.q+1,1);
            for i=1:spestm_obj.q
                A_k = A;
                A_k(:,i) = -b;
                b_b(i+1)=det(A_k)./detA;
            end
            
            b=b_b;
        end
        
        function [a,sig2]=arp_minnorm(spestm_obj,x,M)
            % Autocorrelation Matrix
            Rxx=spestm_obj.side_autocorrmat(x,0,M);
            
            % Eigendecomposition
            [v,lambda]=eig(Rxx);
            lambda=diag(lambda);
            [~, c] = sort(lambda,'descend');
            v=v(:,c);

            % Parameter a Estimation
            u1=zeros(M+1,1); u1(1)=1;

            Vn=v(:,(spestm_obj.p+1):M);
            Pn=Vn*ctranspose(Vn);
            lambda=ctranspose(u1)*Pn*u1;

            a=(Pn*u1)./(eps+lambda);

            % Parameter sig2 Estimation
            sig2=sum(a'.*Rxx(1,:));
        end

        %% Side Functions
        function xn = side_awgn(spestm_obj,SNR)
            sigEner = sqrt(sum((spestm_obj.x).^2,1));                    % energy of the signal
            noiseEner = sigEner./(10^(SNR/10));        % energy of noise to be added
            noiseVar = noiseEner./(size(spestm_obj.x,1)-1);     % variance of noise to be added
            noiseStd = sqrt(noiseVar);                   % std. deviation of noise to be added
            noise = noiseStd.*randn(size(spestm_obj.x));           % noise
            xn = spestm_obj.x+noise;                        % noisy signal
        end
        
        function Rxx=side_autocorrmat(spestm_obj,x,pshift,M)
            Nx=size(x,1);
            x=x-mean(x,1);
            rxx=xcorr(x,x)/spestm_obj.fs;
            rxx=circshift(rxx,-pshift);
            Rxx=zeros(M+1,M+1);

            for j=1:(M+1)
                Rxx(j,j:(M+1))=rxx(Nx:(Nx+M+1-j));
                Rxx(j,1:(j-1))=conj(rxx((Nx-j+1):(Nx-1)));
            end
        end
        
        function Cxx=side_modcovmat(spestm_obj,x,pshift,M)           
            Cxx = zeros( M+1, M+1  );
            
            xconj=conj(x);
            for j=0:M
                for k=0:M
                    kk=mod(k+pshift,M+1);
                    n1=(M+1):spestm_obj.N;
                    n2=1:(spestm_obj.N-M);
                    Cxx(j+1,k+1) = sum(xconj(n1-j).*x(n1-kk)) + sum(x(n2+j).*xconj(n2+kk));
                end
            end
            
            Cxx=Cxx./(spestm_obj.fs);
        end
    end
end

