% In this test, we verify the CPQR by approximating RatioCut for arbitraty
% k unequal sized-clusters 
clear
close all
clc
rng(1)
k = 3;
% m =  [70 80 90 100 110 120 130];
%170 200 230 280
m = [300 300 300]; 
n = sum(m);
noise = 0;
Na = 20;
Nb = 5;
success = zeros(Na,Nb);
% a_grid =  zeros(Na,Nb);
% b_grid =  zeros(Na,Nb);
% p_grid =  zeros(Na,Nb);
% q_grid =  zeros(Na,Nb);
a_range = linspace(2,20,Na);
b_range = linspace(1,10,Nb);
% a = 90;
% b= 1;
trials = 9;
for t = 1:trials
    t
    a_count = 0;
    for a = a_range
        a_count = a_count+1;
        b_count = 0;
        for b = b_range
            b_count = b_count+1;
            
            if k == 2
                p = a*log(n)/n;
                q = b*log(n)/n;
            else
                p = a*log(m(1))/m(1);
                q = b*log(m(1))/m(1);
            end
           
            % only when this makes sense (just for testing faster)
            if p > q
                valid = 0;
                while ~valid
                    A = tril(rand(n) < q);
                    truth = [];
                    idx_start = 1;
                    for j = 1:k
                        idx_end = idx_start + m(j)-1;
                        idx = idx_start:idx_end;                
                        A(idx,idx) = tril(rand(m(j)) < p);
                        truth = [truth j*ones(1,m(j))];
                        idx_start = idx_end+1;
                    end
                 
                    A = A + A';
                    
                    if noise > 0
                        An = tril(rand(n)<noise);
                        An = An + An';
                        A(A==0) = An(A==0);
                    end
                    
                    
                    A = A - diag(diag(A));
                    
                    
                    
                    d = sum(A);
                    if min(d) > 0
                        valid = 1;
                    end
                end
                
                Ds = diag(sqrt(1./d));
                Lhat = Ds*A*Ds;

                
                [Vhat e] = eigs(Lhat,k);
                
                if norm(Vhat'*Vhat - eye(k)) > 1e-12
                    disp('eigs issue');
                    [Vhat e] = eig(Lhat);
                    Vhat = Vhat(:,n-k+1:end);
                end
                
                
                [Vloc piv] = lrcol_rand(Vhat,k,1,5);
                [throw, set] = max(abs(Vloc'));
%                 sucess = 0; 
                if verify_true(set,truth)
%                     sucess = 1; 
                    success(a_count,b_count) = success(a_count,b_count)+1;
                end
                
       
                % Now we have a similar matrix A, next build the graph
                % Laplacian matrix based on L = D - A; 
                vol = sum(d);
                D = diag(d);
                L = D-A; 
                [V,E] = eig(L);
%               multiplier = 1/V(1,1);
%               V = multiplier * V;

                % Now do the sanity checks for L to see if it satisfies
                % othe properties of graph Laplacian.
               
                
                % Normalize the graph Laplacian with the random walk
                % format.
                Lrw = D^(-1)*L;
                [Vrw,Erw] = eig(Lrw);
                
                % Want to minimize RatioCut/Ncut
                % First do RatioCut
                % For k > 2 problem, we are doing the
                % optimization problem indicated by min_H Tr(H'LH) subject
                % to H'H = I where the columns of H (n*k) are the k indicator 
                % vectors h_j = (h_{1,j},...h_{n,j}) where h_{i,j} are defined 
                % as in the equation 5. 
                % By closed form solution, we know the optimal H contains
                % the first k eigenvectors of L as columns. We can see that
                % the matrix H is in fact the matrix U used in the
                % unnormalized spectral clustering algo. 
                
                H = V(:,1:k);
                scheckOrthog = H'*H -eye(k,k);
                idx = kmeans(H,k);
                % Construct answers of similar format from CPQR solution.
                % Notice in this case different clusters have the same
                % size (namely 150)
                [row,col] = size(H);
                Hcpqr = ones(row,col);
                for i= 1:row
                    for j = 1:col
                        if(set(i)==j) 
                                Hcpqr(i,j) = 1/sqrt(m(j));
%                             end
                        else
                            Hcpqr(i,j) = 0;
                        end 
                    end 
                end 
%                 scheck2= Hcpqr' * Hcpqr - eye(k,k)
                trace1 = trace(H'*L*H);
                trace2 = trace(Hcpqr'*L*Hcpqr);
                FINAL = (trace1 - trace2)/trace1
            end
            a_grid(a_count,b_count) = a;
            b_grid(a_count,b_count) = b;
            p_grid(a_count,b_count) = p;
            q_grid(a_count,b_count) = q;
    end       
end 
                % Now check the difference between this closed
                % form solution and the one computed by our CPQR algo. 
               
%                 group1diff = length(groupA)-length(groupA);
%                 
%                 group2diff = length(groupB)-length(groupB);
%                 
                % Next do Ncut (the optimization is Eq(9) in the paper)
                % where the solution is achieved at the second eignevector
                % of Lrw
                frw = Vrw(:,2);
                checkOrthogRW = sum(D*frw.*ones(n,1));
                diffRW= frw'*D*frw-vol;
                

end

successRate = success/trials;
% success2 = success2/trials;
% fiedler = fiedler/trials;
% qr_kmeans = qr_kmeans/trials;
% comp_true = comp_true/trials;
% MLE = MLE/trials;
% theory(theory<1) = 0;
% theory(theory>1) = 1;
% 
% save backup_block_test_unequal.mat a_grid b_grid success fiedler success2 qr_kmeans theory qr_kmeans p_grid q_grid k m n trials comp_true MLE
% fname = [num2str(k) '_block_test_unequal.mat'];
% save(fname,'a_grid','b_grid','success','fiedler','comp_true','MLE',...
%     'success2','qr_kmeans','theory','qr_kmeans','p_grid','q_grid','k','m','n','trials');
save testEasy.mat k m n L Lrw A 
fname = [num2str(k) '_block_test_sherlock.mat'];
