% In this test, we verify the CPQR by approximating RatioCut for arbitraty
% k with unequal sized-clusters 
clear
close all
clc
rng(1)
k = 3;
m = [120 140 170]; 
n = sum(m);
noise = 0;
Na = 20;
success = zeros(Na,1);
% a = 90;
b= 5;
a_range = linspace(5,100,Na);
rel_error = zeros(Na,1);
trials = 9;
a_count = 0;
for a = a_range
        a_count = a_count+1;
     for t = 1:trials
            t
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
%                 Lhat = Ds*A*Ds;                
                [Vhat e] = eigs(A,k);
                
%                 if norm(Vhat'*Vhat - eye(k)) > 1e-12
%                     disp('eigs issue');
%                     [Vhat e] = eig(Lhat);
%                     Vhat = Vhat(:,n-k+1:end);
%                 end
                
                
                [Vloc piv] = lrcol_rand(Vhat,k,1,5);
                [throw, set] = max(abs(Vloc'));
%                 sucess = 0; 
                if verify_true(set,truth)
%                     sucess = 1; 
                    success(a_count,1) = success(a_count,1)+1;
                end
                
       
                % Now we have a similar matrix A, next build the graph
                % Laplacian matrix based on L = D - A; 
%                 vol = sum(d);
                D = diag(d);
                L = D-A; 
                [V,E] = eig(L);
%                 % Now do the sanity checks for L to see if it satisfies
%                 % othe properties of graph Laplacian.
%                 
%                 % Normalize the graph Laplacian with the random walk
%                 % format.
%                 Lrw = D^(-1)*L;
%                 [Vrw,Erw] = eig(Lrw);

                % Divide into k groups
                
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
                        else
                            Hcpqr(i,j) = 0;
                        end 
                    end 
                end
                scheck2= Hcpqr' * Hcpqr - eye(k,k);
                trace1 = trace(H'*L*H);
                trace2 = trace(Hcpqr'*L*Hcpqr);
                rel_error(a_count) = rel_error(a_count)+(trace1 - trace2)/trace1;
      end 
   end 
end
 
rel_error = rel_error/trials; 
figure
plot(a_range/b,rel_error); 

save RatioCut_k_unequal.mat k m n L A rel_error idx success set truth V E  
fname = [num2str(k) '_block_test_sherlock.mat'];
