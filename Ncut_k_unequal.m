% In this test, we verify the CPQR by approximating NCut for arbitraty
% k with equal sized-clusters 
clear
close all
clc
rng(1)
k = 4;
m = [120 140 170 210]; 
n = sum(m);
noise = 0;
b= 5;
Na = 20; 
success = zeros(Na,1);
a_range = linspace(5,100,Na);
rel_error = zeros(Na,1);
trials = 9;
a_count = 0;
for a = a_range
        a_count = a_count+1;

    for t = 1:trials
        t
            p = a*log(m(1))/m(1);
            q = b*log(m(1))/m(1);
           
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
                
%                 Ds = diag(sqrt(1./d));
%                 Lhat = Ds*A*Ds;
                  [Vhat e] = eigs(A,k);
%                 [Vhat e] = eigs(Lhat,k);
                
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
                    success(a_count) = success(a_count)+1;
                end
                 
                % Now we have a similar matrix A, next build the graph
                % Laplacian matrix based on L = D - A; 
%                 vol = sum(d);
                D = diag(d);
                L = D-A; 
                [V,E] = eig(L);

                % Now do the sanity checks for L to see if it satisfies
                % othe properties of graph Laplacian.
                
              
                % Normalize the graph Laplacian with the random walk
                % format.
                Lrw = D^(-1)*L;
                [Vrw,Erw] = eig(Lrw);
                H = Vrw(:,1:k); 
                H = D^(-1/2)*H; 
                T = D^(1/2)*H; 
                sanitiyCheck = T'*T - eye(k); 
                obj = trace(T'*D^(-1/2)*L*D^(-1/2)*T); 
                
                % our algo
                vol = [sum(d(1:120));sum(d(121:260));sum(d(261:430));sum(d(431:640))];
            
                H_our = zeros(n,k);
                for i= 1:n
                    for j = 1:k
                        if(set(i)==j) 
                                H_our(i,j) = 1/sqrt(vol(j));
                        end 
                    end 
                end 
                T_our = D^(1/2)*H_our; 
                sanitiyCheck_our = T_our'*T_our - eye(k); 
                obj_our = trace(T_our'*D^(-1/2)*L*D^(-1/2)*T_our); 
                rel_error(a_count) = rel_error(a_count)+(obj - obj_our)/obj;
            end
    end
end
rel_error = rel_error/trials; 
figure
plot(a_range/b,rel_error); 
save NCut_k_unequal4.mat k m n L A rel_error idx success set truth V E Vrw Erw vol 
fname = [num2str(k) '_block_test_sherlock.mat'];
