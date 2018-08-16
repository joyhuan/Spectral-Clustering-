% In this test, we verify the CPQR by approximating NCut for arbitraty
% k with equal sized-clusters 
clear
close all
clc
rng(1)
k = 2;
m = [100 130]; 
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
                sucess = 0; 
                if verify_true(set,truth)
                    sucess = 1;
%                     success(a_count) = success(a_count)+1;
                end
                 
                % Now we have a similar matrix A, next build the graph
                % Laplacian matrix based on L = D - A; 
                vol1 = sum(d);
                D = diag(d);
                L = D-A; 
                [V,E] = eig(L);

                % Now do the sanity checks for L to see if it satisfies
                % othe properties of graph Laplacian.
                
              
                % Normalize the graph Laplacian with the random walk
                % format.
                Lsym = D^(-1/2)*L*D^(-1/2); 
                [Vsym,Esym] = eig(Lsym);
%                 Lrw = D^(-1)*L;
%                 [Vrw,Erw] = eig(Lrw);
%                 H = Vrw(:,1:k); 
%                 H = D^(-1/2)*H; 
%                 T = D^(1/2)*H; 
                g = Vsym(:,2); 
                sanityCheck1 = dot(g,D^(1/2)*ones(n,1)); 
%                 sanityCheck2 = (norm(g))^2 - vol1; 
%                 T = Vsym(:,1:k);
%                 sanitiyCheck = T'*T - eye(k); 
%                 obj = trace(T'*D^(-1/2)*L*D^(-1/2)*T); 
                obj = g'*D^(-1/2)*L*D^(-1/2)*g*vol1;
                % our algo
                vol = [sum(d(1:100));sum(d(101:200))];
                f = zeros(n,1);
                for i= 1:n
                    
                        if(set(i)==1) 
                            f(i) = sqrt(vol(2)/vol(1));
                        else
                            f(i) = -sqrt(vol(1)/vol(2));
                        end 
                end 
                sanityCheck1_our = dot(D*f,ones(n,1)); 
                sanityCheck2_our = f'*D*f- vol1; 
                
                obj_our = f'*L*f; 
%                 Final = (obj - obj_our)/obj;
                rel_error(a_count) = rel_error(a_count)+(obj - obj_our)/obj;
            end
    end
end
rel_error = rel_error/trials; 
figure
plot(a_range/b,rel_error); 
save NCut_sym_2_unequal.mat k m n L A rel_error idx success set truth V E 
fname = [num2str(k) '_block_test_sherlock.mat'];
