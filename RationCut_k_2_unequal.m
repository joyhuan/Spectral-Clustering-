% In this test, we verify the CPQR by approximating RatioCut for two  
% unequal-sized clusters
clear
close all
clc
rng(1)
k = 2;
% m =  [70 80 90 100 110 120 130];
m = [90 110]; 
n = sum(m);
noise = 0;

a = 60;
b= 15;
trials = 10;
for t = 1:trials
    t
%     a_count = 0;
%     for a = a_range
%         a_count = a_count+1;
%         b_count = 0;
%         for b = b_range
%             b_count = b_count+1;
            
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
                sucess = 0; 
                if verify_true(set,truth)
                    sucess = 1; 
%                     success(a_count,b_count) = success(a_count,b_count)+1;
                end
                
                % Now we have a similar matrix A, next build the graph
                % Laplacian matrix based on L = D - A; 
                vol = sum(d);
                D = diag(d);
                L = D-A; 
                [V,E] = eig(L);
                
                % Normalize the graph Laplacian with the random walk
                % format.
                Lrw = D^(-1)*L;
                [Vrw,Erw] = eig(Lrw);

                secondG = find(set-1);
                firstG = find(-set+2); 
                
                % First test in a Graph Cut Point of View, i.e.
                cut = 0; 
                for i=secondG
                    for j = firstG
                        cut = cut + A(i,j);
                    end 
                end 
                
                Rcut = cut/length(firstG) + cut/length(secondG);
                vol1 = 0;
                for i = firstG
                        vol1 = vol1 + D(i,i);
                end 
                vol2 = 0;
                for j = secondG
                        vol2 = vol2 + D(j,j);
                end 
                Ncut = cut/vol1 + cut/vol2; 
               
                % Want to minimize RatioCut/Ncut
                % First do RatioCut

                % By closed form solution, we know the optimal 
                % f is the eigenvector corresponding to the second smallest
                % eigenvalue.
                f = V(:,2);
                scheckOrthog = sum(f.*ones(n,1));
                scheckdiff = norm(f)-sqrt(n);
                groupA = [];
                groupB = [];
                for i = 1:n
                    if f(i)>= 0 
                        groupA = [groupA i];
                    else
                        groupB = [groupB i];
                    end 
                end 
                
                
                for i = 1:n
                    if f(i)>= 0 
                        f(i)=2;
                    else
                        f(i)=1;
                    end 
                end 
                % Now check the difference between this closed
                % form solution and the one computed by our CPQR algo. 
               
                group1diff = length(firstG)-length(groupA);
                group2diff = length(secondG)-length(groupB);
                FINAL = verify_true(f,truth)

                % Next do Ncut (the optimization is Eq(9) in the paper)
                % where the solution is achieved at the second eignevector
                % of Lrw
                frw = Vrw(:,2);
                checkOrthogRW = sum(D*frw.*ones(n,1));
                diffRW= frw'*D*frw-vol;
%                 
%                 % this only makes sense in the two block case, Fiedler
%                 % vector test
%                 if k == 2
%                     F = V(:,2);
%                     F(F>=0) = 1;
%                     F(F<0) = 2;
%                     if verify_true(F',truth)
%                         fiedler(a_count,b_count) = fiedler(a_count,b_count)+1;
%                     end
%                 else
%                     F = kmeans(V,k);
%                     if verify_true(F',truth)
%                         fiedler(a_count,b_count) = fiedler(a_count,b_count)+1;
%                     end
%                 end
%                 
%                 
%                 [V e] = eigs(A,k);
%                 if norm(V'*V - eye(k)) > 1e-12
%                     disp('eigs issue');
%                     [V e] = eig(A);
%                     V = V(:,n-k+1:end);
%                 end
%                 
%                 [Vloc piv] = lrcol_rand(V,k,1,5);
%                 [throw, set] = max(abs(Vloc'));
%                 if verify_true(set,truth)
%                     success2(a_count,b_count) = success2(a_count,b_count)+1;
%                 end
%                 
%                 if max(comp_cut(A,set)) <= max(comp_cut(A,truth))
%                     comp_true(a_count,b_count) = comp_true(a_count,b_count) + 1;
%                 end
%                 if comp_mle(A,set) >= comp_mle(A,truth)
%                     MLE(a_count,b_count) = MLE(a_count,b_count) + 1;
%                 end
%                 
%                 
%             end
%             a_grid(a_count,b_count) = a;
%             b_grid(a_count,b_count) = b;
%             p_grid(a_count,b_count) = p;
%             q_grid(a_count,b_count) = q;
%             if k == 2
%                 theory(a_count,b_count) = (a+b)/2-sqrt(a*b);
%             else
%                 theory(a_count,b_count) = sqrt(a)-sqrt(b);
%             end
%  
%         end
    end
end
% success = success/trials;
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
save testEasy.mat k m n firstG secondG L Lrw A 
fname = [num2str(k) '_block_test_sherlock.mat'];
