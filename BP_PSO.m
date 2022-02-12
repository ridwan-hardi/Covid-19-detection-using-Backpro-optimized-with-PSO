clc; clear all; close all

%% inisiasi data
%data = xlsread ('E:\COVID 19\hasil program\dataset GLCM.xlsx');
[data,txt,raw] = xlsread('E:\COVID 19\hasil program\dataset GLCM.xlsx');

% variabel = txt(:,2:5);            %sudut 0 derajat
% data1 = data(:,2:5);
% target_latih = data(:,1)';

% variabel = txt(:,13:17);           %sudut 45 derajat
% data1 = data(:,14:17);
% target_latih = data(:,13)';

% variabel = txt(:,29:32);             %sudut 90 derajat
% data1 = data(:,29:32);
% target_latih = data(:,28)';

variabel = txt(:,44:47);             %sudut 180 derajat
data1 = data(:,44:47);
target_latih = data(:,43)';


percobaan = [4];
data_latih = data1(:,percobaan)';    % 1: contrast 2: correlation 3: Energy 4: Homogenity

[~,N] = size(data_latih);

%% pHasarameter PSO

nvar        = 1;                    %banyaknya variabel
varsize     = [1 nvar];             %matriks besar variabel
lb          = [5 0.1];              %lower
ub          = [40 0.9];             %upper
maxit       = 1;                   %maximum iterasi
nPop        = 1;                   %banyaknya populasi
w           = 1;                    %inersia
wdamp       = 0.99;                 %damping ratio of inertia coeff

namafile = ['Hasil PSO_Homogenity 135 derajat_',num2str(nPop)];

kappa = 1;
phi1 = 2.05;
phi2 = 2.05;
phi = phi1 + phi2;
chi = 2*kappa/abs(2-phi-sqrt(phi^2-4*phi));

c1          = chi*phi1;                    %koefisien 1
c2          = chi*phi2;                    %koefisien 2
MaxVelocity = 0.2.*(ub-lb);
MinVelocity = -MaxVelocity;
%% inisiasi
empty_particle.position         = [];
empty_particle.velocity         = [];
empty_particle.cost             = [];
empty_particle.best.position    = [];
empty_particle.best.cost        = [];

particle = repmat (empty_particle, nPop, 1);

%inisiasi Gbest
GlobalBest.cost = inf;
%% membentuk posisi dan kecepatan awal serta menentukan Pbest dan Gbest awal
for i=1:nPop
    %generate random solution
    particle(i).position        = unifrnd (lb,ub,varsize); %uniform random fungtion
    particle(i).position(1)     = round(particle(i).position(1));
   
    % inisiasi velocity
    particle(i).velocity        = zeros(varsize);
    %% code BP
    net_train = bp_co (data_latih, target_latih, particle(i).position(2), particle(i).position(1))
    
    %% menentukan hasil cost
    %particle(i).cost            = fgsiobjek(particle(i).position)
    particle(i).cost            = net_train.error;
    % update
    particle(i).best.position   = particle(i).position;
    particle(i).best.cost       = particle(i).cost;
    
    %update Gbest
    if particle(i).best.cost < GlobalBest.cost
        GlobalBest = particle(i).best;
    end
end
Bestcost= zeros(maxit,1);
costgoal= 0.01;
    
%%main loop 
for it=1:maxit 
    if GlobalBest.cost <= costgoal
        BestCost(it) = GlobalBest.cost;
        break
    end
        for i=1:nPop 
            particle(i).velocity = w*particle(i).velocity ...             
                +c1*rand(varsize).*(particle(i).best.position -  particle(i).position)... 
                +c2*rand(varsize).*(GlobalBest.position -  particle(i).position);

            particle(i).velocity = max(particle(i).velocity, MinVelocity);
            particle(i).velocity = min(particle(i).velocity, MaxVelocity);

            particle(i).position = particle(i).position + particle(i).velocity; 
            particle(i).position(1) = round(particle(i).position(1));
            
            particle(i).position = max(particle(i).position, lb);
            particle(i).position = min(particle(i).position, ub);

            %% code bp
            net = bp_co (data_latih, target_latih, particle(i).position(2), particle(i).position(1));
            

            %% menentukan hasil cost
            particle(i).cost = net.error;


            if particle(i).cost < particle(i).best.cost  
                particle(i).best.position   = particle(i).position; 
                particle(i).best.cost       = particle(i).cost;

               if particle(i).best.cost < GlobalBest.cost  
                   GlobalBest = particle(i).best;
                    
               end
            end
        end 
        BestCost(it) = GlobalBest.cost; 
        

    %% display hasil cost
        disp(['iteration: ',num2str(it) '  Best Cost: ', num2str(BestCost(it))]); 
        % damping inertia coeff
        w=w*wdamp;
        
end    


    %% plot HASIL
BestCost
GlobalBest
figure ()
% plot(BestCost, 'LineWidth', 2);
semilogy(BestCost, 'LineWidth', 2);
title('Hasil PSO')
grid on
xlabel ('iteration')
ylabel ('Best Cost')


%%
% data2 = data(1:50,8:11);     % 0 derajat
% target_uji = data(1:50,7)';

% data2 = data(1:50,21:24);    %45 derajat
% target_uji = data(1:50,20)';

% data2 = data(1:50,36:39);    %90 derajat
% target_uji = data(1:50,35)';

data2 = data(1:50,52:55);      %180 derajat
target_uji = data(1:50,51)';

data_uji = data2(1:50,percobaan)'; % 1: contrast 2: correlation 3: Energy 4: Homogenity

net_keluaran_mainloop = net.Net;

% target_uji= kesimpulan_glcm([1:50,181:220],1)';
% data_uji = glcmsudut_135([1:50,181:220],1:3)';

[~,M] = size(data_uji);
hasil_uji = round(sim(net_keluaran_mainloop,data_uji));
[u,v] = find(hasil_uji==target_uji);

akurasi = (sum(u)/M)*100

[a,aa] = find(hasil_uji == 0);
negatif_benar = 0;
for i = 1:length(aa)
    if aa(i) > 25
        negatif_benar = negatif_benar+1;
    end
end
negatif_palsu   = length(a)-negatif_benar;

[b,bb] = find(hasil_uji == 1);
positif_benar = 0;
for i = 1:length(bb)
    if bb(i) <= 25
        positif_benar = positif_benar+1;
    end
end
positif_palsu = length(bb) - positif_benar;
spesifitas      = (negatif_benar/(positif_palsu + negatif_benar))*100
sensitivitas    = (positif_benar/(positif_benar + negatif_palsu))*100
%% menyimpan gambar dan workspace
savefig (['E:\COVID 19\hasil program\',namafile,'.fig'])
save(['E:\COVID 19\hasil program\',namafile,'.mat'])
