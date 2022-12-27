clear
close all

param_path = 'C:\Users\adlogan\Dropbox\Frequency Conversion\Lumerical\Projects\Genetic Gratings\';
save_path = 'C:\Users\adlogan\Dropbox\Frequency Conversion\Lumerical\Projects\Genetic Gratings\';
save_file = '2019_06_18_SHG_250';

save_id = fopen([save_path save_file '.txt'],'a+');

obj_fn = @gratingObjectiveGeneral;

max_steps = 200;
num_expand = 1e3;       
num_trials = 12;             % ~99.9% chance of seeing one of best 0.001% neighbors

% 'Scatterer' refers to any feature that can scatter light upwards...

param_file = 'S_param_SHG0_h250_str.mat';   
%Should contain d, end_type, lambda, n_eff, r, s, t for the wavelength set
%used for the objective function

load([param_path param_file],"gr_spec");

s_types = length(gr_spec{1}.d);
end_type = gr_spec{1}.end;

s_range = 3;
% s_range = -3:3;
N_min = 8;             %Indicates the minimum number of standard scatterers
N_max = 12;             %Indicates the maximum number of standard scatterers
L_min = 5e-2;           %Minimum distance between two scatterers (um)
L_max = 0.8;            %Maximum initial distance betwen two scatterers
L_step = 0.005;         %'Resolution' of lengths (um)
L_range = 40;


num_runs = 20;

 sigma = [1.0;
 %         1.5;
          2.0]./2.35;           %Parameters for the objective function
%sigma = 1.0/2.35;           %Parameters for the objective function
     
     
% 
 base_center =  [2.4;
 %                3;
                 2.4];
%base_center =  1.5;
var_center = 0.6;

fprintf(save_id,'%8s %8s %8s %8s\n','Lambda','N_eff','Center','Sigma');
for il=1:length(gr_spec)
    fprintf(save_id,'%8.3f %8.3f %8.3f %8.4f\n',[gr_spec{il}.lambda,gr_spec{il}.neff base_center(il) sigma(il)]);
end
fprintf(save_id,'%s\n','Scatterer Widths');
for ip=1:s_types
    fprintf(save_id,'%6.3f ', gr_spec{1}.d(ip));
end
if end_type
    fprintf(save_id,'%6s\n', 'End');
else
    fprintf(save_id,'%s\n', '');
end
fprintf(save_id,'%s\n', '');
fprintf(save_id,'%s\n', '');


% pool = parpool;
L_final = cell(1,num_runs);
S_final = cell(1,num_runs);
center_final = cell(1,num_runs);
N_final = cell(1,num_runs);
obj_final = cell(1,num_runs);

parfor kk = 1:num_runs
    
%     kk
    fprintf('%s %i\n','Starting run',kk);
    center_rand = rand(length(gr_spec),1);
    center =  base_center + (center_rand-0.5)*var_center;
    
    N = randi([N_min N_max])+end_type;  %Number of scatterers
    if(end_type)
        S =  [randi([1 s_types],1,N-1) s_types+1];
    else
        S =  randi([1 s_types],1,N);
    end
    
    L = round(((rand(1,N-1).*(L_max-L_min))+L_min)/L_step)*L_step;
    
    halt = 0;
    step = 1;
    obj_opt = obj_fn(gr_spec,S,L,center,sigma)
    S_opt = S;
    L_opt = L;
    
    L_test = cell(1,num_expand);
    S_test = cell(1,num_expand);
    obj_test = cell(1,num_expand);
    
    
    while(and(step < max_steps,~halt))
        fprintf('%s %i %s %i %s %e\n','Run',kk,'starting step',step,'with objective',obj_opt(1));
        step = step + 1;
        trial = 0;
        update = 0;
        
        while (and(~update,trial<num_trials))
            trial = trial+1;
            for ix=1:num_expand
                L_test{ix} = L_opt;
                S_test{ix} = S_opt;
                for iy = 1:N-end_type
                    L_test{ix}(iy) = L_opt(iy) + L_step*randi([max(-L_range,ceil((L_min-L_opt(iy))/L_step)),min(L_range,floor((L_max-L_opt(iy))/L_step))]);
                    S_test{ix}(iy) = randi([max(1,S_opt(iy)-s_range),min(s_types,S_opt(iy)+s_range)]);
                end
                obj_test{ix} = obj_fn(gr_spec, S_test{ix},L_test{ix},center,sigma);
            end
            
            for iz = 1:num_expand
                if(obj_test{iz}(1) > obj_opt(1))
                    obj_opt = obj_test{iz};
                    L_opt = L_test{iz};
                    S_opt = S_test{iz};
                    update = 1;
                end
            end
        end
        if update
%             obj_opt
        else
%             'Halting'
%             step
            halt = 1;
        end
    end
    
    if(~halt)
        'Reached max step number'
    end
    
%     S
%     L*1000
    
    L_final{kk} = L;
    S_final{kk} = S;
    center_final{kk} = center;
    N_final{kk} = N;
    obj_final{kk} = obj_opt;
    
    fprintf('%s %i %s %e\n','Run',kk,'finished, objective:',obj_opt(1));
end

for kk = 1:num_runs
    fprintf(save_id,'Run %3i\n',kk);
    fprintf(save_id,'Scatterers: %3i\n', N_final{kk}-end_type);
    fprintf(save_id,'Centers: %8.3f %8.3f %8.3f\n', center_final{kk});
    fprintf(save_id,'Objective: ');
    for io = 1:length(obj_final{kk})
        fprintf(save_id,'%.4e ', obj_final{kk}(io));
    end
    fprintf(save_id,'\n');
    
    for is = 1:length(S_final{kk})
        fprintf(save_id,'%6i', S_final{kk}(is));
    end
    fprintf(save_id,'%s\n','');
    for il = 1:length(L_final{kk})
        fprintf(save_id,'%6i', round(L_final{kk}(il)*1000));
    end
    fprintf(save_id,'%s\n','');
    fprintf(save_id,'%s\n','');
end

save([save_path save_file '.mat'],'gr_spec','L_final','S_final','obj_final');
% delete(pool);

fclose(save_id);
