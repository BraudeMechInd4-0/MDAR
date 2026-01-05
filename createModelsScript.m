% Test using the chebysev polynomials as model
clc;clear;
%% general parameters
% This section is required to compute r and v from the data downloaded from
% NORAD. If you got the r0 and v0, you can just input them manually and
% delete this section
addpath(genpath("ParseGP"))
addpath(genpath("Satellites"))
addpath(genpath("fundamentals-of-astrodynamics"))
if exist('mu')
    path = which('mu');
    rmpath(path(1:end-4));
end

%% we need this section only whene generating lists
fprintf('loading...');
whichconsts = 84;
Cd = 2.2;
mu = 398600.5;
Re = 6378.137;
%    J1?     J2            J3                   J4              J5                J6
J = [0, 0.01082635854, -0.2532435346e-5, -0.1619331205e-5, -0.2277161016e-6 ,0.5396484906e-6];

delta = 16;

%% generate list of objects and parameters:
[m,A,epoch,r,v,satrecList] = generateList("Satellites\Iridium.xml",datetime("2025-04-03"), duration("15:00:00:00"), 60*24*60*60,whichconsts,false);
%load('Satellites/Iridium33.mat');

t_max = 2*365*24*60*60; %two years

options.AbsTol = 1e-7;
options.RelTol = 1e-5;
options.N = 16;
fprintf('done!\n');
%% create all the models
create_model = tic;
wb = waitbar(0,"creating model");
for i = length(m):-1:1 %going backwards so that on the first run the memory allocation will be made
    a = (2/norm(r(i,:)) -norm(v(i,:))^2/mu)^-1;
    T = sqrt(a^3/mu)*2*pi;
    options.Sec=T/delta;
    [~,~,~,model(i)]=odeMPCI(@(t,x)orbit_eq_J6_drag(t,x,mu,Cd,A(i),m(i),Re,J),[0 t_max],[r(i,:) v(i,:)],options);
    waitbar((length(m)-i+1)/length(m),wb,"creating model");
end
create_model_time = toc(create_model);
waitbar(1,wb,"finished");
close(wb);
fprintf("created model, it took %.2f s\n",create_model_time);
%%
save('Satellites/Iridium33Model.mat', 'model','-v7.3')
%% now compare propagating the satellites to 20 radom timepoints somewhere
%between 0 s and 2 years.

wb = waitbar(0,"comparing runtime");

times = sort(rand(20,1)*t_max);

runtimes_ode = zeros(length(m),length(times));
runtimes_models = runtimes_ode;
posDifference = runtimes_ode;
velDifference = runtimes_ode;

for i = 1:length(m)
    for j = 1:length(times)
        runtime_ = tic;
        a = (2/norm(r(i,:)) -norm(v(i,:))^2/mu)^-1;
        T = sqrt(a^3/mu)*2*pi;
        options.Sec=T/delta;
        [~,xode] = odeMPCI(@(t,x)orbit_eq_J6_drag(t,x,mu,Cd,A(i),m(i),Re,J),[0 times(j)],[r(i,:) v(i,:)],options);
        xode = xode(end,:);

        runtimes_ode(i,j)=toc(runtime_);

        runtime_ = tic;
        xmodel= propCheb(times(j),model(i));
        runtimes_models(i,j)=toc(runtime_);

        posDifference(i,j) = norm(xode(1:3)-xmodel(1:3));
        velDifference(i,j) = norm(xode(4:6)-xmodel(4:6));
        waitbar(((i-1)*length(times)+j)/(length(m)*length(times)),wb,"comparing runtime");
    end
end
waitbar(1,wb,"finished");
close(wb);

%%Now generate plots and data about the difference in time and values
Average_runtimes_ode = sum(runtimes_ode)/length(m);
Average_runtimes_mod = sum(runtimes_models)/length(m);
yyaxis left
plot(times,Average_runtimes_ode)
title("average runtimes")
xlabel("propagation time");
ylabel("runtime ode");
yyaxis right
plot(times,Average_runtimes_model)
ylabel("runtime model");

average_posDiff = sum(posDifference)/length(m);
average_velDiff = sum(velDifference)/length(m);

figure
yyaxis left
plot(times,average_posDiff)
ylabel("position diff")
yyaxis right
plot(times,average_velDiff);
ylabel("velocity diff")
xlabel("propagation time")
title("Errors")