%% 1984 LINEAR RECIPROCAL INNERVATION OCULOMOTORMODEL Simulation Code
%% Setting  parameters for moving object
step = 0.001;
finish = 10;
distance = 10;
car_time = transpose(linspace(0,finish,finish/step+1));
%Setting up finishing and beggining positions of moving object
car_value = transpose(linspace(0.7,+0.1,finish/step+1));
%Creating time series for moving object
car = timeseries(car_value,car_time);

fps = 150;%Assuming Human eye fps is 150
%Calculating the degrees for both left and right eye
velocity = 0.6/0.001;
right=round(velocity/fps);%right eye angle
left = atand(tand(right)+0.01);%
left = round(left*10)/10;%left eye angle
angle = [right left];

for eye = 1:2
    tsum = 0;
    tsum1 = 0;
    tsum2 = 0;
    sum1 = 0;
    sum2 = 0;
   
    %% SETTING UP PARAMETERS
Kse = 1.8;
K = 0.86;
B = 0.018;
J = 4.3 * 10e-5;
AG_AC= 13e-3;

%Radian to Degree Parameter
r = pi/180;

step1 = 16;
step2 = 16;
pulse1 = 16;

%% Angle Separator
% This section seperate the coming angle into the 30,20,10,5,1,0.1 which
% are required degrees for  1984 LINEAR RECIPROCAL INNERVATION OCULOMOTORMODEL

%Constants of 30,20,10,5,1,0.1 values
c1 =0;
c2 =0;
c3 =0;
c4 =0;
c5 =0;
c6 =0;
c7 =0;

sum1 = [0 1 2 3 4];
Theta = angle(eye);
Q = Theta;

while(Theta-(c1*30+c2*20+c3*10+c4*5+c5*1+c6*0.5+c7*0.1)~=0)
    if(Q-30>=0)
        while(Q-30>=0)
            Q = Q-30;
            c1 = c1+1;
        end
    end
    if (Q-20>=0)
        while(Q-20>=0)
            Q = Q-20;
            c2 = c2+1;
        end
    end
    if (Q-10>=0)
        while(Q-10>=0)
            Q = Q-10;
            c3 = c3+1;
        end
    end
    if (Q-5>=0)
        while(Q-5>=0)
            Q = Q-5;
            c4 = c4+1;
        end
    end
    if (Q-1>=0)
        while(Q-1>=0)
            Q = Q-1;
            c5 = c5+1;
        end
    end      
    if (Q-0.5>=0)
        while(Q-0.5>=0)
            Q = Q-0.5;
            c6 = c6+1;
        end
    end
    if (Q>=0)
        while(Q-0.1>=-0.0001)
            Q = Q-0.1;
            c7 = c7+1;
        end
    end
end 
      
a = [c1*30 c2*20];

if(c3 > 0)
    c3_a(1) = 10;  
for i = 1:c3   
    c3_a(i) = 10;    
end
a = horzcat(a,c3_a);
end
if(c4 > 0)
    c4_a(1) = 5;  
for i = 1:c4
   
    c4_a(i) = 5;
end
a = horzcat(a,c4_a);
end
if(c5 > 0)
    c5_a(1) = 1;  
for i = 1:c5
    c5_a(i) = 1;
end
a = horzcat(a,c5_a);
end
if(c6 > 0)
    c6_a(1) = 0.5;  
for i = 1:c6
    c6_a(i) = 0.5;
end
a = horzcat(a,c6_a);
end
if(c7 > 0)
    c7_a(1) = 0.1;  
for i = 1:c7
    c7_a(i) = 0.1;
end
a = horzcat(a,c7_a);%a holds the angle values with zeros
end

z = 1;

%Getting rid of the zeros of a
for i = 1:length(a)
    if(a(i)~= 0)
        b(z) = a(i);%b holds the angle values without zeros
        z = z +1;
    end
end
constant = (velocity/fps)/length(b);

%% Simulation
% DESCRIPTIVE TEXT

x = 0;
for i=1:size(transpose(b))
dTheta = b(i);
disp(i);
zero=10e-3;
switch dTheta
case 0.1, PH=17.6; PW=10e-3;
case 0.5, PH=20; PW=10e-3;
case 1, PH=22; PW=11e-3;
case 5, PH=53; PW=15e-3;
case 10, PH=87; PW=20e-3;
case 20, PH=124; PW=31e-3;
case 30, PH=155; PW=40e-3;
end
%Pulse and Step Values
%In this part previous values for pulses and steps are recording

%AG
N_AG_Pulse = PH;
N_AG_Step = step1 + 0.8*(dTheta);
step1 = N_AG_Step;

%ANT
N_AN_Pulse= 0.5 + pulse1*exp((-dTheta)/2.5);
pulse1 = N_AN_Pulse;
N_AN_Step = step2-0.06*(dTheta);
step2 = N_AN_Step;
TauAG_AC =(13-0.1*(dTheta))*1e-3;

%Setting up time parameters for eyemodel simulation
set_param('eyemodel','StartTime','constant*(i-1)/(velocity/fps)')
model = sim('eyemodel',i*constant/(velocity/fps)); % simulate with command script.
% disp("time");
% disp(i*constant/(velocity/20))

%Assigning the parameter values 
if(eye == 1)
if(i == 1)
sum = model.theta;
tsum = model.t;
tsum1 = tsum(length(tsum));
dTheta1 = dTheta;
else
tsum = vertcat(tsum,tsum1+model.t);
sum = vertcat(sum,dTheta1+model.theta);
dTheta1 = dTheta + dTheta1;
tsum1 = tsum(length(tsum));
end
else
if(i == 1)
sum = model.theta;
tsum = model.t;
tsum2 = tsum(length(tsum));
dTheta1 = dTheta;
else
tsum = vertcat(tsum,tsum2+model.t);
sum = vertcat(sum,dTheta1+model.theta);
dTheta1 = dTheta + dTheta1;
tsum2 = tsum(length(tsum));
    end
end
if(eye == 1)
%Right eye data
sum1 = sum;
time_sum1 = tsum;
answer1 = timeseries(sum1,time_sum1);

else
%Left eye data
sum2 = sum;
time_sum2 = 10*tsum/3.5;
answer2 = timeseries(sum2,time_sum2);

end
end
end
time_sum1 = time_sum2;
%% 3D Simulation for eye model
set_param('simulation','StartTime','0');
model = sim('simulation',10);




