
function [Eopt,b3,z3]=Bisection_method_E_nondim_b_u(Rm)
    if Rm<0
        Eopt=0
        b3=0
        z3=0
    else
        Emin=2/3^(3/2);
        Emax=100000000;
        ERRORmax=10^(-3);
        ERROR=1;
        counter=1;
        while ERROR>ERRORmax
                counter=1+counter;
                b0 = 1;
                zspan = [0 3];
                Emid=(Emin+Emax)/2;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%% MIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Popt=[Emin Rm];
                options = odeset('Events',@event_function);
                [z1,b1,zEmin,bEmin] = ode45(@myode,zspan,b0,options,Popt);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Popt=[Emax Rm];
                options = odeset('Events',@event_function);
                [z2,b2,zEmax,bEmax] = ode45(@myode,zspan,b0,options,Popt);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%% MID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Popt=[Emid Rm];
                options = odeset('Events',@event_function);
                [z3,b3,zEmid,bEmid] = ode45(@myode,zspan,b0,options,Popt);
                if zEmid>1
                    Emin=Emid;
                elseif zEmid<1
                    Emax=Emid;
                end
                if counter>10 && Emid<((2/3^(3/2))+0.00001)
                    ERROR=ERRORmax-0.00001;
                    else
                    ERROR=abs(zEmid-1);
                end 
        end
    Eopt=Emid
    ERROR
    end 
end
function dbdz = myode(z,b,Popt)
Rm=Popt(2);
E=Popt(1);
dbdz=-Rm*(E-b*(1-b^2));
end
function [value,isterminal,direction] = event_function(z,b,Popt)
value = b;  % when value = 0, an event is triggered
isterminal = 1; % terminate after the first event
direction = 0;  % get all the zeros
end