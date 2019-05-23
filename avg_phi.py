from math import *
#this function takes in any floating point value of phi and spits out the average value of the phi sector in binary 
phi0 = -pi;
def avg_phi(phi): #phi is the float representation of the given phi value of each track
    for i in range(9):
        if (phi >= phi0+(2*i*pi)/9) and (phi < phi0+(2*(i+1)*pi)/9):
            avg_phi = phi0 + (i*pi)/9 + ((i+1)*pi)/9;
            
    avg_phi=round((avg_phi/pi)*(pow(2,16)-1));
    avg_phi=bin(int(str(avg_phi)));
    avg_phi=avg_phi.split('b')[1];
    
    #this is to make sure the length of avg_phi is 17 bits
    S=len(avg_phi);
    diff=17-S;
    
    if diff>=1:
        zero='';
        for i in range(diff-1):
            zero=zero+'0';
    else: 
        zero='0';
        
    #this part is to assign the sign of avg_phi. "0" for positive or zero and "1" for negative
    if phi<0:
        avg_phi='1'+zero+avg_phi;
    elif phi==0:
        avg_phi='0'+zero+avg_phi;
    else:
        avg_phi='0'+zero+avg_phi;
    return avg_phi;

def num_lines(name):
    length=len(open(name).readlines());
    return length;

for i in range(27):
    data=open('phi'+str(i)+'.dat','r+');
    for j in range(num_lines('phi'+str(i)+'.dat')):
        x=data.readline();
        x=int(x,16);
        x=bin(x);
        x=x.split('b')[1];
        
        S=len(x);
        diff=100-S;
        if diff>=1:
            zero='';
            for i in range(diff):
                zero=zero+'0';
        else: 
            zero='0';
        x=zero+x;
        if x[17]==0:
            phi_bin=x[18:35];
            sign=1;
        else:
            phi_bin=x[18:35];
            sign=-1;
        phi=sign*pi*(int(phi_bin,2)/(pow(2,16)-1));
        y=avg_phi(phi);
        x=x.replace(x[17:34],y);
        x=hex(int(x,2));
        data.write(x);
    data.close();