import math
import os
#this function takes in any floating point value of phi and spits out the average value of the phi sector in binary 
phi0 = -math.pi;
def avg_phi(phi): #phi is the float representation of the given phi value of each track
    for i in range(9):
        if (phi >= phi0+(2*i*math.pi)/9) and (phi < phi0+(2*(i+1)*math.pi)/9):
            avg_phi = phi0 + (i*math.pi)/9 + ((i+1)*math.pi)/9;
            
    avg_phi=round((avg_phi/math.pi)*(pow(2,16)-1));
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

#This script will convert a coe file (in binary format) to a usable hexadecimal file for input to C++ emulation.
nfibers = 12 #how many fibers have data in the coe file
nphi = 27 #how many total phi sectors there are (zero out unused ones)
wordlength = 100 #the number of bits per track
tpe = 24  #tracks per event (for each fiber)
nevents = 161 #how many events to make inputs for
coe = open("vcu118_input_patterns.coe", "r")
#read info lines at top of file
for i in range(10):
	coe.readline()

for event in range(nevents):
	for track in range(tpe):
		all_tracks = coe.readline()
		for phi in range(nphi):
			if phi < nfibers:
				data = ""
				for i in range(phi*wordlength, (phi+1)*wordlength, 4):
					hexnum = hex(int(all_tracks[i:i+4], 2))
					#print "hexnum: " + str(hexnum)
					data = data + str(hexnum)[2]
			else:
				zlist = ['0' for i in range(wordlength)]
				#print str(zlist)
				data = "".join(zlist)
			fname = "phi" + str(phi) + ".dat"
			if event == 0 and track == 0:
				let = 'w'
			else:
				let = 'a' #if not first time opening file, append to it (don't overwrite)
			phifile = open(fname, let)
			phifile.write("0x" + data + "\n")
			phifile.close()
	#now write 0s to signify end of event.
	for phi in range(nphi):
		fname = "phi" + str(phi) + ".dat"
		phifile = open(fname, 'a')
		phifile.write("0x" + "".join(['0' for i in range(int(math.ceil(wordlength/4)))]) + "\n")

print("Data written successfully")
coe.close()

def num_lines(name):
    length=len(open(name).readlines());
    return length;

for i in range(nphi):
    data=open('phi'+str(i)+'.dat','r');
    data_w=open('phi'+str(i)+'_mod.dat','w')
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
        phi=sign*math.pi*(int(phi_bin,2)/(pow(2,16)-1));
        y=avg_phi(phi);
        x=x.replace(x[17:34],y);
        x=hex(int(x,2));
        x=x.split('x')[1];
        if len(x)<25:
            zlist=['0' for k in range(25-len(x))];
            zeros=''.join(zlist);
        else:
            zeros='';
        x='0x'+zeros+x;
        data_w.write(x+"\n");
    data.close();
    data_w.close();
for i in range(nphi):
    os.remove('phi'+str(i)+'.dat');  
for i in range(nphi):
    os.rename('phi'+str(i)+'_mod.dat','phi'+str(i)+'.dat');