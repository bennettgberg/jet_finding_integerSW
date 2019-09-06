import math
import os

#This script will convert a coe file (in binary format) to a usable hexadecimal file for input to C++ emulation.
nfibers = 12 #how many fibers have data in the coe file
nphi = 27 #how many total phi sectors there are (zero out unused ones)
wordlength = 96 #the number of bits per track
tpe = 112  #tracks per event (for each fiber)
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
coe.close()

print("Data written successfully");
