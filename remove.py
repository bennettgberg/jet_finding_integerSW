import os
for i in range(27):
	x = 'phi'+str(i)+'.dat';
	os.remove(x);
for i in range(9):
	x = 'phi'+str(i)+'_p.dat';
	os.remove(x);
for i in range(9):
	x = 'phi'+str(i)+'_n.dat';
	os.remove(x);

