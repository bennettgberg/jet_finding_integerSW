#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <cmath>
#include "hexbin.h"

using namespace std;
const int nphibits = 12;
const float phi_bin_width = (2*M_PI)/27;
float int_to_float(int int_num){
	float floatval = (int_num/(pow(2,nphibits-1)-1))*M_PI;
	return floatval;
}

int main(){
	const int nphibins = 27; //number of phi bins
	const int nphisectors = 9; //number of phi sectors
	const int tpe = 224; //track per event
	ifstream in_tracks[2][nphisectors];
	int ntrks[nphibins];
	for (int i = 0; i < nphisectors; ++i) {
		string fname1, fname2; //"1" is for positive eta files, "2" is for negative eta files.
		ifstream file1, file2;
		string s;
		stringstream out;
		out << i;
		s = out.str();
		fname1 = "phi" + s + "_n.dat";
		fname2 = "phi" + s + "_p.dat";
		in_tracks[0][i].open(fname1.c_str());
		in_tracks[1][i].open(fname2.c_str());
	}
	for (int j = 0; j < nphibins; ++j){
		ntrks[j] = 0;
	}
//initialize new input files
	ofstream in_tracks_final[nphibins];
	for (int i = 0; i < nphibins; ++i){
		string fname;
		string s;
		stringstream out;
		out << i;
		s = out.str();
		fname = "phi" + s + ".dat";
		in_tracks_final[i].open(fname.c_str());
	}

	for (int i = 0; i < nphisectors; ++i){
		string line1, line2;
            	//fill each phi file, from 0 to 26, one event at a time
		while (getline(in_tracks[0][i],line1) && getline(in_tracks[1][i],line2)){
			string phi_bin1 = hex_to_bin(line1.substr(2,24),96).substr(15,12);
			string phi_bin2 = hex_to_bin(line2.substr(2,24),96).substr(15,12); 
			int int1 = bin_to_int(phi_bin1);
			int int2 = bin_to_int(phi_bin2);
			float phi_float1 = int_to_float(int1);
			float phi_float2 = int_to_float(int2); 
			if (line1 != "0x000000000000000000000000" && line2 != "0x000000000000000000000000"){ 
				if (phi_float1 < -M_PI/27){
					in_tracks_final[3*i] << line1 << endl;
					ntrks[3*i] += 1;
					if (phi_float2 < -M_PI/27){
						in_tracks_final[3*i] << line2 << endl;	
						ntrks[3*i] += 1;
					}
					else if (phi_float2 >= -phi_bin_width/2 && phi_float2 < phi_bin_width/2){	
						in_tracks_final[3*i+1] << line2 << endl;
						ntrks[3*i+1] += 1;
					}
					else if (phi_float2 >= phi_bin_width/2){
						in_tracks_final[3*i+2] << line2 << endl;
						ntrks[3*i+2] += 1;
					}
				}
				else if (phi_float1 >= -phi_bin_width/2 && phi_float1 < phi_bin_width/2){
					in_tracks_final[3*i+1] << line1 << endl;
					ntrks[3*i+1] += 1;
					if (phi_float2 < -M_PI/27){
						in_tracks_final[3*i] << line2 << endl;	
						ntrks[3*i] += 1;
					}
					else if (phi_float2 >= -phi_bin_width/2 && phi_float2 < phi_bin_width/2){	
						in_tracks_final[3*i+1] << line2 << endl;
						ntrks[3*i+1] += 1;
					}
					else if (phi_float2 >= phi_bin_width/2){
						in_tracks_final[3*i+2] << line2 << endl;
						ntrks[3*i+2] += 1;
					}
				}
				else if (phi_float1 >= phi_bin_width/2){
					in_tracks_final[3*i+2] << line1 << endl;
					ntrks[3*i+2] += 1;
					if (phi_float2 < -M_PI/27){
						in_tracks_final[3*i] << line2 << endl;	
						ntrks[3*i] += 1;
					}
					else if (phi_float2 >= -phi_bin_width/2 && phi_float2 < phi_bin_width/2){	
						in_tracks_final[3*i+1] << line2 << endl;
						ntrks[3*i+1] += 1;
					}
					else if (phi_float2 >= phi_bin_width/2){
						in_tracks_final[3*i+2] << line2 << endl;
						ntrks[3*i+2] += 1;
					}
				}
			}
			else if (line1 != "0x000000000000000000000000" && line2 == "0x000000000000000000000000"){
				if (phi_float1 < -M_PI/27){
					in_tracks_final[3*i] << line1 << endl;
					ntrks[3*i] += 1;
				}
				else if (phi_float1 >= -phi_bin_width/2 && phi_float1 < phi_bin_width/2){
					in_tracks_final[3*i+1] << line1 << endl;
					ntrks[3*i+1] += 1;	
				}
				else if (phi_float1 >= phi_bin_width/2){
					in_tracks_final[3*i+2] << line1 << endl;
					ntrks[3*i+2] += 1;	
				}
			}
			else if (line1 == "0x000000000000000000000000" && line2 != "0x000000000000000000000000"){
				if (phi_float2 < -M_PI/27){
					in_tracks_final[3*i] << line2 << endl;	
					ntrks[3*i] += 1;
				}
				else if (phi_float2 >= -phi_bin_width/2 && phi_float2 < phi_bin_width/2){	
					in_tracks_final[3*i+1] << line2 << endl;
					ntrks[3*i+1] += 1;
				}
				else if (phi_float2 >= phi_bin_width/2){
					in_tracks_final[3*i+2] << line2 << endl;
					ntrks[3*i+2] += 1;
				}
			}	
		}
		if (ntrks[3*i] < tpe){
			for (int k = 0; k < (tpe - ntrks[3*i]); ++k){
				in_tracks_final[3*i] << "0x000000000000000000000000" << endl;
			}
			
		}
		if (ntrks[3*i+1] < tpe){
			for (int k = 0; k < (tpe - ntrks[3*i+1]); ++k){
				in_tracks_final[3*i+1] << "0x000000000000000000000000" << endl;
			}
			
		}
		if (ntrks[3*i+2] < tpe){
			for (int k = 0; k < (tpe - ntrks[3*i+2]); ++k){
				in_tracks_final[3*i+2] << "0x000000000000000000000000" << endl;
			}
			
		}

		in_tracks_final[3*i].close();
		in_tracks_final[3*i+1].close();
		in_tracks_final[3*i+2].close();
		in_tracks[0][i].close();
		in_tracks[1][i].close();
	}

	return 0;
}
