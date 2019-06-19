#include "hexbin.h"
#include "jet_header.h"

/* Software to emulate the hardware 2-layer jet-finding algorithm (integers). *Layers 1 and 2*           
 *
 * 2019
 *
 * Authors: Bennett Greenberg, Felipe Paucar-Velasquez, Yuri Gershtein
 * Rutgers, the State University of New Jersey
 *  Revolutionary for 250 years
 */

//Holds data from tracks, converted from their integer versions.
int nzbins = 6;
int main(int argc, char ** argv){
	nzbins = 6;
	int eventstart = 0;            
	int eventend = 999999;         
	if(argc == 2){
		nzbins = atoi(argv[1]);                                                 
	}
	else if(argc == 3){
		eventstart = atoi(argv[1]);
		eventend = atoi(argv[2]);
	}
	else if(argc == 4){
		nzbins = atoi(argv[1]);
		eventstart = atoi(argv[2]);
		eventend = atoi(argv[3]);
	}
	else if(argc > 4){
		cout << "Unable to understand arguments. Please include only nzbins, or eventstart eventend, or nzbins eventstart eventend." << endl;
		exit(0);
	}
	cout << "Running with " << nzbins << " zbins. " << endl;
	string fname = "phi"; // "/home/mapsa/fpga/emulation/phi";  
     //abs(pT) must be >2, so abs(1/pT) will be <0.5
//	const float pTinvmax = 0.5;
     //number of bits the firmware will use for pT
	const int pTbits = 15;
	const int ntbits = 16; //this is for eta
	const int nzbits = 12;
     //value of one bit of pT
//	const float pTstep = pTmax / (pow(2, pTbits)-1);
     //the three parts of the incoming track
	int pTinverse; //as of now, it is actually just pT
	int t; //really just eta now
	int z0;
     //float values of these
//	double pTinvf;
	double pT;
//	double tf;
	int counter;
	string data_in;
	string bin_data;
	string filename;
//#  This version should be consistent with the interface specified in the technical proposal:
//#	15 bits of pT* (note: NOT 1/pT)
//#	12 bits of phi0
//#	13 bits of d0
//#	12 bits of z0*
//#	16 bits of t* (note: actually just eta)
//#	4 bits of chi2
//#     3 bits of bend-chi2
//#     7 bits of hit mask
//#     3 bits of track quality MVA
//#     6 bits for space for two specialized MVA selections
//#	5 bits of spare
//# For a total of a 96 bit track. For our purposes all will be 0 except the ones with *s.
	//Open input file, read data to tracks array
	//number of tracks in each phi bin is 224 for now.
	ifstream in_tracks[nphibins];
	int ntrks[nphibins];
	for(int i = 0; i < nphibins; ++i){
		string s;
		stringstream out;
		out << i;
		s = out.str();
		filename = fname + s + ".dat";
		in_tracks[i].open(filename.c_str());
		ntrks[i] = 0;
	}
	
	ofstream out_clusts;
	string outname = "int_em_out.txt";
	out_clusts.open(outname.c_str());
	int zint;
	int etaint;
	int pTint;
	int phi;
	string data;
	string bin_z;
	string bin_eta;
	string bin_pT;
	string bin_phi;
	string bin_nt;
	string bin_nx;
	int ntracks;
	int maxntrk = 0;
	int last_maxntrk;
	int nevents = 0;
	int nreal_events = 0;
	for(nevents = 0; nevents <= eventend; ++nevents){
	    ntracks = 0;
	    last_maxntrk = maxntrk;
	    maxntrk = 0;
	    vector<track_data> tracks(0);
         //read from one phi slice at a time
	    for(int pslice = 0; pslice < nphibins; ++pslice){
		getline(in_tracks[pslice], data_in);
		if(data_in == "") {
			    goto data_read;
		}
	//Get all the lines of 0s of the last event (consequence of all phisectors having same number of tracks for each event).
		for(int tk=ntrks[pslice]; tk < last_maxntrk; ++tk){
			getline(in_tracks[pslice], data_in);	
		}
		ntrks[pslice] = 0;
		while(true) {
			//if data is 0, event has ended. 
			if(data_in == "0x000000000000000000000000"){                           
				if(nevents < eventstart){
					ntracks = 0;
					break;   
				}
				break;                                                       
			}//end data is 0; event has ended
			if(data_in == "") { //if end of file reached, we're all done reading in data.
				    goto data_read;
			}
			ntrks[pslice]++;
			bin_data = hex_to_bin(data_in, 96); 
//pTinverse-->pT_actual
			pTinverse = bin_to_int(bin_data.substr(0, 15));
			t = bin_to_int(bin_data.substr(27, 16)); //this is actually eta
			z0 = bin_to_int(bin_data.substr(43, 12));
			track_data trkd;
			if(bin_data.substr(95, 1) == "1"){
				trkd.xbit = true;
			}
			else {
				trkd.xbit = false;
			}
			//Convert integers pTinverse, t, and z0 to their correct
			// floating point values, store in track_data struct
//			if(1.0*pTinverse < pow(2, nptbits-1)){
//				pTinverse += (int)pow(2,nptbits-1);
//			}
//			else {
//				pTinverse -= (int)pow(2, nptbits-1);
//			}
	//		pTinvf = -1.0* /*pTinvmax*/ pTmax + pTinverse * 2 * pTmax / (pow(2,nptbits)-1);
			//taking only the absolute value of pT.
//			pT = fabs(/*1.0 / */ pTinvf);
			pT = pTinverse*1.0;
			if(pT > 511) pT = 511;
			trkd.pT = (int)pT;
		//	tracks[ntracks].pT = (int)round(pT / pTstep);
			//cout << data_in << " pT: " << trkd.pT << endl;
			trkd.eta = t*1.0;
			//if(1.0 * t < pow(2, ntbits-1)){
			//	t += (int)pow(2, ntbits-1);
			//}
			//else {
			//	t -= (int)pow(2, ntbits-1);
			//}
			
			//tf = -maxt + t * 2.0 * maxt / (pow(2, ntbits)-1);
			//trkd.eta = -1.0 * log(sqrt(tf*tf + 1.0) - tf);
			if(trkd.eta > maxeta) {
				trkd.eta = maxeta;
			}
			else if(trkd.eta < -1.0 * maxeta) {
				trkd.eta = -1.0 * maxeta;
			}
			
			if(1.0 * trkd.eta < pow(2, ntbits-1)){
				trkd.eta += (int)pow(2, ntbits-1);
			}
			else {
				trkd.eta -= (int)pow(2, ntbits-1);
			}
			
			if(1.0 * z0 < pow(2, nzbits-1)){
				z0 += (int)pow(2, nzbits-1);
			}
			else {
				z0 -= (int)pow(2, nzbits-1);
			}
			trkd.z = -1.0*maxz + z0 * 2.0 * maxz /( pow(2, nzbits)-1);
			trkd.phi = -1.0*maxphi + pslice * phistep + (phistep / 2);  
			trkd.eta = -1.0*maxeta + trkd.eta*2.0*maxeta/(pow(2,ntbits)-1);
			trkd.bincount = 0;
			++ntracks;                                         
			tracks.push_back(trkd);
			getline(in_tracks[pslice], data_in);
		}
		if(maxntrk < ntrks[pslice]) maxntrk = ntrks[pslice];
	}
			
	//Call L2cluster
	//
	//Then print to output file all clusters from the most energetic zbin 
data_read:
		if(ntracks == 0){ continue;}
		cout << "****EVENT " << nreal_events << " ****" << endl;
		nreal_events++;
		maxzbin mzb = L2_cluster(tracks, nzbins, ntracks); //left off here
		if(mzb.isEmpty == true) { 
			continue;
		}
	for(int kk = 0; kk < nzbins-1; ++kk){	
	        if(kk != mzb.znum) continue;
		vector<etaphibin> clusters = all_zbins[kk].clusters;
		for(int k = 0; k < all_zbins[kk].nclust; ++k){
			//Convert z, eta, and pT to their integer equivalents.
			//Measure from the middle of each bin to avoid errors.
			//zbins will go from 0 to 4.
			counter = 0;
			zint = all_zbins[kk].znum;
			//etabins will go from 0 to 23, not -12 to 11.
			counter = 0;
			while(clusters[k].eta > -1 * maxeta + etastep + counter * etastep) {
				++counter;
			}
			etaint = counter;
			counter = 0;
			while(clusters[k].phi > -1 * maxphi + phistep + counter * phistep) {
				++counter;
			}
			phi = counter;
			//pT is in absolute value form.
			if(clusters[k].pTtot < pow(2, pTbits)){
				pTint = clusters[k].pTtot;
			}
			else {
				pTint = (1 << pTbits) - 1;
			}
			int nt = clusters[k].numtracks;
			int nx = clusters[k].nx_tracks;
			//Convert integers to binary strings.
//			cout << "nclust: " << k << " pT: " << pTint << endl;
			bin_z = int_to_bin(zint, 4);
			bin_eta = int_to_bin(etaint, 5);
			bin_pT = int_to_bin(pTint, pTbits);
			bin_phi = int_to_bin(phi, 5);
			bin_nt = int_to_bin(nt, 5);
			bin_nx = int_to_bin(nx, 4);
			//Concatenate into one binary string, convert to hex, write to output file.
			data = bin_nt + bin_nx + bin_z + bin_eta + bin_phi + bin_pT;
			data = bin_to_hex(data);
			out_clusts << data << endl;		
		}
         } //for each zbin
		if(mzb.ht == 0) cout << "WARNING: HT = 0 (Event " << nevents << ")" << endl;
		out_clusts << "00000000" << endl;
    }
	for(int it = 0; it < nphibins; ++it) {
        	in_tracks[it].close();
	}
	out_clusts.close();
	cout << "Data written to " << outname << endl;
	return 0;
}
