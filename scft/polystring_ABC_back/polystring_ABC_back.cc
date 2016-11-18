/*
 * =====================================================================================
 *
 *       Filename:  polystring.cc
 *
 *    Description:  for string method
 *
 *        Version:  1.0
 *        Created:  09/30/2016 08:23:57 PM
 *       Revision:  none
 *       Compiler:  icc
 *
 *         Author:  Junqing Song
 *        Company:  FDU
 *
 * =====================================================================================
 */

#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <thread>
#include <mutex>

#include "spline.h" //for cublic spline interpolation, refer to "https://github.com/ttk592/spline".

#include "common.h"
#include "Grid.h"
#include "Config.h"
#include "scft.h"
#include "models.h"
#include "armadillo"
#include <blitz/array.h>

using namespace std;
using namespace arma;
using namespace blitz;

void run_string(string);
void initialize_string();
void parameterize_string();
void redistribute_stringBeads();
void save_data();
void run_scftInString(string);
void thread_scft(Model *, string, Array<double, 4>, Array<double, 4>, Array<double, 4>, int);
Array<double, 3> input_data3(string, string);
Array<double, 4> input_data4(string, string);


int m; //number of beads on string.
int maxT; //max iteration of string.
Array<double, 1> S; //string nodes
Array<double, 1> H; //free energy of each node.
Array<double, 1> F; //whole free energy of one string
Array<double, 4> wa, wb, wc, phia, phib, phic; //field and density data in each string node.
int Nxx;    //grid size
int Nyy;
int Nzz;
bool is_readString; //initialization from an existing string
Model *pmodel;
std::mutex mtx; //need to below protected variable

int main(int argc, char* argv[]){
	// the default configuration file is in the current working directory
    string config_file = "paramString.ini";
    Config cfg(config_file.c_str());
    m = cfg.get_integer("string","num_of_beads");
    maxT = cfg.get_integer("string","string_iteration");
    is_readString = cfg.get_bool("string", "is_readString");
    S.resize(m);
    H.resize(m);
    F.resize(maxT);
    // Specify configuration file through input arguments.
    if(argc > 1){
        config_file = argv[1];
	}
    run_string(config_file);

    return 0;
}

void run_string(string config_file) {
	blitz::Range all = blitz::Range::all();
	Config cfg(config_file.c_str());
	if(cfg.model() == ModelType::AB_C) {
        pmodel = new Model_AB_C(config_file);  
    } //do not construct ETDRK4 frequently in string iterations, because it is time-consuming.
    int st = 0; //index of string iteration.
	Nxx = cfg.Lx();
	Nyy = cfg.Ly();
	Nzz = cfg.Lz();
	wa.resize(m, Nxx, Nyy, Nzz);
	wb.resize(m, Nxx, Nyy, Nzz);
    wc.resize(m, Nxx, Nyy, Nzz);
	phia.resize(m, Nxx, Nyy, Nzz);
	phib.resize(m, Nxx, Nyy, Nzz);
    phic.resize(m, Nxx, Nyy, Nzz);

	initialize_string();
    do {
        cout << "**************************************************" << endl;
        cout << "This is the " << st << "th string iteration!" << endl;
        cout << "**************************************************" << endl;
        run_scftInString(config_file);
        F(st) = mean(H);
        cout << "H = " << H << endl;
        cout << "F = " << F(st) << endl;
        parameterize_string();
        save_data();
        cout << "S = " << S << endl;
        redistribute_stringBeads();
        //save_data();
        st += 1;
    }
    while(abs(F(st)-F(st-1))>1.0e-7 && st<maxT);
    
	
}

Array<double, 3> input_data3(string str1, string str2) {  //read 3D matrix
	string filename;
    string keyname;
    filename = str2;
    keyname = str1;
    Array<double, 3> data(Nxx, Nyy, Nzz, blitz::fortranArray);
    CMatFile mat;
    mat.matInit(filename.c_str(), "r");
    if(!mat.queryStatus()) {   //be necessary for different dimensions
    	mat.matGetArray(keyname, data.data(), data.size()*sizeof(double));
    	mat.matRelease();
    	return data;
	}
}

Array<double, 4> input_data4(string str1, string str2) {  //read 4D matrix
	string filename;
    string keyname;
    filename = str2;
    keyname = str1;
    Array<double, 4> data(m, Nxx, Nyy, Nzz, blitz::fortranArray);
    CMatFile mat;
    mat.matInit(filename.c_str(), "r");
    if(!mat.queryStatus()) {   //be necessary for different dimensions
    	mat.matGetArray(keyname, data.data(), data.size()*sizeof(double));
    	mat.matRelease();
    	return data;
	}
}

void initialize_string() {
	blitz::Range all = blitz::Range::all();
    if(is_readString) {
    	cout << "Initialization from an existing string!" << endl;
    	wa = input_data4("wa", "stringData.mat");
    	wb = input_data4("wb", "stringData.mat");
        wc = input_data4("wc", "stringData.mat");
    	phia = input_data4("phia", "stringData.mat");
    	phib = input_data4("phib", "stringData.mat");
        phic = input_data4("phic", "stringData.mat");
    }
    else {
    	cout << "Inilializing string by two ends!" << endl;
		/**************************initialize the starting bead and the end bead************************/
		wa(0, all, all, all) = input_data3("wA", "defect.mat");
    	wb(0, all, all, all) = input_data3("wB", "defect.mat");
        wc(0, all, all, all) = input_data3("wC", "defect.mat");
    	phia(0, all, all, all) = input_data3("phiA", "defect.mat");
    	phib(0, all, all, all) = input_data3("phiB", "defect.mat");
        phic(0, all, all, all) = input_data3("phiC", "defect.mat");
    	wa(m-1, all, all, all) = input_data3("wA", "vLam.mat");
    	wb(m-1, all, all, all) = input_data3("wB", "vLam.mat");
        wc(m-1, all, all, all) = input_data3("wC", "vLam.mat");
    	phia(m-1, all, all, all) = input_data3("phiA", "vLam.mat");
    	phib(m-1, all, all, all) = input_data3("phiB", "vLam.mat");
        phic(m-1, all, all, all) = input_data3("phiC", "vLam.mat");

    	/*************************initialize string for spinodal transition*********************/
    	for (int i=0; i<m; i++) {
    		wa(i, all, all, all) = wa(0, all, all, all) + 1.0*i/(m-1)* (wa(m-1, all, all, all) - wa(0, all, all, all));
    		wb(i, all, all, all) = wb(0, all, all, all) + 1.0*i/(m-1)* (wb(m-1, all, all, all) - wb(0, all, all, all));
            wc(i, all, all, all) = wc(0, all, all, all) + 1.0*i/(m-1)* (wc(m-1, all, all, all) - wc(0, all, all, all));
        	phia(i, all, all, all) = phia(0, all, all, all) + 1.0*i/(m-1)* (phia(m-1, all, all, all) - phia(0, all, all, all));
        	phib(i, all, all, all) = phib(0, all, all, all) + 1.0*i/(m-1)* (phib(m-1, all, all, all) - phib(0, all, all, all));
            phic(i, all, all, all) = phic(0, all, all, all) + 1.0*i/(m-1)* (phic(m-1, all, all, all) - phic(0, all, all, all));
    	}
    }
}

void save_data() {
	CMatFile mat;
    mat.matInit("result.mat","w");
    mwSize mm = (mwSize)m;
    mwSize Nx = (mwSize)Nxx;
    mwSize Ny = (mwSize)Nyy;
    mwSize Nz = (mwSize)Nzz;
    mwSize dim_array[4]={mm, Nx, Ny, Nz};
    blitz::Array<double, 4> data(m, Nx, Ny, Nz, fortranArray);
    data = wa;
    mat.matPut("wa",data.data(),data.size()*sizeof(double),4,dim_array,mxDOUBLE_CLASS,mxREAL);
    data = wb;
    mat.matPut("wb",data.data(),data.size()*sizeof(double),4,dim_array,mxDOUBLE_CLASS,mxREAL);
    data = wc;
    mat.matPut("wc",data.data(),data.size()*sizeof(double),4,dim_array,mxDOUBLE_CLASS,mxREAL);
    data = phia;
    mat.matPut("phia",data.data(),data.size()*sizeof(double),4,dim_array,mxDOUBLE_CLASS,mxREAL);
    data = phib;
    mat.matPut("phib",data.data(),data.size()*sizeof(double),4,dim_array,mxDOUBLE_CLASS,mxREAL);
    data = phic;
    mat.matPut("phic",data.data(),data.size()*sizeof(double),4,dim_array,mxDOUBLE_CLASS,mxREAL);
    mwSize dims1[1] = {(mwSize)maxT};
    mwSize dims2[1] = {(mwSize)m};
    mat.matPut("F", F.data(), F.size()*sizeof(double), 1, dims1, mxDOUBLE_CLASS, mxREAL);
    mat.matPut("S", S.data(), S.size()*sizeof(double), 1, dims2, mxDOUBLE_CLASS, mxREAL);
    mat.matPut("H", H.data(), H.size()*sizeof(double), 1, dims2, mxDOUBLE_CLASS, mxREAL);
    
    mat.matRelease();
}

void run_scftInString(string config_file) {  //scft calcuation along one string
	blitz::Range all = blitz::Range::all();
	//Config cfg(config_file.c_str());
	//Model *pmodel;
    //if(cfg.model() == ModelType::AB) {
    //    pmodel = new Model_AB(config_file);
    //}
	/*for (int s=0; s<m; s++) {
		cout << endl;
		cout << "**************************************************" << endl;
		cout << "scft running of " << s+1 << "th bead" << endl;
		cout << "**************************************************" << endl;
		cout << endl;
		pmodel->input_AField(wa(s,all,all,all));  //field initialization of single scft calculation
    	pmodel->input_BField(wb(s,all,all,all));
        pmodel->input_CField(wc(s,all,all,all));
    	pmodel->init_data_field();
    	scft sim(config_file, pmodel);
    	sim.run(); 
        H(s) = pmodel->H();

    	cout << endl;
    	cout << "saving data now ..." << endl;
    	cout << endl;
    	blitz::Array<double, 4> tmp(pmodel->output_data()); 
    	wa(s, all, all, all) = tmp(0, all, all, all);  //write the  scft results after 5 iterations into arrays
    	wb(s, all, all, all) = tmp(1, all, all, all);
        wc(s, all, all, all) = tmp(2, all, all, all);
    	phia(s, all, all, all) = tmp(3, all, all, all);
    	phib(s, all, all, all) = tmp(4, all, all, all);
        phic(s, all, all, all) = tmp(5, all, all, all);
	}*/
	Config cfg(config_file.c_str());
    Model *pmodel0;
    if(cfg.model() == ModelType::AB_C) {
        pmodel0 = new Model_AB_C(config_file);
    }
    Model *pmodel1;
    if(cfg.model() == ModelType::AB_C) {
        pmodel1 = new Model_AB_C(config_file);
    }
    Model *pmodel2;
    if(cfg.model() == ModelType::AB_C) {
        pmodel2 = new Model_AB_C(config_file);
    }
    Model *pmodel3;
    if(cfg.model() == ModelType::AB_C) {
        pmodel3 = new Model_AB_C(config_file);
    }
    Model *pmodel4;
    if(cfg.model() == ModelType::AB_C) {
        pmodel4 = new Model_AB_C(config_file);
    }
    Model *pmodel5;
    if(cfg.model() == ModelType::AB_C) {
        pmodel5 = new Model_AB_C(config_file);
    }
    Model *pmodel6;
    if(cfg.model() == ModelType::AB_C) {
        pmodel6 = new Model_AB_C(config_file);
    }
    Model *pmodel7;
    if(cfg.model() == ModelType::AB_C) {
        pmodel7 = new Model_AB_C(config_file);
    }
    Model *pmodel8;
    if(cfg.model() == ModelType::AB_C) {
        pmodel8 = new Model_AB_C(config_file);
    }
    Model *pmodel9;
    if(cfg.model() == ModelType::AB_C) {
        pmodel9 = new Model_AB_C(config_file);
    }
    Model *pmodel10;
    if(cfg.model() == ModelType::AB_C) {
        pmodel10 = new Model_AB_C(config_file);
    }
    Model *pmodel11;
    if(cfg.model() == ModelType::AB_C) {
        pmodel11 = new Model_AB_C(config_file);
    }
    Model *pmodel12;
    if(cfg.model() == ModelType::AB_C) {
        pmodel12 = new Model_AB_C(config_file);
    }
    Model *pmodel13;
    if(cfg.model() == ModelType::AB_C) {
        pmodel13 = new Model_AB_C(config_file);
    }
    Model *pmodel14;
    if(cfg.model() == ModelType::AB_C) {
        pmodel14 = new Model_AB_C(config_file);
    }
    Model *pmodel15;
    if(cfg.model() == ModelType::AB_C) {
        pmodel15 = new Model_AB_C(config_file);
    }
    Model *pmodel16;
    if(cfg.model() == ModelType::AB_C) {
        pmodel16 = new Model_AB_C(config_file);
    }
    Model *pmodel17;
    if(cfg.model() == ModelType::AB_C) {
        pmodel17 = new Model_AB_C(config_file);
    }
    Model *pmodel18;
    if(cfg.model() == ModelType::AB_C) {
        pmodel18 = new Model_AB_C(config_file);
    }
    Model *pmodel19;
    if(cfg.model() == ModelType::AB_C) {
        pmodel19 = new Model_AB_C(config_file);
    }
    Model *pmodel20;
    if(cfg.model() == ModelType::AB_C) {
        pmodel20 = new Model_AB_C(config_file);
    }

    string config0 = config_file;
    string config1 = config_file;
    string config2 = config_file;
    string config3 = config_file;
    string config4 = config_file;
    string config5 = config_file;
    string config6 = config_file;
    string config7 = config_file;
    string config8 = config_file;
    string config9 = config_file;
    string config10 = config_file;
    string config11 = config_file;
    string config12 = config_file;
    string config13 = config_file;
    string config14 = config_file;
    string config15 = config_file;
    string config16 = config_file;
    string config17 = config_file;
    string config18 = config_file;
    string config19 = config_file;
    string config20 = config_file;
    std::thread th[21];
    th[0] = std::thread(thread_scft, pmodel0, config0, wa, wb, wc, 0);
    th[1] = std::thread(thread_scft, pmodel1, config1, wa, wb, wc, 1);
    th[2] = std::thread(thread_scft, pmodel2, config2, wa, wb, wc, 2);
    th[3] = std::thread(thread_scft, pmodel3, config3, wa, wb, wc, 3);
    th[4] = std::thread(thread_scft, pmodel4, config4, wa, wb, wc, 4);
    th[5] = std::thread(thread_scft, pmodel5, config5, wa, wb, wc, 5);
    th[6] = std::thread(thread_scft, pmodel6, config6, wa, wb, wc, 6);
    th[7] = std::thread(thread_scft, pmodel7, config7, wa, wb, wc, 7);
    th[8] = std::thread(thread_scft, pmodel8, config8, wa, wb, wc, 8);
    th[9] = std::thread(thread_scft, pmodel9, config9, wa, wb, wc, 9);
    th[10] = std::thread(thread_scft, pmodel10, config10, wa, wb, wc, 10);
    th[11] = std::thread(thread_scft, pmodel11, config11, wa, wb, wc, 11);
    th[12] = std::thread(thread_scft, pmodel12, config12, wa, wb, wc, 12);
    th[13] = std::thread(thread_scft, pmodel13, config13, wa, wb, wc, 13);
    th[14] = std::thread(thread_scft, pmodel14, config14, wa, wb, wc, 14);
    th[15] = std::thread(thread_scft, pmodel15, config15, wa, wb, wc, 15);
    th[16] = std::thread(thread_scft, pmodel16, config16, wa, wb, wc, 16);
    th[17] = std::thread(thread_scft, pmodel17, config17, wa, wb, wc, 17);
    th[18] = std::thread(thread_scft, pmodel18, config18, wa, wb, wc, 18);
    th[19] = std::thread(thread_scft, pmodel19, config19, wa, wb, wc, 19);
    th[20] = std::thread(thread_scft, pmodel20, config20, wa, wb, wc, 20);
    th[0].join();
    th[1].join();
    th[2].join();
    th[3].join();
    th[4].join();
    th[5].join();
    th[6].join();
    th[7].join();
    th[8].join();
    th[9].join();
    th[10].join();
    th[11].join();
    th[12].join();
    th[13].join();
    th[14].join();
    th[15].join();
    th[16].join();
    th[17].join();
    th[18].join();
    th[19].join();
    th[20].join();
}

void thread_scft(Model *pmodel, string config_file, Array<double, 4> wa, Array<double, 4> wb, Array<double, 4> wc, int s) {
	blitz::Range all = blitz::Range::all();
	mtx.lock();
    pmodel->input_AField(wa(s, all, all, all));
    pmodel->input_BField(wb(s, all, all, all));
    pmodel->input_CField(wc(s, all, all, all));
    pmodel->init_data_field();
    scft sim(config_file, pmodel);
    mtx.unlock();
    sim.run();
    mtx.lock();
    H(s) = pmodel->H(); 
    blitz::Array<double, 4> tmp(pmodel->output_data()); 
    wa(s, all, all, all) = tmp(0, all, all, all);  //write the  scft results after 5 iterations into arrays
    wb(s, all, all, all) = tmp(1, all, all, all);
    wc(s, all, all, all) = tmp(2, all, all, all);
    phia(s, all, all, all) = tmp(3, all, all, all);
    phib(s, all, all, all) = tmp(4, all, all, all);
    phic(s, all, all, all) = tmp(5, all, all, all);
    mtx.unlock();
}

void parameterize_string() {
	blitz::Range all = blitz::Range::all();
	Config cfg("paramString.ini");
	UnitCell uc(cfg);
	double arcLength;
	//double sum = blitz::sum(H) - H(0)/2.0 - H(m)/2.0;
	//Array<double, 1> w(m-1);
	S(0) = 0;
	for (int s=1; s<m; s++) {
		blitz::Array<double, 3> tmp((wa(s, all, all, all)-wa(s-1, all, all, all)) * (wa(s, all, all, all)-wa(s-1, all, all, all)));
		Grid g(uc, Nxx, Nyy, Nzz, tmp);
		arcLength = sqrt(g.quadrature());
		//w(s) = (H(s-1)+H(s))/2.0/sum; //weighting function of arc length: w(s) = 2H/sum(2H), where 2H=(H(s-1)+H(s))/2.0
		//arcLength *= w(s);
        cout << "arcLength = " << arcLength << endl;
		S(s) = S(s-1) + arcLength;
	}
	S = 1.0*S/S(m-1);
}

void redistribute_stringBeads() {
	std::vector<double> X(m), Ya(m), Yb(m), Yc(m);  //for cublic spline interpolation
    tk::spline cs;  
    for (int i=0; i<Nxx; i++)
      for (int j=0; j<Nyy; j++)
        for (int k=0; k<Nzz; k++) {
          for (int z=0; z<m; z++) {
              X[z] = S(z);        //must converted to "std::vector<double>" stype.
              Ya[z] = wa(z, i, j, k); 
              Yb[z] = wb(z, i, j, k);
              Yc[z] = wc(z, i, j, k);
          }

          cs.set_points(X, Ya);
          for (int z=0; z<m; z++) {
            wa(z, i, j, k) = cs(1.0*z/(m-1));
          } 
          cs.set_points(X, Yb);
          for (int z=0; z<m; z++) {
            wb(z, i, j, k) = cs(1.0*z/(m-1));
          }
          cs.set_points(X, Yc);
          for (int z=0; z<m; z++) {
            wc(z, i, j, k) = cs(1.0*z/(m-1));
          }
        }
          
}