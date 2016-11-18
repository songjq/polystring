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
#include <vector>

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
//void thread_scft(string, Array<double, 3>, Array<double, 3>, int);
void thread_scft(Model *, string, Array<double, 3>, Array<double, 3>, int);
Array<double, 3> input_data3(string, string);
Array<double, 4> input_data4(string, string);

int m; //number of beads on string.
int maxT; //max iteration of string.
Array<double, 1> S; //string nodes.
Array<double, 1> H; //free energy of each node.
Array<double, 1> F; //whole free energy of one string
Array<double, 4> wa, wb, phia, phib; //field and density data in each string node.
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
	if(cfg.model() == ModelType::AB) {
        pmodel = new Model_AB(config_file);  
    } //do not construct ETDRK4 frequently in string iterations, because it is time-consuming.
    int st = 0; //index of string iteration.
	Nxx = cfg.Lx();
	Nyy = cfg.Ly();
	Nzz = cfg.Lz();
	wa.resize(m, Nxx, Nyy, Nzz);
	wb.resize(m, Nxx, Nyy, Nzz);
	phia.resize(m, Nxx, Nyy, Nzz);
	phib.resize(m, Nxx, Nyy, Nzz);

	initialize_string();
    do {
        cout << "**************************************************" << endl;
        cout << "This is the " << st << "th string iteration!" << endl;
        cout << "**************************************************" << endl;       
        run_scftInString(config_file);
        //cout << "Master thread: " << std::this_thread::get_id() << endl;
        //F(st) = mean(H);
        //cout << "H = " << H << endl;
        //cout << "F = " << F(st) << endl;
        //parameterize_string();
        //save_data();
        //cout << "S = " << S << endl;
        //redistribute_stringBeads();
        //save_data();
        st += 1;
    }
    //while(abs(F(st)-F(st-1))>1.0e-7 && st<maxT);
    while(st<maxT);
    
	
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
    	phia = input_data4("phia", "stringData.mat");
    	phib = input_data4("phib", "stringData.mat");
    }
    else {
    	cout << "Inilializing string by two ends!" << endl;
		/**************************initialize the starting bead and the end bead************************/
		wa(0, all, all, all) = input_data3("wA", "defect.mat");
    	wb(0, all, all, all) = input_data3("wB", "defect.mat");
    	phia(0, all, all, all) = input_data3("phiA", "defect.mat");
    	phib(0, all, all, all) = input_data3("phiB", "defect.mat");
    	wa(m-1, all, all, all) = input_data3("wA", "vLam.mat");
    	wb(m-1, all, all, all) = input_data3("wB", "vLam.mat");
    	phia(m-1, all, all, all) = input_data3("phiA", "vLam.mat");
    	phib(m-1, all, all, all) = input_data3("phiB", "vLam.mat");

    	/*************************initialize string for spinodal transition*********************/
    	for (int i=0; i<m; i++) {
    		wa(i, all, all, all) = wa(0, all, all, all) + 1.0*i/(m-1)* (wa(m-1, all, all, all) - wa(0, all, all, all));
    		wb(i, all, all, all) = wb(0, all, all, all) + 1.0*i/(m-1)* (wb(m-1, all, all, all) - wb(0, all, all, all));
        	phia(i, all, all, all) = phia(0, all, all, all) + 1.0*i/(m-1)* (phia(m-1, all, all, all) - phia(0, all, all, all));
        	phib(i, all, all, all) = phib(0, all, all, all) + 1.0*i/(m-1)* (phib(m-1, all, all, all) - phib(0, all, all, all));
    	}
    }
}

void save_data() {
    cout << "saving data now ..." << endl;
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
    data = phia;
    mat.matPut("phia",data.data(),data.size()*sizeof(double),4,dim_array,mxDOUBLE_CLASS,mxREAL);
    data = phib;
    mat.matPut("phib",data.data(),data.size()*sizeof(double),4,dim_array,mxDOUBLE_CLASS,mxREAL);
    cout << "saving field finished ..." << endl;
    mwSize dims1[1] = {(mwSize)maxT};
    mwSize dims2[1] = {(mwSize)m};
    mat.matPut("F", F.data(), F.size()*sizeof(double), 1, dims1, mxDOUBLE_CLASS, mxREAL);
    cout << "saving mean energy finished ..." << endl;
    mat.matPut("S", S.data(), S.size()*sizeof(double), 1, dims2, mxDOUBLE_CLASS, mxREAL);
    cout << "saving string finished ..." << endl;
    mat.matPut("H", H.data(), H.size()*sizeof(double), 1, dims2, mxDOUBLE_CLASS, mxREAL);
    cout << "saving data finished ..." << endl;
    
    mat.matRelease();
}

/*void run_scftInString(string config_file) {  //scft calcuation along one string
	blitz::Range all = blitz::Range::all();
	//Config cfg(config_file.c_str());
	//Model *pmodel;
    //if(cfg.model() == ModelType::AB) {
    //    pmodel = new Model_AB(config_file);
    //}
	for (int s=0; s<m; s++) {
		cout << endl;
		cout << "**************************************************" << endl;
		cout << "scft running of " << s+1 << "th bead" << endl;
		cout << "**************************************************" << endl;
		cout << endl;
		pmodel->input_AField(wa(s,all,all,all));  //field initialization of single scft calculation
    	pmodel->input_BField(wb(s,all,all,all));
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
    	phia(s, all, all, all) = tmp(2, all, all, all);
    	phib(s, all, all, all) = tmp(3, all, all, all);
	}
}*/
void run_scftInString(string config_file) {
    /*std::vector<std::thread> th;

    for (int i=0; i<m; i++) {
        th.push_back(std::thread(thread_scft, config_file, wa, wb, phia, phib, i));
    }

    for (auto &t : th) {
        t.join();
    }*/
    /*int m = 2;
    std::thread th[m];
    for (int i=0; i<m; i++) {
        th[i] = std::thread(thread_scft, config_file, wa, wb, phia, phib, i);
    }
    for (int i=0; i<m; i++) {
        th[i].join();
    }*/
    blitz::Range all = blitz::Range::all();
    Array<double, 3> wA0(wa(0,all,all,all));
    Array<double, 3> wB0(wb(0,all,all,all));
    Array<double, 3> wA1(wa(1,all,all,all));
    Array<double, 3> wB1(wb(1,all,all,all));
    Array<double, 3> wA2(wa(2,all,all,all));
    Array<double, 3> wB2(wb(2,all,all,all));
    Array<double, 3> wA3(wa(3,all,all,all));
    Array<double, 3> wB3(wb(3,all,all,all));
    Array<double, 3> wA4(wa(4,all,all,all));
    Array<double, 3> wB4(wb(4,all,all,all));
    Array<double, 3> wA5(wa(5,all,all,all));
    Array<double, 3> wB5(wb(5,all,all,all));
    Array<double, 3> wA6(wa(6,all,all,all));
    Array<double, 3> wB6(wb(6,all,all,all));
    Array<double, 3> wA7(wa(7,all,all,all));
    Array<double, 3> wB7(wb(7,all,all,all));
    Array<double, 3> wA8(wa(8,all,all,all));
    Array<double, 3> wB8(wb(8,all,all,all));
    Array<double, 3> wA9(wa(9,all,all,all));
    Array<double, 3> wB9(wb(9,all,all,all));
    Array<double, 3> wA10(wa(10,all,all,all));
    Array<double, 3> wB10(wb(10,all,all,all));
    Array<double, 3> wA11(wa(11,all,all,all));
    Array<double, 3> wB11(wb(11,all,all,all));
    Array<double, 3> wA12(wa(12,all,all,all));
    Array<double, 3> wB12(wb(12,all,all,all));
    Array<double, 3> wA13(wa(13,all,all,all));
    Array<double, 3> wB13(wb(13,all,all,all));
    Array<double, 3> wA14(wa(14,all,all,all));
    Array<double, 3> wB14(wb(14,all,all,all));
    Array<double, 3> wA15(wa(15,all,all,all));
    Array<double, 3> wB15(wb(15,all,all,all));
    Array<double, 3> wA16(wa(16,all,all,all));
    Array<double, 3> wB16(wb(16,all,all,all));
    Array<double, 3> wA17(wa(17,all,all,all));
    Array<double, 3> wB17(wb(17,all,all,all));
    Array<double, 3> wA18(wa(18,all,all,all));
    Array<double, 3> wB18(wb(18,all,all,all));
    Array<double, 3> wA19(wa(19,all,all,all));
    Array<double, 3> wB19(wb(19,all,all,all));
    Array<double, 3> wA20(wa(20,all,all,all));
    Array<double, 3> wB20(wb(20,all,all,all));
    Config cfg(config_file.c_str());
    Model *pmodel0;
    if(cfg.model() == ModelType::AB) {
        pmodel0 = new Model_AB(config_file);
    }
    Model *pmodel1;
    if(cfg.model() == ModelType::AB) {
        pmodel1 = new Model_AB(config_file);
    }
    Model *pmodel2;
    if(cfg.model() == ModelType::AB) {
        pmodel2 = new Model_AB(config_file);
    }
    Model *pmodel3;
    if(cfg.model() == ModelType::AB) {
        pmodel3 = new Model_AB(config_file);
    }
    Model *pmodel4;
    if(cfg.model() == ModelType::AB) {
        pmodel4 = new Model_AB(config_file);
    }
    Model *pmodel5;
    if(cfg.model() == ModelType::AB) {
        pmodel5 = new Model_AB(config_file);
    }
    Model *pmodel6;
    if(cfg.model() == ModelType::AB) {
        pmodel6 = new Model_AB(config_file);
    }
    Model *pmodel7;
    if(cfg.model() == ModelType::AB) {
        pmodel7 = new Model_AB(config_file);
    }
    Model *pmodel8;
    if(cfg.model() == ModelType::AB) {
        pmodel8 = new Model_AB(config_file);
    }
    Model *pmodel9;
    if(cfg.model() == ModelType::AB) {
        pmodel9 = new Model_AB(config_file);
    }
    Model *pmodel10;
    if(cfg.model() == ModelType::AB) {
        pmodel10 = new Model_AB(config_file);
    }
    /*Model *pmodel11;
    if(cfg.model() == ModelType::AB) {
        pmodel11 = new Model_AB(config_file);
    }
    Model *pmodel12;
    if(cfg.model() == ModelType::AB) {
        pmodel12 = new Model_AB(config_file);
    }
    Model *pmodel13;
    if(cfg.model() == ModelType::AB) {
        pmodel13 = new Model_AB(config_file);
    }
    Model *pmodel14;
    if(cfg.model() == ModelType::AB) {
        pmodel14 = new Model_AB(config_file);
    }
    Model *pmodel15;
    if(cfg.model() == ModelType::AB) {
        pmodel15 = new Model_AB(config_file);
    }
    Model *pmodel16;
    if(cfg.model() == ModelType::AB) {
        pmodel16 = new Model_AB(config_file);
    }
    Model *pmodel17;
    if(cfg.model() == ModelType::AB) {
        pmodel17 = new Model_AB(config_file);
    }
    Model *pmodel18;
    if(cfg.model() == ModelType::AB) {
        pmodel18 = new Model_AB(config_file);
    }
    Model *pmodel19;
    if(cfg.model() == ModelType::AB) {
        pmodel19 = new Model_AB(config_file);
    }
    Model *pmodel20;
    if(cfg.model() == ModelType::AB) {
        pmodel20 = new Model_AB(config_file);
    }*/
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
    /*string config11 = config_file;
    string config12 = config_file;
    string config13 = config_file;
    string config14 = config_file;
    string config15 = config_file;
    string config16 = config_file;
    string config17 = config_file;
    string config18 = config_file;
    string config19 = config_file;
    string config20 = config_file;*/
    std::thread th[11];
    th[0] = std::thread(thread_scft, pmodel0, config0, wA0, wB0, 0);
    th[1] = std::thread(thread_scft, pmodel0, config1, wA1, wB1, 1);
    th[2] = std::thread(thread_scft, pmodel0, config2, wA2, wB2, 2);
    th[3] = std::thread(thread_scft, pmodel0, config3, wA3, wB3, 3);
    th[4] = std::thread(thread_scft, pmodel0, config4, wA4, wB4, 4);
    th[5] = std::thread(thread_scft, pmodel0, config5, wA5, wB5, 5);
    th[6] = std::thread(thread_scft, pmodel0, config6, wA6, wB6, 6);
    th[7] = std::thread(thread_scft, pmodel0, config7, wA7, wB7, 7);
    th[8] = std::thread(thread_scft, pmodel0, config8, wA8, wB8, 8);
    th[9] = std::thread(thread_scft, pmodel0, config9, wA9, wB9, 9);
    th[10] = std::thread(thread_scft, pmodel0, config10, wA10, wB10, 10);
    /*th[11] = std::thread(thread_scft, pmodel11, config11, wA11, wB11, 11);
    ==th[12] = std::thread(thread_scft, pmodel12, config12, wA12, wB12, 12);
    th[13] = std::thread(thread_scft, pmodel13, config13, wA13, wB13, 13);
    th[14] = std::thread(thread_scft, pmodel14, config14, wA14, wB14, 14);
    th[15] = std::thread(thread_scft, pmodel15, config15, wA15, wB15, 15);
    th[16] = std::thread(thread_scft, pmodel16, config16, wA16, wB16, 16);
    th[17] = std::thread(thread_scft, pmodel17, config17, wA17, wB17, 17);
    th[18] = std::thread(thread_scft, pmodel18, config18, wA18, wB18, 18);
    th[19] = std::thread(thread_scft, pmodel19, config19, wA19, wB19, 19);
    th[20] = std::thread(thread_scft, pmodel20, config20, wA20, wB20, 20);*/
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
    /*th[11].join();
    th[12].join();
    th[13].join();
    th[14].join();
    th[15].join();
    th[16].join();
    th[17].join();
    th[18].join();
    th[19].join();
    th[20].join();*/
}

void thread_scft(Model *pmodel, string config_file, Array<double, 3> wA, Array<double, 3> wB, int s) {
    blitz::Range all = blitz::Range::all();
    /*Config cfg(config_file.c_str());
    Model *pmodel;
    if(cfg.model() == ModelType::AB) {
        pmodel = new Model_AB(config_file);
    }*/
    
    cout << endl;
    cout << "**************************************************" << endl;
    cout << "scft running of " << s+1 << "th bead" << endl;
    cout << "**************************************************" << endl;
    cout << endl;
    //pmodel->input_AField(wa(s,all,all,all));  //field initialization of single scft calculation
    //pmodel->input_BField(wb(s,all,all,all));
    
    mtx.lock();
    pmodel->input_AField(wA);
    pmodel->input_BField(wB);
    pmodel->init_data_field();
    scft sim(config_file, pmodel);
    sim.run(); 
    mtx.unlock();

     

    mtx.lock();
    //std::lock_guard<std::mutex> guard(mtx);
    H(s) = pmodel->H(); 
    //mtx.unlock();
    blitz::Array<double, 4> tmp(pmodel->output_data());  
    wa(s, all, all, all) = tmp(0, all, all, all);  //write the  scft results after 5 iterations into arrays
    wb(s, all, all, all) = tmp(1, all, all, all);
    phia(s, all, all, all) = tmp(2, all, all, all);
    phib(s, all, all, all) = tmp(3, all, all, all);
    mtx.unlock();
}

void parameterize_string() {
    cout << "parameterize string now ..." << endl;
	blitz::Range all = blitz::Range::all();
	//Config cfg("paramString.ini");
	//UnitCell uc(cfg);
	double arcLength;
	S(0) = 0;
	for (int s=1; s<m; s++) {
		blitz::Array<double, 3> tmp((wa(s, all, all, all)-wa(s-1, all, all, all)) * (wa(s, all, all, all)-wa(s-1, all, all, all)));
        arcLength = sqrt(mean(tmp));
		//Grid g(uc, Nxx, Nyy, Nzz, tmp);
		//arcLength = sqrt(g.quadrature());
		cout << "arcLength = " << arcLength << endl;
		S(s) = S(s-1) + arcLength;
	}
	S = 1.0*S/S(m-1);
	cout << "S = " << S << endl;
}

void redistribute_stringBeads() {
    cout << "natural cublic spline interpolation now ..." << endl;
    cout << "S is " << S << endl;
	std::vector<double> X(m), Ya(m), Yb(m);  //for cublic spline interpolation
    tk::spline cs;  
    cout << "interpolate field" << endl;
    for (int i=0; i<Nxx; i++)
      for (int j=0; j<Nyy; j++)
        for (int k=0; k<Nzz; k++) {
          for (int z=0; z<m; z++) {
              X[z] = S(z);        //must converted to "std::vector<double>" stype.
              Ya[z] = wa(z, i, j, k); 
              Yb[z] = wb(z, i, j, k);
          }
  
          cs.set_points(X, Ya);
          for (int z=0; z<m; z++) {
            wa(z, i, j, k) = cs(1.0*z/(m-1));
          } 
          cs.set_points(X, Yb);
          for (int z=0; z<m; z++) {
            wb(z, i, j, k) = cs(1.0*z/(m-1));
          }
        }
          
}