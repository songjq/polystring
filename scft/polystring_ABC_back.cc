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
void parameterize_string(int);
void redistribute_stringBeads(int);
void save_data();
void run_scftInString(int, string);
Array<double, 3> input_data3(string, string);
Array<double, 4> input_data4(string, string);


int m; //number of beads on string.
int maxT; //max iteration of string.
Array<double, 2> S; //string nodes
Array<double, 2> H; //free energy of each node.
Array<double, 2> dH; //energy difference between two string iteration of each notes
Array<double, 1> F; //whole free energy of one string
Array<double, 4> wa, wb, wc, phia, phib, phic; //field and density data in each string node.
int Nxx;    //grid size
int Nyy;
int Nzz;
bool is_readString; //initialization from an existing string
double thresh_string_H; //thresh of energy to stop string iteration
Model *pmodel;

int main(int argc, char* argv[]){
	// the default configuration file is in the current working directory
    string config_file = "paramString.ini";
    Config cfg(config_file.c_str());
    m = cfg.get_integer("string","num_of_beads");
    maxT = cfg.get_integer("string","string_iteration");
    is_readString = cfg.get_bool("string", "is_readString");
    thresh_string_H = cfg.get_double("string", "thresh_string_H");
    S.resize(maxT, m);
    H.resize(maxT, m);
    dH.resize(maxT, m);
    dH = 10.0;
    F.resize(maxT);
    // Specify configuration file through input arguments.
    if(argc > 1){
        config_file = argv[1];
	}
    run_string(config_file);

    return 0;
}

void run_string(string config_file) {
	cout << "polystring_ABC_85\n";
	blitz::Range all = blitz::Range::all();
	Config cfg(config_file.c_str());
	cout << "polystring_ABC_88\n";
	if(cfg.model() == ModelType::AB_C) {
        pmodel = new Model_AB_C(config_file);  
    } //do not construct ETDRK4 frequently in string iterations, because it is time-consuming.
    int st = 0; //index of string iteration.
    cout << "polystring_ABC_93\n";
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
        run_scftInString(st, config_file);
        F(st) = mean(H(st, all));
        if(st>0)
            dH(st, all) = abs(H(st, all)-H(st-1, all));
        cout << "H = " << H(st, all) << endl;
        cout << "dH = " << dH(st, all) << endl;
        cout << "F = " << F(st) << endl;
        parameterize_string(st);
        save_data();
        cout << "S = " << S(st, all) << endl;
        redistribute_stringBeads(st);
        st += 1;
    }
    while(max(dH(st-1, all)) > thresh_string_H && st<maxT);
    
	
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
    mwSize dims2d[2] = {(mwSize)maxT, (mwSize)m};
    blitz::Array<double, 2> data2d(maxT, m, fortranArray);
    data2d = H;
    mat.matPut("H", data2d.data(), data2d.size()*sizeof(double), 2, dims2d, mxDOUBLE_CLASS, mxREAL);
    data2d = S;
    mat.matPut("S", data2d.data(), data2d.size()*sizeof(double), 2, dims2d, mxDOUBLE_CLASS, mxREAL);
    data2d = dH;
    mat.matPut("dH", data2d.data(), data2d.size()*sizeof(double), 2, dims2d, mxDOUBLE_CLASS, mxREAL);
    
    mat.matRelease();
}

void run_scftInString(int st, string config_file) {  //scft calcuation along one string
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
		cout << "polystring_ABC_245\n";
		pmodel->input_AField(wa(s,all,all,all));  //field initialization of single scft calculation
		cout << "polystring_ABC_247\n";
    	pmodel->input_BField(wb(s,all,all,all));
        pmodel->input_CField(wc(s,all,all,all));
    	pmodel->init_data_field();
    	cout << "polystring_ABC_251\n";
    	scft sim(config_file, pmodel);
    	cout << "polystring_ABC_253\n";
    	sim.run(); 
    	cout << "polystring_ABC_255\n";
        H(st, s) = pmodel->H();

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
	}
}

void parameterize_string(int st) {
	blitz::Range all = blitz::Range::all();
	Config cfg("paramString.ini");
	UnitCell uc(cfg);
	double arcLength;
	//double sum = blitz::sum(H) - H(0)/2.0 - H(m)/2.0;
	//Array<double, 1> w(m-1);
	S(st, 0) = 0;
	for (int s=1; s<m; s++) {
		blitz::Array<double, 3> tmp((wa(s, all, all, all)-wa(s-1, all, all, all)) * (wa(s, all, all, all)-wa(s-1, all, all, all)));
		Grid g(uc, Nxx, Nyy, Nzz, tmp);
		arcLength = sqrt(g.quadrature());
		//w(s) = (H(s-1)+H(s))/2.0/sum; //weighting function of arc length: w(s) = 2H/sum(2H), where 2H=(H(s-1)+H(s))/2.0
		//arcLength *= w(s);
        cout << "arcLength = " << arcLength << endl;
		S(st, s) = S(st, s-1) + arcLength;
	}
	S(st, all) = 1.0*S(st, all)/S(st, m-1);
}

void redistribute_stringBeads(int st) {
	std::vector<double> X(m), Ya(m), Yb(m), Yc(m);  //for cublic spline interpolation
    tk::spline cs;  
    for (int i=0; i<Nxx; i++)
      for (int j=0; j<Nyy; j++)
        for (int k=0; k<Nzz; k++) {
          for (int z=0; z<m; z++) {
              X[z] = S(st, z);        //must converted to "std::vector<double>" stype.
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