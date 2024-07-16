//Change if you need
//#include "src/dofits/DoFit_Jpsi_Run.h"
//#include "src/dofits/DoFit_Jpsi_Run_2xGaus.h"
#include "src/dofits/DoFit_Jpsi_MC.h"
#include "src/dofits/DoFit_Jpsi_MC_2xGaus.h"

#include "src/create_folder.h"
#include "src/create_TH2D.h"
#include "src/get_efficiency_TH2D.h"
#include "src/yields_n_errs_to_TH2Ds_bin.h"

//Which Muon Id do you want to study?
string MuonId   = "trackerMuon";
//string MuonId   = "standaloneMuon";
//string MuonId   = "globalMuon";

// Bins to study
string xquantity = "Pt";
double xbins[] = {0.0, 3.4, 4.0, 5.0, 6.0, 8.0, 40.};
//double xbins[] = {0.0, 0.2125, 0.425, 0.6375, 0.85, 1.0625, 1.275, 1.4875, 1.7, 1.9125, 2.125, 2.3375, 2.55, 2.7625, 2.975, 3.1875, 3.4, 3.4375, 3.475, 3.5125, 3.55, 3.5875, 3.625, 3.6625, 3.7, 3.7375, 3.775, 3.8125, 3.85, 3.8875, 3.925, 3.9625, 4.0, 4.0625, 4.125, 4.1875, 4.25, 4.3125, 4.375, 4.4375, 4.5, 4.5625, 4.625, 4.6875, 4.75, 4.8125, 4.875, 4.9375, 5.0, 5.0625, 5.125, 5.1875, 5.25, 5.3125, 5.375, 5.4375, 5.5, 5.5625, 5.625, 5.6875, 5.75, 5.8125, 5.875, 5.9375, 6.0, 6.125, 6.25, 6.375, 6.5, 6.625, 6.75, 6.875, 7.0, 7.125, 7.25, 7.375, 7.5, 7.625, 7.75, 7.875, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0};

//double xbins[] = {4.0, 5.0, 6.0};
string yquantity = "Eta";
double ybins[] = {0.0, 0.4, 0.6, 0.95, 1.2, 1.4, 2.4};
//double ybins[] = {0.6, 0.95};

//Note: the y axis is absolute!

void plot_sys_efficiency_2d()
{
	//Path where is going to save fit results png for every bin 
	string path_bins_fit_folder = string("results/bins_fit/systematic_2D/") + output_folder_name + "/"+ MuonId + "/";
	create_folder(path_bins_fit_folder.c_str(), true);


	//Path where is going to save the efficiency results
	string directoryToSave = string("results/efficiencies/systematic_2D/") + output_folder_name + "/";
	create_folder(directoryToSave.c_str());

	//Get number of bins
	const int nbinsy = sizeof(ybins)/sizeof(*ybins) - 1;
	const int nbinsx = sizeof(xbins)/sizeof(*xbins) - 1;


	string file_path = directoryToSave + yquantity + "_" + xquantity + "_" + MuonId + ".root";
	TFile* generatedFile = new TFile(file_path.c_str(),"recreate");
	generatedFile->mkdir("histograms/");
	generatedFile->   cd("histograms/");

	TH2D *hist_all_nominal     = create_TH2D("all_nominal"   ,  "All Nominal",     xquantity, yquantity, nbinsx, nbinsy, xbins, ybins);
	TH2D *hist_all_2gauss      = create_TH2D("all_2xGauss"   ,  "All 2xGauss",     xquantity, yquantity, nbinsx, nbinsy, xbins, ybins);
	TH2D *hist_all_massup      = create_TH2D("all_MassUp"    ,  "All MassUp",      xquantity, yquantity, nbinsx, nbinsy, xbins, ybins);
	TH2D *hist_all_massdown    = create_TH2D("all_MassDown"  ,  "All MassDown",    xquantity, yquantity, nbinsx, nbinsy, xbins, ybins);
	TH2D *hist_all_binup       = create_TH2D("all_BinUp"     ,  "All BinUp",       xquantity, yquantity, nbinsx, nbinsy, xbins, ybins);
	TH2D *hist_all_bindown     = create_TH2D("all_BinDown"   ,  "All BinDown",     xquantity, yquantity, nbinsx, nbinsy, xbins, ybins);
	TH2D *hist_all_systematic  = create_TH2D("all_systematic",  "All Systematic",  xquantity, yquantity, nbinsx, nbinsy, xbins, ybins);
	TH2D *hist_all_final       = create_TH2D("all_final"     ,  "All Final",       xquantity, yquantity, nbinsx, nbinsy, xbins, ybins);

	TH2D *hist_pass_nominal    = create_TH2D("pass_nominal"   , "Pass Nominal",    xquantity, yquantity, nbinsx, nbinsy, xbins, ybins);
	TH2D *hist_pass_2gauss     = create_TH2D("pass_2xGauss"   , "Pass 2xGauss",    xquantity, yquantity, nbinsx, nbinsy, xbins, ybins);
	TH2D *hist_pass_massup     = create_TH2D("pass_MassUp"    , "Pass MassUp",     xquantity, yquantity, nbinsx, nbinsy, xbins, ybins);
	TH2D *hist_pass_massdown   = create_TH2D("pass_MassDown"  , "Pass MassDown",   xquantity, yquantity, nbinsx, nbinsy, xbins, ybins);
	TH2D *hist_pass_binup      = create_TH2D("pass_BinUp"     , "Pass BinUp",      xquantity, yquantity, nbinsx, nbinsy, xbins, ybins);
	TH2D *hist_pass_bindown    = create_TH2D("pass_BinDown"   , "Pass BinDown",    xquantity, yquantity, nbinsx, nbinsy, xbins, ybins);
	TH2D *hist_pass_systematic = create_TH2D("pass_systematic", "Pass Systematic", xquantity, yquantity, nbinsx, nbinsy, xbins, ybins);
	TH2D *hist_pass_final      = create_TH2D("pass_final"     , "Pass Final",      xquantity, yquantity, nbinsx, nbinsy, xbins, ybins);


	//Loop and fits
	for (int i = 0; i < nbinsx; i++)
	{
		for (int j = 0; j < nbinsy; j++)
		{
			//Creates conditions
			string conditions = string(    "ProbeMuon_" + xquantity + ">=" + to_string(xbins[i]  ));
			conditions +=       string(" && ProbeMuon_" + xquantity + "< " + to_string(xbins[i+1]));
			conditions +=       string(" && abs(ProbeMuon_" + yquantity + ")>=" + to_string(ybins[j]  ));
			conditions +=       string(" && abs(ProbeMuon_" + yquantity + ")< " + to_string(ybins[j+1]));

			const double default_min = _mmin;
			const double default_max = _mmax;
			string mmin_string;
			string mmax_string;
			double* yields_n_errs = NULL;
			double  yields_n_errs_systematic[4] = {0};
			double  yields_n_errs_final[4] = {0};

			//Nominal
			cout << "Nominal calculation -----\n";
			_mmin = default_min;
			_mmax = default_max;
			fit_bins = 100;
			prefix_file_name = "nominal_";
			yields_n_errs = doFit(conditions, MuonId, string(path_bins_fit_folder + prefix_file_name).c_str());
			yields_n_errs_systematic[0] =  yields_n_errs[0];
			yields_n_errs_systematic[1] =  yields_n_errs[1];
			yields_n_errs_final[0] =  yields_n_errs[0];
			yields_n_errs_final[1] =  yields_n_errs[1];
			yields_n_errs_final[2] =  yields_n_errs[2];
			yields_n_errs_final[3] =  yields_n_errs[3];
			yields_n_errs_to_TH2Ds_bin(hist_all_nominal, hist_pass_nominal, i+1, j+1, yields_n_errs);

			//2Gauss
			cout << "2xGassians calculation -----\n";
			_mmin = default_min;
			_mmax = default_max;
			fit_bins = 100;
			prefix_file_name = "2xgaus_";
			yields_n_errs = doFit2xGaus(conditions, MuonId, string(path_bins_fit_folder + prefix_file_name).c_str());
			yields_n_errs_systematic[2] += pow(yields_n_errs[2], 2);
			yields_n_errs_systematic[3] += pow(yields_n_errs[3], 2);
			yields_n_errs_to_TH2Ds_bin(hist_all_2gauss, hist_pass_2gauss, i+1, j+1, yields_n_errs);

			//MassUp
			cout << "MassUp calculation -----\n";
			_mmin = default_min - 0.05;
			_mmax = default_max + 0.05;
			fit_bins = 100;
			mmin_string = to_string(_mmin);
			mmax_string = to_string(_mmax);
			replace(mmin_string.begin(), mmin_string.end(), '.', 'p');
			replace(mmax_string.begin(), mmax_string.end(), '.', 'p');
			prefix_file_name  = string("mass_") + mmin_string.substr(0, mmin_string.length()-4) + string("_");
			prefix_file_name +=                   mmax_string.substr(0, mmax_string.length()-4) + string("_");
			yields_n_errs = doFit(conditions, MuonId, string(path_bins_fit_folder + prefix_file_name).c_str());
			yields_n_errs_systematic[2] += pow(yields_n_errs[2], 2);
			yields_n_errs_systematic[3] += pow(yields_n_errs[3], 2);
			yields_n_errs_to_TH2Ds_bin(hist_all_massup, hist_pass_massup, i+1, j+1, yields_n_errs);

			//MassDown
			cout << "MassDown calculation -----\n";
			_mmin = default_min + 0.05;
			_mmax = default_max - 0.05;
			fit_bins = 100;
			mmin_string = to_string(_mmin);
			mmax_string = to_string(_mmax);
			replace(mmin_string.begin(), mmin_string.end(), '.', 'p');
			replace(mmax_string.begin(), mmax_string.end(), '.', 'p');
			prefix_file_name  = string("mass_") + mmin_string.substr(0, mmin_string.length()-4) + string("_");
			prefix_file_name +=                   mmax_string.substr(0, mmax_string.length()-4) + string("_");
			yields_n_errs = doFit(conditions, MuonId, string(path_bins_fit_folder + prefix_file_name).c_str());
			yields_n_errs_systematic[2] += pow(yields_n_errs[2], 2);
			yields_n_errs_systematic[3] += pow(yields_n_errs[3], 2);
			yields_n_errs_to_TH2Ds_bin(hist_all_massdown, hist_pass_massdown, i+1, j+1, yields_n_errs);

			//BinUp
			cout << "BinUp calculation -----\n";
			_mmin = default_min;
			_mmax = default_max;
			fit_bins = 105;
			prefix_file_name = "binfit105_";
			yields_n_errs = doFit(conditions, MuonId, string(path_bins_fit_folder + prefix_file_name).c_str());
			yields_n_errs_systematic[2] += pow(yields_n_errs[2], 2);
			yields_n_errs_systematic[3] += pow(yields_n_errs[3], 2);
			yields_n_errs_to_TH2Ds_bin(hist_all_binup, hist_pass_binup, i+1, j+1, yields_n_errs);

			//BinDown
			cout << "BinDown calculation -----\n";
			_mmin = default_min;
			_mmax = default_max;
			fit_bins = 95;
			prefix_file_name = "binfit95_";
			yields_n_errs = doFit(conditions, MuonId, string(path_bins_fit_folder + prefix_file_name).c_str());
			yields_n_errs_systematic[2] += pow(yields_n_errs[2], 2);
			yields_n_errs_systematic[3] += pow(yields_n_errs[3], 2);
			yields_n_errs_to_TH2Ds_bin(hist_all_bindown, hist_pass_bindown, i+1, j+1, yields_n_errs);



			//Make the systematic calculations
			yields_n_errs_systematic[2] = sqrt(yields_n_errs_systematic[2]);
			yields_n_errs_systematic[3] = sqrt(yields_n_errs_systematic[3]);
			yields_n_errs_to_TH2Ds_bin(hist_all_systematic, hist_pass_systematic, i+1, j+1, yields_n_errs_systematic);


			//Make the statistic+systematic calculations
			yields_n_errs_final[2] = sqrt(pow(yields_n_errs_final[2], 2) + pow(yields_n_errs_systematic[2], 2));
			yields_n_errs_final[3] = sqrt(pow(yields_n_errs_final[3], 2) + pow(yields_n_errs_systematic[3], 2));
			yields_n_errs_to_TH2Ds_bin(hist_all_final, hist_pass_final, i+1, j+1, yields_n_errs_final);


			delete yields_n_errs;
		}
	}

	generatedFile->cd("/");
	get_efficiency_TH2D(hist_all_nominal,    hist_pass_nominal,    xquantity, yquantity, MuonId, "Nominal"   );
	get_efficiency_TH2D(hist_all_2gauss,     hist_pass_2gauss,     xquantity, yquantity, MuonId, "2xGauss"   );
	get_efficiency_TH2D(hist_all_massup,     hist_pass_massup,     xquantity, yquantity, MuonId, "MassUp"    );
	get_efficiency_TH2D(hist_all_massdown,   hist_pass_massdown,   xquantity, yquantity, MuonId, "MassDown"  );
	get_efficiency_TH2D(hist_all_binup,      hist_pass_binup,      xquantity, yquantity, MuonId, "BinUp"     );
	get_efficiency_TH2D(hist_all_bindown,    hist_pass_bindown,    xquantity, yquantity, MuonId, "BinDown"   );
	get_efficiency_TH2D(hist_all_systematic, hist_pass_systematic, xquantity, yquantity, MuonId, "Systematic");
	get_efficiency_TH2D(hist_all_final,      hist_pass_final,      xquantity, yquantity, MuonId, "Final"     );

	generatedFile->Write();

	cout << "\n------------------------\n";
	cout << "Output: " << file_path;
	cout << "\n------------------------\n";
}