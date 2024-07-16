#include <iostream>
#include <TTree.h>
#include <vector>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TH1S.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooExponential.h>
#include <RooFitResult.h>
#include <TPaveText.h>

using namespace std;
using namespace RooFit;
using namespace RooStats;

const double LUMINOSITY = 1000; // 1pb-1 = 1nb-1
const double EFFICIENCY = 0.9856;
const double ACCEPTANCE = 0.230165;

void fitMassHistogram(TH1S* hMass, double& yield, double& yieldError) {
    RooRealVar mass("mass", "#mu^{+}#mu^{-} invariant mass", 2.9, 3.3, "GeV/c^{2}");
    RooDataHist dh("dh", "dh", mass, Import(*hMass));

    // Define background model (exponential) and its parameter
    RooRealVar lambda("lambda", "lambda", -0.3, -4., 0.);
    RooExponential background("background", "background", mass, lambda);

    // Define signal model (Gaussian) and its parameters
    RooRealVar mean("mean", "mean", 0.5 * (2.9 + 3.3), 2.9, 3.3);
    RooRealVar sigma("sigma", "sigma", 0.1 * (3.3 - 2.9), 0., 0.5 * (3.3 - 2.9));
    RooGaussian gaussian("gaussian", "gaussian", mass, mean, sigma);

    // Define signal model (Crystal Ball) and its parameters
    RooRealVar sigma_cb("sigma_cb", "sigma_cb", 0.06, 0.01, 0.2);
    RooRealVar alpha("alpha", "alpha", 2.5, 0.0, 6.0);
    RooRealVar n("n", "n", 3.98, 1, 10);
    RooCBShape crystalball("crystalball", "crystalball", mass, mean, sigma_cb, alpha, n);

    // Combine Gaussian and Crystal Ball models
    RooRealVar frac1("frac1", "frac1", 0.3);
    RooAddPdf signal("signal", "signal", RooArgList(gaussian, crystalball), RooArgList(frac1));

    // Define variables for number of signal and background events
    double n_signal_initial = 0.8 * dh.sumEntries();
    double n_back_initial = 0.2 * dh.sumEntries();
    RooRealVar n_signal("n_signal", "n_signal", n_signal_initial, 0., dh.sumEntries());
    RooRealVar n_back("n_back", "n_back", n_back_initial, 0., dh.sumEntries());

    // Sum signal and background models
    RooAddPdf model("model", "model", RooArgList(signal, background), RooArgList(n_signal, n_back));

    model.fitTo(dh);
    yield = n_signal.getVal();
    yieldError = n_signal.getError();
}

void calculateDifferentialCrossSection(const vector<double>& yields, const vector<double>& pt_bins, vector<double>& diff_cross_sections, vector<double>& errors) {
    cout << "Seção de choque diferencial (nb/GeV) por intervalo de pt:" << endl;
    for (size_t i = 0; i < yields.size(); ++i) {
        double delta_pt = pt_bins[i + 1] - pt_bins[i];
        double differential_cross_section = (yields[i] / (LUMINOSITY * EFFICIENCY * ACCEPTANCE * delta_pt)); 
        double error = (sqrt(yields[i]) / (LUMINOSITY * EFFICIENCY * ACCEPTANCE * delta_pt)); 
        diff_cross_sections.push_back(differential_cross_section);
        errors.push_back(error);
        cout << "Intervalo de pt " << pt_bins[i] << " - " << pt_bins[i + 1] << " GeV: "
             << differential_cross_section << " ± " << error << " nb/GeV" << endl;
    }
}

void plotDifferentialCrossSection(const vector<double>& pt_bins, const vector<double>& diff_cross_sections, const vector<double>& errors) {
    TGraphErrors *graph = new TGraphErrors(diff_cross_sections.size());
    for (size_t i = 0; i < diff_cross_sections.size(); ++i) {
        double pt_center = (pt_bins[i] + pt_bins[i + 1]) / 2.0;
        double delta_pt = pt_bins[i + 1] - pt_bins[i];
        graph->SetPoint(i, pt_center, diff_cross_sections[i]);
        graph->SetPointError(i, delta_pt / 2.0, errors[i]);
    }

    TCanvas *c1 = new TCanvas("c_diff_cross_section", "Seção de Choque Diferencial", 800, 600);
    graph->SetTitle(";pt (GeV);d#sigma/dpt (nb/GeV)");
    graph->GetXaxis()->SetTitleSize(0.04);
    graph->GetYaxis()->SetTitleSize(0.04);
    graph->Draw("AP");

    c1->SaveAs("differential_cross_section_gauss.png");

    //delete graph;
    //delete c;
}

void analiseCSGauss2() {
    Double_t event;   
    TLorentzVector *dimuon_p4 = nullptr;
    TLorentzVector *muonP_p4 = nullptr;
    TLorentzVector *muonN_p4 = nullptr;
   
    TBranch *b_event = nullptr;
    TBranch *b_dimuon_p4 = nullptr;
    TBranch *b_muonP_p4 = nullptr;
    TBranch *b_muonN_p4 = nullptr;
    
    // Histograma principal
    TH1S *hJpsiMassCut = new TH1S("JpsiMass", "", 200, 2.9, 3.3); 


    TChain* fChain = new TChain("oniaTree;4", "");
    fChain->Add("Skim4.root");

    fChain->SetBranchAddress("event", &event, &b_event);
    fChain->SetBranchAddress("dimuon_p4", &dimuon_p4, &b_dimuon_p4);
    fChain->SetBranchAddress("muonP_p4", &muonP_p4, &b_muonP_p4);
    fChain->SetBranchAddress("muonN_p4", &muonN_p4, &b_muonN_p4);

    int nentries = fChain->GetEntries();	
    int n_pass_cuts = 0;

    cout << "Numero de entradas " << nentries << endl;

    // Vetores para armazenar os yields e erros por intervalo de pt
    vector<double> yields(11, 0.0);
    vector<double> yieldErrors(11, 0.0);
    vector<double> pt_bins = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60};
 
    // Contadores para diferentes faixas de pt
    vector<int> pt_counts(11, 0);
    vector<string> pt_labels = {"0-5", "5-10", "10-15", "15-20", "20-25", "25-30", "30-35", "35-40", "40-45", "45-50", "50+"};

    // Criar histogramas para cada intervalo de pt
    vector<TH1S*> hists;
    for (int i = 0; i < 11; ++i) {
        hists.push_back(new TH1S(Form("hJpsiMassCut_pt_%s", pt_labels[i].c_str()), "", 200, 2.9, 3.3));
    }

    for(int i = 0; i < nentries; i++) {
        fChain->GetEntry(i);
		
        double etacut = abs(dimuon_p4->Eta());
        double mass = dimuon_p4->M();
        double pt1 = muonN_p4->Pt();
        double pt2 = muonP_p4->Pt();
        double pt3 = dimuon_p4->Pt();

        if(pt1 > 1 && pt2 > 1 && etacut < 2.4 && mass >= 2.9 && mass <= 3.3) {
            n_pass_cuts++;
            hJpsiMassCut->Fill(mass); // Preencher o histograma principal

            // Preencher o histograma correto para o intervalo de pt
            if (pt3 >= 0 && pt3 < 5) {
                hists[0]->Fill(mass);
                pt_counts[0]++;
            } else if (pt3 >= 5 && pt3 < 10) {
                hists[1]->Fill(mass);
                pt_counts[1]++;
            } else if (pt3 >= 10 && pt3 < 15) {
                hists[2]->Fill(mass);
                pt_counts[2]++;
            } else if (pt3 >= 15 && pt3 < 20) {
                hists[3]->Fill(mass);
                pt_counts[3]++;
            } else if (pt3 >= 20 && pt3 < 25) {
                hists[4]->Fill(mass);
                pt_counts[4]++;
            } else if (pt3 >= 25 && pt3 < 30) {
                hists[5]->Fill(mass);
                pt_counts[5]++;
            } else if (pt3 >= 30 && pt3 < 35) {
                hists[6]->Fill(mass);
                pt_counts[6]++;
            } else if (pt3 >= 35 && pt3 < 40) {
                hists[7]->Fill(mass);
                pt_counts[7]++;
            } else if (pt3 >= 40 && pt3 < 45) {
                hists[8]->Fill(mass);
                pt_counts[8]++;
            } else if (pt3 >= 45 && pt3 < 50) {
                hists[9]->Fill(mass);
                pt_counts[9]++;
            } else if (pt3 >= 50) {
                hists[10]->Fill(mass);
                pt_counts[10]++;
            }
        }
    }

    // Ajuste da distribuição de massa para cada intervalo de pt
    for (int i = 0; i < 11; i++) {
        if (pt_counts[i] > 0) {
            double yield, yieldError;
            fitMassHistogram(hists[i], yield, yieldError);
            yields[i] = yield;
            yieldErrors[i] = yieldError;
        }
    }

    // Calcular a seção de choque diferencial
    vector<double> diff_cross_sections;
    vector<double> errors;
    calculateDifferentialCrossSection(yields, pt_bins, diff_cross_sections, errors);

    // Plotar a seção de choque diferencial
    plotDifferentialCrossSection(pt_bins, diff_cross_sections, errors);

    // Verificar se o arquivo associado ao TChain foi carregado corretamente
    TFile* f = fChain->GetFile();
    if (f) {
        f->Close();
    } else {
        cerr << "Error: Could not close the file because it was not loaded properly." << endl;
    }

    TCanvas *c = new TCanvas("c","c",5,5,1000,1000);
    c->SetLeftMargin(0.15); // Aumenta a margem esquerda para evitar o corte do eixo Y

    hJpsiMassCut->GetXaxis()->SetTitle("m_{J/#psi} (GeV)");
    hJpsiMassCut->GetYaxis()->SetTitle("Events");
    hJpsiMassCut->SetMarkerColor(1);
    hJpsiMassCut->SetMarkerStyle(2);

    hJpsiMassCut->Draw();
 


    RooRealVar mass("mass", "#mu^{+}#mu^{-} invariant mass ", 2.9, 3.3, "GeV/c^{2}");

    // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'mass'
    RooDataHist dh("dh", "dh", mass, Import(*hJpsiMassCut));

    // Define background model (exponential) and its parameter
    RooRealVar lambda("lambda", "lambda", -0.3, -4., 0.);
    RooExponential background("background", "background", mass, lambda);

    // Define signal model (Gaussian) and its parameters
    RooRealVar mean("mean", "mean", 0.5 * (2.9 + 3.3), 2.9, 3.3);
    RooRealVar sigma("sigma", "sigma", 0.1 * (3.3 - 2.9), 0., 0.5 * (3.3 - 2.9));
    RooGaussian gaussian("gaussian", "gaussian", mass, mean, sigma);

    // Define signal model (Crystal Ball) and its parameters
    RooRealVar sigma_cb("sigma_cb", "sigma_cb", 0.038);
    RooRealVar alpha("alpha", "alpha", 1.71);
    RooRealVar n("n", "n", 3.96);
    RooCBShape crystalball("crystalball", "crystalball", mass, mean, sigma_cb, alpha, n);

    // Combine Gaussian and Crystal Ball models 3
    RooRealVar frac1("frac1", "frac1", 0.5);
    RooAddPdf signal("signal", "signal", RooArgList(gaussian, crystalball), RooArgList(frac1));

    // Define variables for number of signal and background events
    double n_signal_initial = 0.8 * dh.sumEntries();
    double n_back_initial = 0.2 * dh.sumEntries();
    RooRealVar n_signal("n_signal", "n_signal", n_signal_initial, 0., dh.sumEntries());
    RooRealVar n_back("n_back", "n_back", n_back_initial, 0., dh.sumEntries());

    // Sum signal and background models
    RooAddPdf model("model", "model", RooArgList(signal, background), RooArgList(n_signal, n_back));

    model.fitTo(dh);

    RooFitResult* fitresult = model.fitTo(dh, Extended(true), Save());

    RooPlot* frame = mass.frame();
    frame->SetTitle(""); // Remove o título da RooPlot
    dh.plotOn(frame);
    //model.plotOn(frame);

    double n_param, reduced_chi_square;
    n_param = fitresult->floatParsFinal().getSize();
    reduced_chi_square = frame->chiSquare(n_param)/2.86;
    //model.plotOn(frame, Components(signal), LineStyle(kDashed), LineColor(kGreen));
   // model.plotOn(frame, Components(background), LineStyle(kDashed), LineColor(kRed));
    frame->Draw();

  TLegend *legend=new TLegend(0.7,0.65,0.88,0.85);
  legend->SetBorderSize(0);
  legend->SetTextFont(40);
  legend->SetTextSize(0.03);
  legend->AddEntry(hJpsiMassCut,"Data","lep");
  legend->Draw("same");

    // TPaveText* paveText = new TPaveText(0.6, 0.7, 0.8, 0.8, "brNDC"); // Ajuste de posição e tamanho
    // paveText->SetFillColor(kWhite);
    // paveText->SetFillStyle(1001);
    // paveText->SetTextSize(0.015);
    // paveText->AddText(Form("#chi^{2}/ndof = %10f ", reduced_chi_square));

    // // paveText->AddText(Form("Mean_{m_{J/#psi}} = %.5f #pm %.5f GeV", mean.getVal(), mean.getError()));
    // // paveText->AddText(Form("#sigma_{m_{J/#psi}} = %.5f #pm %.5f GeV", sigma.getVal(), sigma.getError()));
    // // paveText->AddText(Form("lambda = %.5f #pm %.5f", lambda.getVal(), lambda.getError()));
    // // paveText->AddText(Form("n_signal = %.5f #pm %.5f", n_signal.getVal(), n_signal.getError()));

    // paveText->AddText(Form("Signal = %0.0f #pm %0.0f", n_signal.getVal(), n_signal.getError()));
    // paveText->AddText(Form("Fundo = %0.0f #pm %0.0f", n_back.getVal(), n_back.getError()));
    // paveText->Draw();

    TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.03); // Tamanho do texto
  latex.SetTextFont(42); // Fonte normal (sem negrito)
  latex.DrawLatex(0.15, 0.9, "#bf{CMS} #it{Open Data}");
  latex.SetTextSize(0.02);
  latex.DrawLatex(0.25,0.73,"Resonance: J/#psi");
  latex.DrawLatex(0.25,0.70,"p_{T} > 1 GeV");
  latex.DrawLatex(0.25,0.67,"|#eta| < 2,4");
  latex.SetTextSize(0.1);
  latex.DrawLatex(0.30, 0.9, ".");








    c->SaveAs("resultado_gaussiana.png");

    // Exibindo os yields por intervalo de pt
    for (int i = 0; i < 11; i++) {
        cout << "Yield em pt " << pt_labels[i] << " GeV: " << yields[i] << " ± " << yieldErrors[i] << endl;
    }

    // Limpeza dos histogramas criados
    for (auto hist : hists) {
        delete hist;
    }
}

int main() {
    analiseCSGauss2();
    return 0;
}
