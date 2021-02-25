// ROOT script to produce an importance sampling spectrum
// Jayson Vavrek, MIT, 2017
// jvavrek@mit.edu
{
	// first open the input brems spectrum in order to get energy limits, e.g.
	//TFile *f = (TFile*) TFile::Open("output126um_1uCnormalized.root");
	TFile *f = (TFile*) TFile::Open("output_2_7MeV_1uCnormalized.root");
	TH1D *ho = (TH1D*) f->Get("EnergyFar");

	double Emin = ho->GetXaxis()->GetXmin();
	double Emax = ho->GetXaxis()->GetXmax();

	// resonance energies in MeV as calculated by G4NRF
	vector<double> Evec;
	// U-238
	Evec.push_back(2.17601067909);
	Evec.push_back(2.20901100546);
	Evec.push_back(2.24501136709);
	// Al-27
	Evec.push_back(2.21210728352);

	double deltaE = 20.0e-6; // width of each important sampling region in MeV
	double binwidth = deltaE/4.0;

	int nbins = (Emax-Emin)/binwidth + 1;

	TH1D *hSample = new TH1D("hSample","hSample",nbins,Emin,Emax);
	TH1D *hBinary = new TH1D("hBinary","hBinary",nbins,Emin,Emax);

	// draw the spectrum
	for (int i = 0; i < nbins; ++i)
	{
		double e = hSample->GetBinCenter(i);
		for (int j = 0; j < Evec.size(); j++)
		{
			if (e < 2.0)
			{
				hSample->SetBinContent(i,0.0001);
				hBinary->SetBinContent(i,0.0);
			}
			else if (e > Evec[j] - deltaE/2.0 && e < Evec[j] + deltaE/2.0)
			{
				hSample->SetBinContent(i,1);
				hBinary->SetBinContent(i,1);
				break;
			}
			else
			{
				hSample->SetBinContent(i,0.01);
				hBinary->SetBinContent(i,0.0);
			}
		}
	}

	// normalize hSample so that its integral is 1
	hSample->Scale(1.0/(hSample->Integral()));

	//hSample->Draw();
	hSample->SaveAs("hSample.root");


	// also create a normalized brems spectrum for weighting
	ho->Smooth(TMath::Power(2,2));
	TH1D *hBremsWeight = new TH1D("hBremsWeight","hBremsWeight",nbins,Emin,Emax);
	for (int i = 1; i <= nbins; i++)
	{
		hBremsWeight->SetBinContent(i,ho->GetBinContent(ho->GetNbinsX()*(i-1)/nbins+1));
	}
	hBremsWeight->Scale(1.0/hBremsWeight->Integral());
	hBremsWeight->SaveAs("brem_normalized.root");

	TCanvas *c0 = new TCanvas();
	c0->cd();
	gPad->SetTicks(1,1);
	gPad->SetLogy();

	hSample->Draw();
	hBremsWeight->SetLineColor(kRed);
	hBremsWeight->Draw("same");

	//hBinary->SetLineColor(kBlack);
	//hBinary->Draw("same");
	hBinary->SaveAs("hBinary.root");

	hSample->GetYaxis()->SetRangeUser(1e-9,1e-1);
	hSample->SetTitle("NRF importance sampling distributions");
	hSample->GetXaxis()->SetTitle("energy #it{E} [MeV]");
	hSample->GetYaxis()->SetTitle(Form("probability per %2.2f eV",hSample->GetBinWidth(1)*1.0e6));
	hSample->SetStats(0);

	c0->SaveAs("weighting.png");

}
