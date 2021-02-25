void peakproc()
{
  TFile *fBrem = TFile::Open("output126um_2_521MeV_1uCnormalized.root");
  TH2D  *hBrem = fBrem->Get("cosangleEnergyFar");
  
  double cosThetaMin = 0.995;
  double cosThetaMax = 1.0;
  int bMin = hBrem->GetYaxis()->FindBin(cosThetaMin);
  int bMax = hBrem->GetYaxis()->FindBin(cosThetaMax);
  int nCosThetaBins = bMax-bMin;

  TH1D *hFix = hBrem->ProjectionX("_hpx1",bMin, bMax);
  hFix->Sumw2();
  hFix->Scale(1.0/(1.0*nCosThetaBins));

  TH1D *hForward = hBrem->ProjectionX("_hpx2",bMax-1,bMax-1);
  hForward->Sumw2();
  hForward->Divide(hFix);
  hForward->Draw();

  vector<double> er;
  er.push_back(2.176);
  er.push_back(2.209);
  er.push_back(2.245);
  er.push_back(2.212107);

  double binsIntegral = 0.0;
  for (size_t i = 0; i < er.size(); ++i){
    int erBin = hBrem->GetXaxis()->FindBin(er[i]);
    binsIntegral += hBrem->Integral(erBin, erBin, bMin, bMax);
  }

  double binsSampled = 1.0 * nCosThetaBins * er.size(); // Number of bins included in the sample range
  double scaleFactor = 100.0e-6 / (1.0*hBrem->GetXaxis()->GetBinWidth(1)); // adjust for the fact that sampling bins are 100 eV while histo bins are 10 keV
  double bav = binsIntegral/binsSampled/scaleFactor;
  cout << binsSampled << endl;
  cout << nCosThetaBins << endl;
  cout << binsIntegral << endl;
  cout << bav << endl;

  double weightFactor = 1.0/bav;
  double inorm = binsIntegral;

  double gammaPerFile = 1.0e9;

  // Find the NRF reaction rate, stripped of efficiencies and outgoing attenuation, as a function of z
  TH1D *hRate2176 = new TH1D("hRate2176","hRate2176",35,0.0,3.5);
  hRate2176->Sumw2();

  // efficiency fits from efficiencyFits.root
  TString configName = "OpenBeam";
  TFile *outfile = new TFile("peakSpectrum_"+configName+"_1e9.root","RECREATE");
  TFile *effFile = TFile::Open("~/ZeroKnowledge/data/efficiencyFits.root");
  for (int detCode = 0; detCode < 3; ++detCode){
    TF1 *eff2176 = (TF1*) effFile->Get(Form("f_U_2176_%d",detCode));
    TF1 *eff2245 = (TF1*) effFile->Get(Form("f_U_2245_%d",detCode));
    TF1 *eff2212 = (TF1*) effFile->Get(Form("f_Al_2212_%d",detCode));

    gEnv->SetValue("TFile.Recover", 0);

    double eps = 1.0e-15;
    TH1D *peakSpectrum = new TH1D(Form("peakSpectrum%d",detCode),"Spectrum generated in detectors;energy [MeV];weighted counts [/mC]",500,2.0-eps,2.5-eps);
    peakSpectrum->Sumw2();

    TH1D *idealEnergy = new TH1D(Form("idealEnergy%d",detCode),"Spectrum incident on Det;energy [MeV];weighted counts [/mC]",2500,0.0-eps,2.5-eps);
    idealEnergy->Sumw2();

    TFile *infile = new TFile(configName+"_1e9.root");

    if (infile->IsZombie()) continue;   //file is unusable
    if (infile->TestBit(TFile::kRecovered)) continue;

    cout<<"Processing "<< detCode <<", Status: "<< infile->GetFd()<<"\n";

    if (infile->GetFd() < 0)
    {
      delete file;
      continue;
    }

    // Set up some way to compute the z-averaged detection probabilities for each line, compare to Brian's results
    double avgEff2176 = 0.0;
    double avgEff2245 = 0.0;
    double avgEff2212 = 0.0;

    double wsum2176 = 0.0;
    double wsum2245 = 0.0;
    double wsum2212 = 0.0;

    TTree *foilTree = (TTree*) infile->Get("foilTree");
    int foilTreeEntries = foilTree->GetEntries();

    double E = 0;
    double w = 0;
    double z = 0;
 
    foilTree->SetBranchAddress("foilEnergy",&E);
    foilTree->SetBranchAddress("foilWeight",&w);
    foilTree->SetBranchAddress("foilZ", &z);

    for (unsigned int entry = 0; entry < foilTreeEntries; ++entry)
    {
      foilTree->GetEntry(entry);

      double eff = 0.0;
      if (E>2.175 && E<2.177) {
        eff = eff2176->Eval(z);
        avgEff2176 += (w*weightFactor*eff); // w: bin weight in PGA sampling; weightFactor: effective bins sampled divided by their integral
        wsum2176   += (w*weightFactor);
        double ofactor = 1.0;//hForward->GetBinContent(hForward->FindBin(E));
        if (detCode == 0) hRate2176->Fill(z-997.2625, w*weightFactor/ofactor); // fill this only on first pass
      }
      else if (E>2.244 && E<2.246) {
        eff = eff2245->Eval(z);
        avgEff2245 += (w*weightFactor*eff);
        wsum2245   += (w*weightFactor);
      }
      else if (E>2.211 && E<2.213) {
        eff = eff2212->Eval(z);
        avgEff2212 += (w*weightFactor*eff);
        wsum2212   += (w*weightFactor);
      }
      else {
        //cout << "zeroing out energy E = " << E << endl;
        eff = 0;
      }
      peakSpectrum->Fill(E,w*weightFactor*eff);
    }

    TTree *idealDetTree = (TTree*) infile->Get("idealDetTree");
    if (idealDetTree){
      int idealDetTreeEntries = idealDetTree->GetEntries();
      double Eid = 0;
      double wid = 0;
      int idxid = 0;
      idealDetTree->SetBranchAddress("idealDetEnergy",&Eid);
      idealDetTree->SetBranchAddress("idealDetWeight",&wid);
      idealDetTree->SetBranchAddress("idealDetIndex",&idxid);
      for (unsigned int ee = 0; ee < idealDetTreeEntries; ++ee){
        idealDetTree->GetEntry(ee);
        if      (idxid==0 && detCode==0) idealEnergy->Fill(Eid,wid*weightFactor);
        else if (idxid==1 && detCode==1) idealEnergy->Fill(Eid,wid*weightFactor);
        else if (idxid==2 && detCode==2) idealEnergy->Fill(Eid,wid*weightFactor);
      }
    }

    delete infile;

    peakSpectrum->Scale(1.0/(gammaPerFile/inorm));
    peakSpectrum->Scale(1e3); // /uC to /mC
    peakSpectrum->Scale(1.0/0.932); // backscatter correction

    outfile->cd();
    peakSpectrum->Write();

    idealEnergy->Scale(1e3/0.932/(gammaPerFile/inorm));
    idealEnergy->Write();

    cout << "writing " << peakSpectrum->GetName() << " to " << outfile->GetName() << endl;
    outfile->Write();
    delete peakSpectrum;

    avgEff2176 /= wsum2176;
    avgEff2245 /= wsum2245;
    avgEff2212 /= wsum2212;

    // in Brian's det convention, 0 is bottom L (my 3), 1 is top L (my 0), and 2 is bottom R (my 2)
    cout << "average detection probabilities for detCode = " << detCode << ":" << endl;
    cout << Form("  2176 keV: %.3e", avgEff2176) << endl;
    cout << Form("  2245 keV: %.3e", avgEff2245) << endl;
    cout << Form("  2212 keV: %.3e", avgEff2212) << endl;


  } // end detCode loop

  new TCanvas();
  hRate2176->Scale(1.0/(gammaPerFile/inorm));
  hRate2176->Scale(1.0/0.932); // backscatter correction
  hRate2176->GetXaxis()->SetTitle("foil depth #it{z} [mm]");
  hRate2176->GetYaxis()->SetTitle(Form("2176 keV NRF reaction rate [/#muC/%1.3f mm bin]",hRate2176->GetBinWidth(1)));
  hRate2176->SetTitle(Form("rx rate comparison, %s: G4NRF (points) vs glass eq (curve)",configName.Data()));
  hRate2176->Draw("e");
  cout << "Integral: " << hRate2176->Integral() << " /uC" << endl;

  TGraph *g0 = new TGraph(Form("predictedRateVsZ%s2176.dat",configName.Data()), "%lg %lg");
  g0->Draw("same");

  // analyze the slopes
  TF1 *fh = new TF1("fh","expo",0.1,3.1);
  TF1 *fg = new TF1("fg","expo",0.1,3.1);

  hRate2176->Fit("fh","BRM");
  g0->Fit("fg","BRM");

  double sh = fh->GetParameter(1);
  double sg = fg->GetParameter(1);
  double dsh = fh->GetParError(1);
  double dsg = fg->GetParError(1);
  cout << "slope discrepancy = " << (sh-sg)/sqrt(dsh*dsh+dsg*dsg) << " sigma (ratio = " << sh/sg << ")" << endl;

  double ch = fh->GetParameter(0);
  double cg = fg->GetParameter(0);
  double dch = fh->GetParError(0);
  double dcg = fg->GetParError(0);
  cout << "const discrepancy = " << (ch-cg)/sqrt(dch*dch+dcg*dcg) << " sigma (ratio = " << ch/cg << ")" << endl;

  outfile->Close();
  return;
}
