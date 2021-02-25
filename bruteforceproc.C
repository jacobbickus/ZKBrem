void bruteforceproc()
{
  // weighting calculations
  TFile *fBrem = (TFile*)TFile::Open("output126um_2_521MeV_1uCnormalized.root");
  TH2D  *hBrem = (TH2D*) fBrem->Get("cosangleEnergyFar");
  
  double cosThetaMin = 0.995; // standard is 0.995; narrow is 0.9999
  double cosThetaMax = 1.0;
  int bMin = hBrem->GetYaxis()->FindBin(cosThetaMin);
  int bMax = hBrem->GetYaxis()->FindBin(cosThetaMax);
  int nCosThetaBins = bMax-bMin;

  vector<double> er;
  er.push_back(2.176);
  er.push_back(2.209);
  er.push_back(2.245);
  er.push_back(2.212107);

  // compute the binsIntegral and correct for narrower beam sampling if necessary
  double binsIntegral = 0.0;
  double origIntegral = 0.0;
  int oMin = hBrem->GetYaxis()->FindBin(0.995);
  for (size_t i = 0; i < er.size(); ++i){
    int erBin = hBrem->GetXaxis()->FindBin(er[i]);
    binsIntegral += hBrem->Integral(erBin, erBin, bMin, bMax);
    origIntegral += hBrem->Integral(erBin, erBin, oMin, bMax);
  }
  double integralRatio = origIntegral/binsIntegral;

  integralRatio = hBrem->ProjectionY()->Integral(oMin,bMax)/hBrem->ProjectionY()->Integral(bMin,bMax);

  double binsSampled = 1.0 * nCosThetaBins * er.size(); // Number of bins included in the sample range
  double scaleFactor = 100.0e-6 / (1.0*hBrem->GetXaxis()->GetBinWidth(1)); // adjust for the fact that sampling bins are 100 eV while histo bins are 10 keV
  double bav = binsIntegral/binsSampled/scaleFactor/integralRatio;
  cout << "binsSampled   = " << binsSampled   << endl;
  cout << "nCosThetaBins = " << nCosThetaBins << endl;
  cout << "binsIntegral  = " << binsIntegral  << endl;
  cout << "integralRatio = " << integralRatio << endl;
  cout << "bav           = " << bav           << endl;

  double weightFactor = 1.0/bav;
  double inorm = binsIntegral;

  double gammaPerFile = 1.0e7;
  int nfiles = 0;

  // I/O
  double eps = 1.0e-15;
  TH1D *peakSpectrum0 = new TH1D("peakSpectrum0","Spectrum generated in Det0;energy [MeV];weighted counts [/mC]",2500,0.0-eps,2.5-eps);
  TH1D *peakSpectrum1 = new TH1D("peakSpectrum1","Spectrum generated in Det1;energy [MeV];weighted counts [/mC]",2500,0.0-eps,2.5-eps);
  TH1D *peakSpectrum2 = new TH1D("peakSpectrum2","Spectrum generated in Det2;energy [MeV];weighted counts [/mC]",2500,0.0-eps,2.5-eps);
  peakSpectrum0->Sumw2();
  peakSpectrum1->Sumw2();
  peakSpectrum2->Sumw2();

  TH1D *ps0uw = new TH1D("ps0uw","Unweighted spectrum in Det0;energy [MeV];counts [/mc]",2500,0.0-eps,2.5-eps);
  TH1D *ps1uw = new TH1D("ps1uw","Unweighted spectrum in Det1;energy [MeV];counts [/mc]",2500,0.0-eps,2.5-eps);
  TH1D *ps2uw = new TH1D("ps2uw","Unweighted spectrum in Det2;energy [MeV];counts [/mc]",2500,0.0-eps,2.5-eps);
  ps0uw->Sumw2();
  ps1uw->Sumw2();
  ps2uw->Sumw2();

  TH1D *fplSpectrum0 = new TH1D("fplSpectrum0","Distribution of foil path lengths in Det0;path length [mm]; weighted counts [/mC]",50,0,10.0);
  fplSpectrum0->Sumw2();

  TH1D *idealEnergy0 = new TH1D("idealEnergy0","Spectrum incident on Det0;energy [MeV];weighted counts [/mC]",2500,0.0-eps,2.5-eps);
  TH1D *idealEnergy1 = new TH1D("idealEnergy1","Spectrum incident on Det1;energy [MeV];weighted counts [/mC]",2500,0.0-eps,2.5-eps);
  TH1D *idealEnergy2 = new TH1D("idealEnergy2","Spectrum incident on Det2;energy [MeV];weighted counts [/mC]",2500,0.0-eps,2.5-eps);

  TH1D *indexDepSpectrum = new TH1D("indexDepSpectrum","Detector indicies;index;weighted counts [/mC]",3,0,3);
  TH1D *indexIdeSpectrum = new TH1D("indexIdeSpectrum","Detector indicies;index;weighted counts [/mC]",3,0,3);

  TH2D *angleOut = new TH2D("angleOut","Angle at NRF creation;cos(#theta);#phi",200,-1,1,200,-TMath::Pi(),TMath::Pi());
  TH2D *pOut     = new TH2D("pOut", "Momenta at NRF creation;px [MeV];py [MeV]",200,-2.6,2.6,200,-2.6,2.6);
  angleOut->Sumw2();
  pOut->Sumw2();

  TString config = "OpenBeam_momenta_isotropic";
  TFile *outfile = new TFile(Form("peakSpectrum_%s_bruteforce.root",config.Data()),"RECREATE");

  for (int i = 0; i < 3000; ++i){

    gEnv->SetValue("TFile.Recover", 0);

    TFile *infile = new TFile(Form("/nobackup1/jvavrek/BruteForce/%s/genuinego2_%d.root",config.Data(),i));

    if (infile->IsZombie()) continue;   //file is unusable
    if (infile->TestBit(TFile::kRecovered)) continue;

    cout<<"Processing file "<< i <<", Status: "<< infile->GetFd()<<"\n";

    if (infile->GetFd() < 0)
    {
      delete infile;
      continue;
    }

    TTree *hpgeTree = (TTree*) infile->Get("hpgeTree1");
    int hpgeTreeEntries = hpgeTree->GetEntries();

    double E = 0;
    double w = 0;
    //double z = 0;
    unsigned int idx = 0;
    double fpl = 0; 
    bool qNRF = false;

    hpgeTree->SetBranchAddress("hpgeEnergy1",&E);
    hpgeTree->SetBranchAddress("weight1",&w);
    //hpgeTree->SetBranchAddress("foilZ", &z);
    hpgeTree->SetBranchAddress("ind1",&idx);
    hpgeTree->SetBranchAddress("foilPathLength1",&fpl);
    hpgeTree->SetBranchAddress("isNRF1",&qNRF);

    for (unsigned int entry = 0; entry < hpgeTreeEntries; ++entry)
    {
      hpgeTree->GetEntry(entry);
      if      (idx==0) {peakSpectrum0->Fill(E,w*weightFactor); ps0uw->Fill(E);}
      else if (idx==1) {peakSpectrum1->Fill(E,w*weightFactor); ps1uw->Fill(E);}
      else if (idx==2) {peakSpectrum2->Fill(E,w*weightFactor); ps2uw->Fill(E);}
    
      indexDepSpectrum->Fill(idx,w*weightFactor);

      if (idx==0 && fpl>0 && qNRF) fplSpectrum0->Fill(fpl, w*weightFactor);
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
        if      (idxid==0) idealEnergy0->Fill(Eid,wid*weightFactor);
        else if (idxid==1) idealEnergy1->Fill(Eid,wid*weightFactor);
        else if (idxid==2) idealEnergy2->Fill(Eid,wid*weightFactor);
        indexIdeSpectrum->Fill(idxid,w*weightFactor);
      }
    }

    TTree *foilTree = (TTree*) infile->Get("foilTree");
    if (foilTree){
      int foilTreeEntries = foilTree->GetEntries();
      double theta = 0;
      double phi = 0;
      double px = 0;
      double py = 0;
      double pz = 0;
      bool idu = false;
      foilTree->SetBranchAddress("foilpX",&px);
      foilTree->SetBranchAddress("foilpY",&py);
      foilTree->SetBranchAddress("foilpZ",&pz);
      foilTree->SetBranchAddress("foilThetaOut",&theta);
      foilTree->SetBranchAddress("foilPhiOut",&phi);
      foilTree->SetBranchAddress("foilIsDU",&idu);
      for (unsigned int ee = 0; ee < foilTreeEntries; ++ee){
        foilTree->GetEntry(ee);
        if (idu) angleOut->Fill(TMath::Cos(theta),phi);
        if (idu) pOut->Fill(px,py);
      }
    }

    delete infile;
    ++nfiles;

  } // end file loop
  cout << "nfiles = " << nfiles << endl;

  peakSpectrum0->Scale(1.0/(nfiles*gammaPerFile/inorm));
  peakSpectrum1->Scale(1.0/(nfiles*gammaPerFile/inorm));
  peakSpectrum2->Scale(1.0/(nfiles*gammaPerFile/inorm));
  peakSpectrum0->Scale(1e3); // /uC to /mC
  peakSpectrum1->Scale(1e3); // /uC to /mC
  peakSpectrum2->Scale(1e3); // /uC to /mC
  peakSpectrum0->Scale(1.0/0.932); // backscatter correction
  peakSpectrum1->Scale(1.0/0.932); // backscatter correction
  peakSpectrum2->Scale(1.0/0.932); // backscatter correction

  fplSpectrum0->Scale(1.0/(nfiles*gammaPerFile/inorm));
  fplSpectrum0->Scale(1e3);
  fplSpectrum0->Scale(1.0/0.932); // backscatter correction

  idealEnergy0->Scale(1e3/0.932/(nfiles*gammaPerFile/inorm));
  idealEnergy1->Scale(1e3/0.932/(nfiles*gammaPerFile/inorm));
  idealEnergy2->Scale(1e3/0.932/(nfiles*gammaPerFile/inorm));

  indexDepSpectrum->Scale(1e3/0.932/(nfiles*gammaPerFile/inorm));
  indexIdeSpectrum->Scale(1e3/0.932/(nfiles*gammaPerFile/inorm));

  outfile->cd();
  peakSpectrum0->Write();
  peakSpectrum1->Write();
  peakSpectrum2->Write();
  fplSpectrum0->Write();
  idealEnergy0->Write();
  idealEnergy1->Write();
  idealEnergy2->Write();
  ps0uw->Write();
  ps1uw->Write();
  ps2uw->Write();
  indexDepSpectrum->Write();
  indexIdeSpectrum->Write();
  angleOut->Write();
  pOut->Write();

  cout << "writing " << peakSpectrum0->GetName() << " et al to " << outfile->GetName() << endl;
  outfile->Write();
  delete peakSpectrum0;
  delete peakSpectrum1;
  delete peakSpectrum2;
  delete ps0uw;
  delete ps1uw;
  delete ps2uw;
  delete fplSpectrum0;
  delete idealEnergy0;
  delete idealEnergy1;
  delete idealEnergy2;
  delete indexDepSpectrum;
  delete indexIdeSpectrum;
  delete angleOut;
  delete pOut;


  //new TCanvas();
  //hRate2176->Scale(1.0/(nfiles*gammaPerFile/inorm));
  //hRate2176->GetXaxis()->SetTitle("foil depth #it{z} [mm]");
  //hRate2176->GetYaxis()->SetTitle(Form("2176 keV NRF reaction rate [/#muC/%1.3f mm bin]",hRate2176->GetBinWidth(1)));
  //hRate2176->Draw("e");
  //cout << "Integral: " << hRate2176->Integral() << " /uC" << endl;

  //TGraph *gPred = new TGraph(Form("predictedRateVsZThinGenu2176.dat"), "%lg %lg");
  //gPred->Draw("same");

  outfile->Close();
  return;
}

