void figure()
{
  TGeoManager::Import("ZKBrem1.gdml");
  gGeoManager->GetTopVolume()->Draw("ogl");

  gGeoManager->GetVolume("lead_roof_det2_L")->SetVisibility(kFALSE);
  gGeoManager->GetVolume("lead_filter_L")->SetVisibility(kFALSE);
  gGeoManager->GetVolume("downbeam_lead_det2_L")->SetVisibility(kFALSE);
  gGeoManager->GetVolume("ss_housing_L7")->SetVisibility(kFALSE);
  gGeoManager->GetVolume("ss_housing_L3")->SetVisibility(kFALSE);
  gGeoManager->GetVolume("ss_housing_L8")->SetVisibility(kFALSE);
  gGeoManager->GetVolume("ss_housing_L4")->SetVisibility(kFALSE);
  gGeoManager->GetVolume("lead_shim_L2")->SetVisibility(kFALSE);

  int tlevel = 70;

  gGeoManager->GetVolume("ss_housing_L1")->SetTransparency(tlevel);
  gGeoManager->GetVolume("extra_filter_L")->SetTransparency(tlevel);
  gGeoManager->GetVolume("downbeam_lead_det1_L")->SetTransparency(tlevel);

  gGeoManager->GetVolume("SeptThinGenuineTarget_7_L")->SetFillColor(kOrange-7);
  gGeoManager->GetVolume("SeptThinGenuineTarget_7_L")->SetLineColor(kOrange-7);
  gGeoManager->GetVolume("SeptThinGenuineTarget_0_L")->SetFillColor(kOrange-7);
  gGeoManager->GetVolume("SeptThinGenuineTarget_0_L")->SetLineColor(kOrange-7);
  gGeoManager->GetVolume("SeptThinGenuineTarget_2_L")->SetFillColor(4);
  gGeoManager->GetVolume("SeptThinGenuineTarget_2_L")->SetLineColor(4);
  gGeoManager->GetVolume("SeptThinGenuineTarget_3_L")->SetFillColor(4);
  gGeoManager->GetVolume("SeptThinGenuineTarget_3_L")->SetLineColor(4);
  gGeoManager->GetVolume("SeptThinGenuineTarget_4_L")->SetFillColor(4);
  gGeoManager->GetVolume("SeptThinGenuineTarget_4_L")->SetLineColor(4);
  gGeoManager->GetVolume("SeptThinGenuineTarget_5_L")->SetFillColor(4);
  gGeoManager->GetVolume("SeptThinGenuineTarget_5_L")->SetLineColor(4);
  gGeoManager->GetVolume("DUFoil_0_L")->SetFillColor(4);
  gGeoManager->GetVolume("DUFoil_0_L")->SetLineColor(4);

  return;
}
