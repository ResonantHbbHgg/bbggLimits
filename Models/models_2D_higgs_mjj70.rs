mgg[100,180];

CMS_hhbbgg_13TeV_mgg_bkg_slope1[0.9,-1000.0, 1000.0];
CMS_hhbbgg_13TeV_mgg_bkg_slope1_cat0[0.9,-1000.0, 1000.0];
CMS_hhbbgg_13TeV_mgg_bkg_slope1_cat1[0.9,-1000.0, 1000.0];


mgg_sig_m0_cat0[125., 122, 127];
mgg_sig_sigma_cat0[2.8, 2.1, 5.0];
mgg_sig_alpha1_cat0[1.5, 1.0, 10.0];
mgg_sig_n1_cat0[2.0, 0.1, 10.0];
mgg_sig_alpha2_cat0[1.5, 1.0, 10.0];
mgg_sig_n2_cat0[2.0, 0.1, 10.0];
mggSig_cat0 = RooDoubleCB(mgg, mgg_sig_m0_cat0, mgg_sig_sigma_cat0, mgg_sig_alpha1_cat0, mgg_sig_n1_cat0, mgg_sig_alpha2_cat0, mgg_sig_n2_cat0);

mgg_sig_m0_cat1[125., 122, 127];
mgg_sig_sigma_cat1[2.8, 2.1, 5.0];
mgg_sig_alpha1_cat1[1.5, 1.0, 10.0];
mgg_sig_n1_cat1[2.0, 0.1, 10.0];
mgg_sig_alpha2_cat1[1.5, 1.0, 10.0];
mgg_sig_n2_cat1[2.0, 0.1, 10.0];
mggSig_cat1 = RooDoubleCB(mgg, mgg_sig_m0_cat1, mgg_sig_sigma_cat1, mgg_sig_alpha1_cat1, mgg_sig_n1_cat1, mgg_sig_alpha2_cat1, mgg_sig_n2_cat1);

mgg_hig_m0_ggh_cat0[124.2, 123, 125];
mgg_hig_sigma_ggh_cat0[2.8, 2.0, 5.0];
mgg_hig_alpha1_ggh_cat0[1.5, 1.0, 10.0];
mgg_hig_n1_ggh_cat0[2.0, 0.1, 10.0];
mgg_hig_alpha2_ggh_cat0[1.5, 1.0, 10.0];
mgg_hig_n2_ggh_cat0[2.0, 0.1, 10.0];
mggHig_ggh_cat0 = RooDoubleCB(mgg, mgg_hig_m0_ggh_cat0, mgg_hig_sigma_ggh_cat0, mgg_hig_alpha1_ggh_cat0, mgg_hig_n1_ggh_cat0, mgg_hig_alpha2_ggh_cat0, mgg_hig_n2_ggh_cat0);

mgg_hig_m0_ggh_cat1[124.2, 123, 125];
mgg_hig_sigma_ggh_cat1[2.8, 2.0, 5.0];
mgg_hig_alpha1_ggh_cat1[1.5, 1.0, 10.0];
mgg_hig_n1_ggh_cat1[2.0, 0.1, 10.0];
mgg_hig_alpha2_ggh_cat1[1.5, 1.0, 10.0];
mgg_hig_n2_ggh_cat1[2.0, 0.1, 10.0];
mggHig_ggh_cat1 = RooDoubleCB(mgg, mgg_hig_m0_ggh_cat1, mgg_hig_sigma_ggh_cat1, mgg_hig_alpha1_ggh_cat1, mgg_hig_n1_ggh_cat1, mgg_hig_alpha2_ggh_cat1, mgg_hig_n2_ggh_cat1);

mgg_hig_m0_tth_cat0[124.2, 123, 125];
mgg_hig_sigma_tth_cat0[2.8, 2.0, 5.0];
mgg_hig_alpha1_tth_cat0[1.5, 1.0, 10.0];
mgg_hig_n1_tth_cat0[2.0, 0.1, 10.0];
mgg_hig_alpha2_tth_cat0[1.5, 1.0, 10.0];
mgg_hig_n2_tth_cat0[2.0, 0.1, 10.0];
mggHig_tth_cat0 = RooDoubleCB(mgg, mgg_hig_m0_tth_cat0, mgg_hig_sigma_tth_cat0, mgg_hig_alpha1_tth_cat0, mgg_hig_n1_tth_cat0, mgg_hig_alpha2_tth_cat0, mgg_hig_n2_tth_cat0);

mgg_hig_m0_tth_cat1[124.2, 123, 125];
mgg_hig_sigma_tth_cat1[2.8, 2.0, 5.0];
mgg_hig_alpha1_tth_cat1[1.5, 1.0, 10.0];
mgg_hig_n1_tth_cat1[2.0, 0.1, 10.0];
mgg_hig_alpha2_tth_cat1[1.5, 1.0, 10.0];
mgg_hig_n2_tth_cat1[2.0, 0.1, 10.0];
mggHig_tth_cat1 = RooDoubleCB(mgg, mgg_hig_m0_tth_cat1, mgg_hig_sigma_tth_cat1, mgg_hig_alpha1_tth_cat1, mgg_hig_n1_tth_cat1, mgg_hig_alpha2_tth_cat1, mgg_hig_n2_tth_cat1);

mgg_hig_m0_vh_cat0[124.2, 123, 125];
mgg_hig_sigma_vh_cat0[2.8, 2.0, 5.0];
mgg_hig_alpha1_vh_cat0[1.5, 1.0, 10.0];
mgg_hig_n1_vh_cat0[2.0, 0.1, 10.0];
mgg_hig_alpha2_vh_cat0[1.5, 1.0, 10.0];
mgg_hig_n2_vh_cat0[2.0, 0.1, 10.0];
mggHig_vh_cat0 = RooDoubleCB(mgg, mgg_hig_m0_vh_cat0, mgg_hig_sigma_vh_cat0, mgg_hig_alpha1_vh_cat0, mgg_hig_n1_vh_cat0, mgg_hig_alpha2_vh_cat0, mgg_hig_n2_vh_cat0);

mgg_hig_m0_vh_cat1[124.2, 123, 125];
mgg_hig_sigma_vh_cat1[2.8, 2.0, 5.0];
mgg_hig_alpha1_vh_cat1[1.5, 1.0, 10.0];
mgg_hig_n1_vh_cat1[2.0, 0.1, 10.0];
mgg_hig_alpha2_vh_cat1[1.5, 1.0, 10.0];
mgg_hig_n2_vh_cat1[2.0, 0.1, 10.0];
mggHig_vh_cat1 = RooDoubleCB(mgg, mgg_hig_m0_vh_cat1, mgg_hig_sigma_vh_cat1, mgg_hig_alpha1_vh_cat1, mgg_hig_n1_vh_cat1, mgg_hig_alpha2_vh_cat1, mgg_hig_n2_vh_cat1);

mgg_hig_m0_vbf_cat0[124.2, 123, 125];
mgg_hig_sigma_vbf_cat0[2.8, 2.0, 5.0];
mgg_hig_alpha1_vbf_cat0[1.5, 1.0, 10.0];
mgg_hig_n1_vbf_cat0[2.0, 0.1, 10.0];
mgg_hig_alpha2_vbf_cat0[1.5, 1.0, 10.0];
mgg_hig_n2_vbf_cat0[2.0, 0.1, 10.0];
mggHig_vbf_cat0 = RooDoubleCB(mgg, mgg_hig_m0_vbf_cat0, mgg_hig_sigma_vbf_cat0, mgg_hig_alpha1_vbf_cat0, mgg_hig_n1_vbf_cat0, mgg_hig_alpha2_vbf_cat0, mgg_hig_n2_vbf_cat0);

mgg_hig_m0_vbf_cat1[124.2, 123, 125];
mgg_hig_sigma_vbf_cat1[2.8, 2.0, 5.0];
mgg_hig_alpha1_vbf_cat1[1.5, 1.0, 10.0];
mgg_hig_n1_vbf_cat1[2.0, 0.1, 10.0];
mgg_hig_alpha2_vbf_cat1[1.5, 1.0, 10.0];
mgg_hig_n2_vbf_cat1[2.0, 0.1, 10.0];
mggHig_vbf_cat1 = RooDoubleCB(mgg, mgg_hig_m0_vbf_cat1, mgg_hig_sigma_vbf_cat1, mgg_hig_alpha1_vbf_cat1, mgg_hig_n1_vbf_cat1, mgg_hig_alpha2_vbf_cat1, mgg_hig_n2_vbf_cat1);

mgg_hig_m0_bbh_cat0[124.2, 123, 125];
mgg_hig_sigma_bbh_cat0[2.8, 2.0, 5.0];
mgg_hig_alpha1_bbh_cat0[1.5, 1.0, 10.0];
mgg_hig_n1_bbh_cat0[2.0, 0.1, 10.0];
mgg_hig_alpha2_bbh_cat0[1.5, 1.0, 10.0];
mgg_hig_n2_bbh_cat0[2.0, 0.1, 10.0];
mggHig_bbh_cat0 = RooDoubleCB(mgg, mgg_hig_m0_bbh_cat0, mgg_hig_sigma_bbh_cat0, mgg_hig_alpha1_bbh_cat0, mgg_hig_n1_bbh_cat0, mgg_hig_alpha2_bbh_cat0, mgg_hig_n2_bbh_cat0);

mgg_hig_m0_bbh_cat1[124.2, 123, 125];
mgg_hig_sigma_bbh_cat1[2.8, 2.0, 5.0];
mgg_hig_alpha1_bbh_cat1[1.5, 1.0, 10.0];
mgg_hig_n1_bbh_cat1[2.0, 0.1, 10.0];
mgg_hig_alpha2_bbh_cat1[1.5, 1.0, 10.0];
mgg_hig_n2_bbh_cat1[2.0, 0.1, 10.0];
mggHig_bbh_cat1 = RooDoubleCB(mgg, mgg_hig_m0_bbh_cat1, mgg_hig_sigma_bbh_cat1, mgg_hig_alpha1_bbh_cat1, mgg_hig_n1_bbh_cat1, mgg_hig_alpha2_bbh_cat1, mgg_hig_n2_bbh_cat1);

mjj[90,190];

mjj_sig_m0_cat0[122.0, 115, 130];
mjj_sig_sigma_cat0[20.0, 10.0, 60.0];
mjj_sig_alpha1_cat0[0.3, 0.2, 2.0];
mjj_sig_n1_cat0[1.0, 0.5, 10.0];
mjj_sig_alpha2_cat0[3, 0.5, 2.0];
mjj_sig_n2_cat0[5.0, 1.0, 100.0];
mjj_sig_gm0_cat0[122.0, 110, 130];
mjj_sig_gsigma_cat0[25.0, 10.0, 60.0];
mjj_sig_frac_cat0[0.1, 0, 1];

mjjGaussSig_cat0 = Gaussian(mjj, mjj_sig_gm0_cat0, mjj_sig_gsigma_cat0);
mjjCBSig_cat0 = CBShape(mjj, mjj_sig_m0_cat0, mjj_sig_sigma_cat0, mjj_sig_alpha1_cat0, mjj_sig_n1_cat0);
mjjSig_cat0 = AddPdf(mjjGaussSig_cat0, mjjCBSig_cat0, mjj_sig_frac_cat0);

mjj_sig_m0_cat1[122.0, 115, 130];
mjj_sig_sigma_cat1[20.0, 10.0, 60.0];
mjj_sig_alpha1_cat1[0.3, 0.2, 2.0];
mjj_sig_n1_cat1[1.0, 0.5, 10.0];
mjj_sig_alpha2_cat1[3, 0.5, 2.0];
mjj_sig_n2_cat1[5.0, 1.0, 100.0];
mjj_sig_gm0_cat1[122.0, 110, 130];
mjj_sig_gsigma_cat1[25.0, 10.0, 60.0];
mjj_sig_frac_cat1[0.1, 0, 1];

mjjGaussSig_cat1 = Gaussian(mjj, mjj_sig_gm0_cat1, mjj_sig_gsigma_cat1);
mjjCBSig_cat1 = CBShape(mjj, mjj_sig_m0_cat1, mjj_sig_sigma_cat1, mjj_sig_alpha1_cat1, mjj_sig_n1_cat1);
mjjSig_cat1 = AddPdf(mjjGaussSig_cat1, mjjCBSig_cat1, mjj_sig_frac_cat1);


CMS_hhbbgg_13TeV_mjj_bkg_slope1[0.9,-1000.0, 1000.0];
CMS_hhbbgg_13TeV_mjj_bkg_slope1_cat0[0.9,-1000.0, 1000.];
CMS_hhbbgg_13TeV_mjj_bkg_slope1_cat1[0.9,-1000.0, 1000.];

mjj_hig_slope1_ggh_cat0[0.1, 0, 10];
mjj_hig_slope2_ggh_cat0[0.1, 0, 10];
mjj_hig_slope3_ggh_cat0[0.1, 0, 10];

mjj_hig_slope1_ggh_cat1[0.1, -1000, 1000];
mjj_hig_slope2_ggh_cat1[0.1, -1000, 1000];
mjj_hig_slope3_ggh_cat1[0.1, -1000, 1000];

mjj_hig_slope1_vbf_cat0[0.1, 0, 10];
mjj_hig_slope2_vbf_cat0[0.1, 0, 10];
mjj_hig_slope3_vbf_cat0[0.1, 0, 10];

mjj_hig_slope1_vbf_cat1[0.1, -1000, 1000];
mjj_hig_slope2_vbf_cat1[0.1, -1000, 1000];
mjj_hig_slope3_vbf_cat1[0.1, -1000, 1000];

mjj_hig_m0_tth_cat0[100, 70, 190];
mjj_hig_sigma_tth_cat0[50, 10, 100];
mjj_hig_alpha1_tth_cat0[1.0, 0.01, 10];
mjj_hig_n1_tth_cat0[1, 0.01, 10];
mjj_hig_alpha2_tth_cat0[1.0, 0.01, 10];
mjj_hig_n2_tth_cat0[1, 0.01, 10];
mjjHig_tth_cat0 = RooDoubleCB(mjj, mjj_hig_m0_tth_cat0, mjj_hig_sigma_tth_cat0, mjj_hig_alpha1_tth_cat0, mjj_hig_n1_tth_cat0, mjj_hig_alpha2_tth_cat0, mjj_hig_n2_tth_cat0);

mjj_hig_m0_tth_cat1[100, 70, 190];
mjj_hig_sigma_tth_cat1[50, 10, 100];
mjj_hig_alpha1_tth_cat1[1.0, 0.01, 10];
mjj_hig_n1_tth_cat1[1, 0.01, 10];
mjj_hig_alpha2_tth_cat1[1.0, 0.01, 10];
mjj_hig_n2_tth_cat1[1, 0.01, 10];
mjjHig_tth_cat1 = RooDoubleCB(mjj, mjj_hig_m0_tth_cat1, mjj_hig_sigma_tth_cat1, mjj_hig_alpha1_tth_cat1, mjj_hig_n1_tth_cat1, mjj_hig_alpha2_tth_cat1, mjj_hig_n2_tth_cat1);

mjj_hig_m0_vh_cat0[100, 70, 190];
mjj_hig_sigma_vh_cat0[50, 10, 100];
mjj_hig_alpha1_vh_cat0[1.0, 0.01, 10];
mjj_hig_n1_vh_cat0[1, 0.01, 10];
mjj_hig_alpha2_vh_cat0[1.0, 0.01, 10];
mjj_hig_n2_vh_cat0[1, 0.01, 10];
mjjHig_vh_cat0 = RooDoubleCB(mjj, mjj_hig_m0_vh_cat0, mjj_hig_sigma_vh_cat0, mjj_hig_alpha1_vh_cat0, mjj_hig_n1_vh_cat0, mjj_hig_alpha2_vh_cat0, mjj_hig_n2_vh_cat0);

mjj_hig_m0_vh_cat1[100, 70, 190];
mjj_hig_sigma_vh_cat1[50, 10, 100];
mjj_hig_alpha1_vh_cat1[1.0, 0.01, 10];
mjj_hig_n1_vh_cat1[1, 0.01, 10];
mjj_hig_alpha2_vh_cat1[1.0, 0.01, 10];
mjj_hig_n2_vh_cat1[1, 0.01, 10];
mjjHig_vh_cat1 = RooDoubleCB(mjj, mjj_hig_m0_vh_cat1, mjj_hig_sigma_vh_cat1, mjj_hig_alpha1_vh_cat1, mjj_hig_n1_vh_cat1, mjj_hig_alpha2_vh_cat1, mjj_hig_n2_vh_cat1);

mjj_hig_m0_bbh_cat0[100, 10, 180];
mjj_hig_sigma_bbh_cat0[50, 1.0, 100];
mjj_hig_alpha1_bbh_cat0[1.0, 0.01, 10];
mjj_hig_n1_bbh_cat0[1, 0.01, 10];
mjj_hig_alpha2_bbh_cat0[1.0, 0.01, 10];
mjj_hig_n2_bbh_cat0[1, 0.01, 10];
mjjHig_bbh_cat0 = RooDoubleCB(mjj, mjj_hig_m0_bbh_cat0, mjj_hig_sigma_bbh_cat0, mjj_hig_alpha1_bbh_cat0, mjj_hig_n1_bbh_cat0, mjj_hig_alpha2_bbh_cat0, mjj_hig_n2_bbh_cat0);

mjj_hig_m0_bbh_cat1[100, 10, 180];
mjj_hig_sigma_bbh_cat1[50, 1.0, 100];
mjj_hig_alpha1_bbh_cat1[1.0, 0.01, 10];
mjj_hig_n1_bbh_cat1[1, 0.01, 10];
mjj_hig_alpha2_bbh_cat1[1.0, 0.01, 10];
mjj_hig_n2_bbh_cat1[1, 0.01, 10];
mjjHig_bbh_cat1 = RooDoubleCB(mjj, mjj_hig_m0_bbh_cat1, mjj_hig_sigma_bbh_cat1, mjj_hig_alpha1_bbh_cat1, mjj_hig_n1_bbh_cat1, mjj_hig_alpha2_bbh_cat1, mjj_hig_n2_bbh_cat1);

mtot[200,1600];

ttHTagger[-1,1];
