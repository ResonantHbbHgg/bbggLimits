void Skim(string sample){

  cout << sample.c_str() << endl;

  int frontier = 490;

  string predir =  string("/eos/cms/store/group/phys_higgs/resonant_HH/RunII/LimitTrees/20180908_MILANO/PhotonLoose_mbb_80-200_MTD_v5/");

  string dir = predir + string("Combine_v3/LimitTrees_HighMass");
  string new_predir =  Form("SplitHigh_v3_Mass_%d/", frontier);

  string mkdir = string("mkdir ") + predir + new_predir;
  gSystem->Exec(mkdir.c_str());
  mkdir = string("mkdir ") + predir + new_predir + "LimitTrees_HighMass";
  gSystem->Exec(mkdir.c_str());
  mkdir = string("mkdir ") + predir + new_predir + "LimitTrees_LowMass";
  gSystem->Exec(mkdir.c_str());
		       

  string filename = dir + "/" + sample;
  string fileout_high = dir + Form("/High_%d", frontier) + sample;
  string fileout_low = dir + Form("/Low_%d", frontier)  + sample;
  string newfileout_high = predir + new_predir + "/LimitTrees_HighMass/" + sample;
  string newfileout_low = predir + new_predir + "/LimitTrees_LowMass/"  + sample;

  TFile *_file0 = TFile::Open(filename.c_str());

  
  TTree *oldtree = (TTree*) _file0->Get("TCVARS;1");    

  TFile *newfile_high = new TFile(fileout_high.c_str(), "recreate");  
  TTree *newtree_high = oldtree->CloneTree(0); 

  newtree_high->SetMakeClass(1);
  oldtree->CopyAddresses(newtree_high);

  float mtot   = 0;
  oldtree->SetBranchAddress("mtot",&mtot);

  for (Long64_t iEntry=0; iEntry<oldtree->GetEntriesFast(); iEntry++)
    { // Main loop
      oldtree->GetEntry(iEntry); // get ith entry

      if(mtot > frontier)
	  newtree_high->Fill();
    }

  newfile_high->Write();
  newfile_high->Close();
  

  TFile *newfile_low = new TFile(fileout_low.c_str(), "recreate");  
  TTree *newtree_low = oldtree->CloneTree(0); 

  newtree_low->SetMakeClass(1);
  oldtree->CopyAddresses(newtree_low);


  for (Long64_t iEntry=0; iEntry<oldtree->GetEntriesFast(); iEntry++)
    { // Main loop
      oldtree->GetEntry(iEntry); // get ith entry

      if(mtot < frontier)
	  newtree_low->Fill();
    }

  newfile_low->Write();
  newfile_low->Close();

  string cp = "cp " + fileout_high + " " + newfileout_high;
  gSystem->Exec(cp.c_str());
  cp = "cp " + fileout_low + " " + newfileout_low;
  gSystem->Exec(cp.c_str());



}

void SkimTree(){

  Skim("LT_DoubleEG.root");
  Skim("LT_output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root");
  Skim("LT_output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph_0.root");
  Skim("LT_output_VBFHToGG_M-125_13TeV_powheg_pythia8.root");
  Skim("LT_output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8.root");
  Skim("LT_output_bbHToGG_M-125_13TeV_amcatnlo.root");
  Skim("LT_output_ttHToGG_M125_13TeV_powheg_pythia8_v2.root");


  


}
