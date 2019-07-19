#include "histSaver.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "observable.h"
histSaver::histSaver(TString _outputfilename) {
  outputfile = new TFile (_outputfilename + ".root", "recreate");
  outputfilename = _outputfilename;
  nvar = 0;
  inputfilename = "hists";
  nregion = 0;
  blinding = 0;
  fweight = NULL;
  dweight = NULL;
  weight_type = 0;
  doROC = 0;
  inputfile = 0;
  lumi = "#it{#sqrt{s}} = 13TeV, 80 fb^{-1}";
  analysis = "FCNC tqH H#rightarrow tautau";
  workflow = "work in progress";
  fromntuple = 1;
  histcount = 0;
  this_region = "nominal";
  read_path = "./" ;
  debug = 1;
  sensitivevariable = "";
  for(Int_t i=0; i<50; i++) {
    nbin[i] = 1; xlo[i] = 0; xhi[i] = 1; var1[i] = 0; var2[i] = 0; MeVtoGeV[i] = 0; var3[i] = 0;
  }
}

histSaver::~histSaver() {
  printf("histSaver::~histSaver()\n");
  for(auto& samp : plot_lib){
    for(auto &reg: samp.second) {
      for (int i = 0; i < nvar; ++i){
        TH1D *target = reg.second[i];
        cout<<"\rdeleting histogram: "<<target->GetName()<<std::flush;
          deletepointer(target);
        cout<<"\rdone deleting histogram"<<std::flush;
      }
    }
  }
  deletepointer(inputfile);
  deletepointer(outputfile);
  printf("histSaver::~histSaver() destructed\n");
}

void histSaver::printyield(TString region){
  printf("Print Yeild: %s\n", region.Data());
  double er;
  for(auto iter: plot_lib){
    TH1D* target = grabhist(iter.first,region,0);
    if(target){
      printf("%s: %4.3f \\pm %4.3f\n", iter.first.Data(), target->IntegralAndError(1,target->GetNbinsX(), er), er);
    }else{
      printf("Warning: histogram not found: %s, %s, %s\n", iter.first.Data(), region.Data(), name[0].Data());
    }
  }
}

TH1D* histSaver::grabhist(TString sample, TString region, int ivar){
  if(plot_lib.find(sample) == plot_lib.end()){
    if(debug) {
      //show();
      printf("histSaver:grabhist  Warning: sample %s not found\n", sample.Data());
    }
    return 0;
  }

  if(plot_lib[sample].find(region) == plot_lib[sample].end()){
    if(debug) printf("histSaver:grabhist  Warning: region %s for sample %s not found\n", region.Data(), sample.Data());
    if(debug) show();
    return 0;
  }
  if(!plot_lib[sample][region][ivar]) if(debug) printf("histSaver:grabhist  WARNING: empty histogram%s\n", name[ivar].Data());
  return plot_lib[sample][region][ivar];
}

TH1D* histSaver::grabhist(TString sample, TString region, TString varname){
  int ivar = -1;
  for (int i = 0; i < nvar; ++i)
  {
    if(varname == name[i]){
      ivar = i;
      break;
    }
  }
  if(plot_lib.find(sample) == plot_lib.end()){
    if(debug) {
      //show();
      printf("histSaver:grabhist  Warning: sample %s not found\n", sample.Data());
    }
    return 0;
  }

  if(plot_lib[sample].find(region) == plot_lib[sample].end()){
    if(debug) printf("histSaver:grabhist  Warning: region %s for sample %s not found\n", region.Data(), sample.Data());
    if(debug) show();
    return 0;
  }
  if(!plot_lib[sample][region][ivar]) if(debug) printf("histSaver:grabhist  WARNING: empty histogram%s\n", name[ivar].Data());
  return plot_lib[sample][region][ivar];
}

TH1D* histSaver::grabbkghist(TString region, int ivar){
  TH1D *hist = 0;
  for(auto iter: stackorder){
    if(find(overlaysamples.begin(),overlaysamples.end(),iter)== overlaysamples.end() && iter != "data"){
      TH1D *target = grabhist(iter,region,ivar);
      if(target){
        if(hist == 0) hist = (TH1D*)target->Clone();
        else hist->Add(target);
      }
    }
  }
  return hist;
}

TH1D* histSaver::grabsighist(TString region, int ivar, TString signal){
  TH1D *hist = 0;
  if(overlaysamples.size()==0) {
    printf("signal not defined, please call overlaysample(\"signal\")\n");
    exit(0);
  }
  return grabhist(signal==""?overlaysamples[0]:signal, region, ivar);
}

TH1D* histSaver::grabdatahist(TString region, int ivar){
  return grabhist("data",region,ivar);
}

void histSaver::add(Int_t nbin_, const Double_t* xbins_, const char* titleX_, const char* name_, Int_t* var_, const char* unit_) {
  if(nvar>=0 && nvar<20) {
    nbin[nvar] = nbin_;
    for(int i=0; i<=nbin_; i++) {
      xbins[nvar][i] = xbins_[i];
    }
    titleX[nvar] = titleX_; var2[nvar] = var_;
    name[nvar] = name_;
    unit[nvar] = unit_;
    ifRebin[nvar] = 1;
    nvar++;
  }
}


void histSaver::add(const char* titleX_, const char* name_, const char* unit_, int _rebin) {
  fromntuple = 0;
  if(nvar>=0 && nvar<50) {
    titleX[nvar] = titleX_;
    name[nvar] = name_ ;
    unit[nvar] = unit_;
    rebin[nvar] = _rebin;
    nvar++;
  }
}

void histSaver::add(Int_t nbin_, const Double_t* xbins_, const char* titleX_, const char* name_, Float_t* var_, Bool_t MeVtoGeV_, const char* unit_) {
  if(nvar>=0 && nvar<20) {
    nbin[nvar] = nbin_;
    for(int i=0; i<=nbin_; i++) {
      xbins[nvar][i] = xbins_[i];
    }
    titleX[nvar] = titleX_; var1[nvar] = var_; MeVtoGeV[nvar] = MeVtoGeV_;
    name[nvar] = name_;
    unit[nvar] = unit_;
    ifRebin[nvar] = 1;
    nvar++;
  }
}

Float_t histSaver::getVal(Int_t i) {
  Float_t tmp = -999999;
  if(i>=0  && i<nvar) {
    if(var1[i])      { if(debug) printf("fill var1\n"); tmp = MeVtoGeV[i] ? *var1[i]/1000 : *var1[i]; }
    else if(var3[i]) { if(debug) printf("fill var3\n"); tmp = MeVtoGeV[i] ? *var3[i]/1000 : *var3[i]; }
    else if(var2[i]) { if(debug) printf("fill var2\n"); tmp = *var2[i]; }
    else printf("error: fill variable failed. no var available\n");
    if(debug == 1) printf("fill value: %f\n", tmp);
  }
  if (!ifRebin[i]){
    if(tmp >= xhi[i]) tmp = xhi[i]*0.999999;
    if(tmp < xlo[i]) tmp = xlo[i];
  }else{
    if(tmp >= xbins[i][nbin[i]]) tmp = xbins[i][nbin[i]]*0.999999;
    if(tmp < xbins[i][nbin[i]]) tmp = xbins[i][nbin[i]];
  }
  return tmp;
}

void histSaver::show(){
  for(auto& iter:plot_lib ){
    printf("histSaver::show()\tsample: %s\n", iter.first.Data());
  }
  for(auto const& region: regions) {
    printf("histSaver::show()\tregion: %s\n", region.Data());
  }
  for (int i = 0; i < nvar; ++i)
  {
    if(var2[i]) printf("histSaver::show()\t%s = %d\n", name[i].Data(), *var2[i]);
    else if(var1[i]) printf("histSaver::show()\t%s = %f\n", name[i].Data(), MeVtoGeV[i] ? *var1[i]/1000 : *var1[i]);
    else if(var3[i]) printf("histSaver::show()\t%s = %f\n", name[i].Data(), MeVtoGeV[i] ? *var3[i]/1000 : *var3[i]);
  }
}

float histSaver::binwidth(int i){
  return (xhi[i]-xlo[i])/nbin[i];
}

void histSaver::merge_regions(TString inputregion1, TString inputregion2, TString outputregion){
  if(debug) printf("histSaver::merge_regions\t %s and %s into %s\n",inputregion1.Data(),inputregion2.Data(),outputregion.Data());
  bool exist = 0;
  bool input1exist = 1;
  bool input2exist = 1;

  for(auto& iter:plot_lib ){
    if(iter.second.find(inputregion1) == iter.second.end()){
      if(debug) printf("histSaver::merge_regions\t inputregion1: %s not found for sample %s\n",inputregion1.Data(), iter.first.Data());
      input1exist = 0;
    }
    if(iter.second.find(inputregion2) == iter.second.end()){
      if(debug) printf("histSaver::merge_regions\t inputregion2: %s not found for sample %s\n",inputregion2.Data(), iter.first.Data());
      input2exist = 0;
    }
    if(input1exist == 0 && input2exist == 0) continue;
    bool outputexist = 0;
    if(iter.second.find(outputregion) != iter.second.end()){
      if(debug) printf("histSaver::merge_regions\t outputregion %s exist, overwrite it\n",outputregion.Data());
      exist = 1;
      for (int i = 0; i < nvar; ++i) deletepointer(iter.second[outputregion][i]);
      iter.second[outputregion].clear();
    }
    for (int i = 0; i < nvar; ++i)
    {
      if(input1exist == 1) iter.second[outputregion].push_back((TH1D*)iter.second[inputregion1][i]->Clone(iter.first+"_"+outputregion+"_"+name[i]));
      else iter.second[outputregion].push_back((TH1D*)iter.second[inputregion2][i]->Clone(iter.first+"_"+outputregion+"_"+name[i]));
      if(input1exist == 1 && input2exist == 1) {
        iter.second[outputregion][i]->Add(iter.second[inputregion2][i]);
        if(debug)
          printf("add %s to %s as %s\n", iter.second[inputregion2][i]->GetName(),iter.second[inputregion1][i]->GetName(),(iter.first+"_"+outputregion+"_"+name[i]).Data());
      }
    }
  }
  //for(auto& iter:plot_lib ){
  //  if(iter.second.find(outputregion) == iter.second.end()){
  //    printf("histSaver::merge_regions merge regions failed\n");
  //    exit(1);
  //  }
  //}
  if(!exist) regions.push_back(outputregion);
}

void histSaver::init_sample(TString samplename, TString histname, TString sampleTitle, enum EColor color){

  outputfile->cd();
  current_sample = samplename;

  if(plot_lib.find(samplename) != plot_lib.end()) return;
  
  if(debug) printf("add new sample: %s\n", samplename.Data());

  for(auto const& region: regions) {
    for (int i = 0; i < nvar; ++i){
      plot_lib[samplename][region].push_back(new TH1D(histname + "_" + region + "_" + name[i],sampleTitle,nbin[i],xlo[i],xhi[i]));
      if (samplename != "data")
      {
        plot_lib[samplename][region][i]->Sumw2();
        plot_lib[samplename][region][i]->SetFillColor(color);
        plot_lib[samplename][region][i]->SetLineWidth(1);
        plot_lib[samplename][region][i]->SetLineColor(kBlack);
        plot_lib[samplename][region][i]->SetMarkerSize(0);
      }
    }
    if(debug == 1) printf("plot_lib[%s][%s]\n", samplename.Data(), region.Data());
  }
  if (samplename == "data") dataref = 1;

  if(debug) printf("finished initializing %s\n", samplename.Data() );
}

void histSaver::read_sample(TString samplename, TString histname, TString sampleTitle, enum EColor color, double norm, TFile *_inputfile){

  if(!inputfile) inputfile = new TFile(inputfilename + ".root", "read");
  TFile *readfromfile;

  if(_inputfile) readfromfile = _inputfile;
  else readfromfile = inputfile;

  if (samplename == "data") dataref = 1;
  for(auto const& region: regions) {
    if (debug == 1)
    {
      printf("read sample %s from %s region\n", samplename.Data(), region.Data());
    }
    ++histcount;
    if(!(TH1D*)(readfromfile->Get(histname+"_"+region+"_"+name[0]))) {
      if(debug) printf("histogram name not found: %s\n", (histname+"_"+region+"_"+name[0]).Data());
      continue;
    }
    if (plot_lib[samplename].find(region) != plot_lib[samplename].end())
    {
      for (int i = 0; i < nvar; ++i)
      {
        if(!(TH1D*)(readfromfile->Get(histname+"_"+region+"_"+name[i]))) {
          if(debug) printf("histogram name not found: %s\n", (histname+"_"+region+"_"+name[i]).Data());
          printf("plot_lib[%s][%s][%d]\n", samplename.Data(), region.Data(), i);
          show();
          exit(1);
        }
        double tmp = ((TH1D*)readfromfile->Get(histname+"_"+region+"_"+name[i]))->Integral();
        if(tmp!=tmp){
          printf("Warning: %s->Integral() is nan, skip\n", (histname+"_"+region+"_"+name[i]).Data());
          continue;
        }

        plot_lib[samplename][region][i]->Add((TH1D*)readfromfile->Get(histname+"_"+region+"_"+name[i]),norm);
      }
    }else{
      for (int i = 0; i < nvar; ++i){
        if(!(TH1D*)(readfromfile->Get(histname+"_"+region+"_"+name[i]))) {
          if(debug) printf("histogram name not found: %s\n", (histname+"_"+region+"_"+name[i]).Data());
          printf("plot_lib[%s][%s][%d]\n", samplename.Data(), region.Data(), i);
          show();
          exit(1);
        }
        double tmp = ((TH1D*)readfromfile->Get(histname+"_"+region+"_"+name[i])->Clone(histname+"_"+region+"_"+name[i]))->Integral();
        if(tmp!=tmp){
          printf("Warning: New hist: %s->Integral() is nan, continue\n", (histname+"_"+region+"_"+name[i]).Data());
          continue;
        }

        plot_lib[samplename][region].push_back((TH1D*)(readfromfile->Get(histname+"_"+region+"_"+name[i])->Clone(histname+"_"+region+"_"+name[i])));

        plot_lib[samplename][region][i]->SetName(samplename+"_"+region+"_"+name[i]);
        plot_lib[samplename][region][i]->Scale(norm);
        plot_lib[samplename][region][i]->SetTitle(sampleTitle);
        plot_lib[samplename][region][i]->SetFillColorAlpha(color,1);
        plot_lib[samplename][region][i]->SetLineWidth(1);
        plot_lib[samplename][region][i]->SetLineColor(kBlack);
        plot_lib[samplename][region][i]->SetMarkerSize(0);
        if(histcount == 1){
          nbin[i] = plot_lib[samplename][region][i]->GetNbinsX();
          xlo[i] = plot_lib[samplename][region][i]->GetXaxis()->GetXmin();
          xhi[i] = plot_lib[samplename][region][i]->GetXaxis()->GetXmax();
        }
      }
    }
  }
}

void histSaver::add_region(TString region){
  regions.push_back(region);
  nregion += 1;
}

void histSaver::fill_hist(TString sample, TString region){
  if (weight_type == 0)
  {
    printf("ERROR: weight not set\n");
  }
  for (int i = 0; i < nvar; ++i){
    double fillval = getVal(i);
    if(fillval!=fillval) {
      printf("Warning: fill val is nan: \n");
      printf("plot_lib[%s][%s][%d]->Fill(%f,%f)\n", sample.Data(), region.Data(), i, fillval, weight_type == 1? *fweight : *dweight);
    }
    if(debug == 1) printf("plot_lib[%s][%s][%d]->Fill(%f,%f)\n", sample.Data(), region.Data(), i, fillval, weight_type == 1? *fweight : *dweight);
    grabhist(sample,region,i)->Fill(fillval,weight_type == 1? *fweight : *dweight);
  }
}

void histSaver::fill_hist(TString region){
  fill_hist(current_sample,region);
}

void histSaver::fill_hist(){
  fill_hist("nominal");
}

void histSaver::write(){
  if(!outputfile) {
    printf("histSaver::write Error: outputfile pointer is empty\n");
    exit(1);
  }
  printf("histSaver::write() Write to file: %s\n", outputfile->GetName());
  for(auto const& region: regions) {
    for(auto& iter : plot_lib){
      double tmp = grabhist(iter.first,region,0)->Integral();
      if(tmp == 0) continue;
      if(tmp != tmp) {
        printf("Warning: hist integral is nan, skip writing for %s\n", grabhist(iter.first,region,0)->GetName());
        continue;
      }
      for (int i = 0; i < nvar; ++i){
        outputfile->cd();
        //if(grabhist(iter.first,region,i)->Integral() == 0) {
        //  printf("Warning: histogram is empty: %s, %s, %d\n", iter.first.Data(),region.Data(),i);
        //}
        grabhist(iter.first,region,i)->Write("",TObject::kWriteDelete);
      }
    }
  }
  printf("histSaver::write() Written\n");
}

void histSaver::write_trexinput(TString NPname, TString writeoption){
  TString trexdir = "trexinputs";
  gSystem->mkdir("trexinputs");
  for (int i = 0; i < nvar; ++i){
    gSystem->mkdir(trexdir + "/" + name[i]);
    for(auto const& region: regions) {
      bool muted = 0;
      for (auto const& mutedregion: mutedregions)
      {
        if(region.Contains(mutedregion))
          muted = 1;
      }
      if(muted) continue;

      gSystem->mkdir(trexdir + "/" + name[i] + "/" + region);
      for(auto& iter : plot_lib){
        if(NPname == "NOMINAL" && iter.first.Contains("data")) continue;
        TString filename = trexdir + "/" + name[i] + "/" + region + "/" + iter.first + ".root";
        TFile outputfile(filename, writeoption);
        if(debug) printf("Writing to file: %s, histoname: %s\n", filename.Data(), NPname.Data());
        TH1D *target = grabhist(iter.first,region,i);
        if(target) target->Write(NPname,TObject::kWriteDelete);
        else if(debug) printf("Warning: histogram plot_lib[%s][%s][%d] not found\n", iter.first.Data(),region.Data(),i);
        outputfile.Close();
      }
    }
  }
}
void histSaver::clearhist(){
  if(debug) printf("histSaver::clearhist()\n");
  for(auto const& region: regions) {
    for (int i = 0; i < nvar; ++i){
      for(auto& iter : plot_lib){
        if(iter.second[region][i]){
          iter.second[region][i]->Reset();
        }else{
          printf("histSaver::Reset() Error: histogram not found: sample: %s, variable: %s, region: %s\n",iter.first.Data(), name[i].Data(),region.Data());
        }
      }
    }
  }
}

void histSaver::overlay(TString _overlaysample){
  overlaysamples.push_back(_overlaysample);
}

double histSaver::templatesample(TString fromregion,string formula,TString toregion,TString newsamplename,TString newsampletitle,enum EColor color, bool scaletogap, double SF){
  istringstream iss(formula);
  vector<string> tokens{istream_iterator<string>{iss},
    istream_iterator<string>{}};
  if(tokens.size()%2) printf("Error: Wrong formula format: %s\nShould be like: 1 data -1 real -1 zll ...", formula.c_str());
  vector<TH1D*> newvec;
  observable scaleto(0,0);
  for (int ivar = 0; ivar < nvar; ++ivar)
  {
    newvec.push_back((TH1D*)grabhist(tokens[1],fromregion,ivar)->Clone(newsamplename+"_"+toregion+name[ivar]));
    newvec[ivar]->Reset();
    newvec[ivar]->SetNameTitle(newsamplename,newsampletitle);
    newvec[ivar]->SetFillColor(color);
  }
  for (int i = 0; i < tokens.size()/2; ++i)
  {
    int icompon = 2*i;
    float numb = 0;
    try{
      numb = stof(tokens[icompon]);
    } 
    catch(const std::invalid_argument& e){
      printf("Error: Wrong formula format: %s\nShould be like: 1 data -1 real -1 zll ...", formula.c_str());
      exit(1);
    }
    if(grabhist(tokens[icompon+1],toregion,0)){
      if(scaletogap) {
        double error = 0;
        observable tmp(grabhist(tokens[icompon+1],toregion,0)->Integral(),gethisterror(grabhist(tokens[icompon+1],toregion,0)));
        scaleto += tmp*numb;
      }
      for (int ivar = 0; ivar < nvar; ++ivar)
      {
        newvec[ivar]->Add(grabhist(tokens[icompon+1],fromregion,ivar),numb);
      }
    }
  }
  observable scalefactor;
  if(scaletogap) {
    observable scalefrom(newvec[0]->Integral(),gethisterror(newvec[0]));
    scalefactor = scaleto/scalefrom;
    printf("scale from %s: %f +/- %f\nto %s: %f +/- %f\nratio: %f +/- %f\n\n",
      fromregion.Data(), scalefrom.nominal, scalefrom.error,
      toregion.Data(), scaleto.nominal, scaleto.error,
      scalefactor.nominal, scalefactor.error);
    for(auto & hists : newvec){
      hists->Scale(scalefactor.nominal);
    }
  }else{
    for(auto & hists : newvec){
      hists->Scale(SF);
    }
  }
  for(int ivar = 0; ivar < nvar; ivar++){
    plot_lib[newsamplename][toregion] = newvec;
  }
  return scalefactor.nominal;
}

double histSaver::gethisterror(TH1* hist){
  double error = 0;
  for (int i = 0; i < hist->GetNbinsX(); ++i)
  {
    error+=pow(hist->GetBinError(i+1),2);
  }
  return sqrt(error);
}

void histSaver::muteregion(TString region){
  mutedregions.push_back(region);
}

void histSaver::unmuteregion(TString region){
  std::vector<TString>::iterator it = find(mutedregions.begin(), mutedregions.end(), region);
  if(it != mutedregions.end()) mutedregions.erase(it);
  else printf("histSaver::unmuteregion WARNING: region %s is not in the mute list\n",region.Data());
}

void histSaver::SetLumiAnaWorkflow(TString _lumi, TString _analysis, TString _workflow){
  lumi = _lumi;
  analysis = _analysis;
  workflow = _workflow;
}
void histSaver::plot_stack(TString outputdir){
  SetAtlasStyle();
  TGaxis::SetMaxDigits(3);
  gSystem->mkdir(outputdir);
  TCanvas cv("cv","cv",600,600);
  TGraph* ROC;
  TH1D *ROC_sig = 0;
  TH1D *ROC_bkg = 0;
  TH1D *datahistblinded = 0;
  vector<TH1D*> buffer;
  if(doROC){
    ROC = new TGraph();
    ROC -> SetName("ROC");
    ROC -> SetTitle("ROC");
  }
  for(auto const& region: regions) {
    bool muted = 0;
    for (auto const& mutedregion: mutedregions)
    {
      if(region.Contains(mutedregion))
        muted = 1;
    }
    if(muted) continue;
    gSystem->mkdir(outputdir + "/" + region);
    for (int i = 0; i < nvar; ++i){
      cv.SaveAs(outputdir + "/" + region + "/" + name[i] + ".pdf[");
      TPad *padlow = new TPad("lowpad","lowpad",0,0,1,0.3);
      TPad *padhi  = new TPad("hipad","hipad",0,0.3,1,1);
      TH1D hmc("hmc","hmc",nbin[i]/rebin[i],xlo[i],xhi[i]);
      TH1D hmcR("hmcR","hmcR",nbin[i]/rebin[i],xlo[i],xhi[i]);
      TH1D hdataR("hdataR","hdataR",nbin[i]/rebin[i],xlo[i],xhi[i]);
      cv.cd();
      padhi->Draw();
//===============================upper pad bkg and unblinded data===============================
      padhi->SetBottomMargin(0.017);
      padhi->SetRightMargin(0.08);
      padhi->SetLeftMargin(0.12);
      padhi->cd();
      hmc.Sumw2();
      THStack *hsk = new THStack(name[i].Data(),name[i].Data());
      TLegend* lg1 = 0;
      lg1 = new TLegend(0.45,0.7,0.90,0.9,"");
      lg1->SetNColumns(2);
      TH1D *histoverlay;
      if(debug) printf("set hists\n");
      for(auto& iter:stackorder ){
        if(iter == "data") continue;
        if(debug) {
          printf("plot_lib[%s][%s][%d]\n", iter.Data(), region.Data(), i);
        }
        if(grabhist(iter,region,i)) buffer.push_back((TH1D*)grabhist(iter,region,i)->Clone());
        else continue;
        if(doROC && sensitivevariable == name[i])
        {
          if(!ROC_bkg) ROC_bkg = (TH1D*) buffer.back()->Clone();
          else ROC_bkg->Add(buffer.back());
        }
        if(rebin[i] != 1) buffer.back()->Rebin(rebin[i]);
        hsk->Add(buffer.back());
        hmc.Add(buffer.back());
        lg1->AddEntry(buffer.back(),buffer.back()->GetTitle(),"F");
      }
      if(!hsk->GetMaximum()){
        printf("ERROR: stack has no entry, continue\n");
        continue;
      }
      double histmax = hmc.GetMaximum() + hmc.GetBinError(hmc.GetMaximumBin());

      TH1D * datahist;
      if(debug) printf("set data\n");
      if (dataref) {
        if(grabhist("data",region,i)) datahist = (TH1D*)grabhist("data",region,i)->Clone();
        datahist->Rebin(rebin[i]);
        if(datahist->Integral() == 0) printf("Warning: data hist is empty\n");
        lg1->AddEntry(datahist,"data","LP");
        datahist->SetMarkerStyle(20);
        datahist->SetMarkerSize(0.4);
        datahist->SetMinimum(0);
        histmax = max(histmax, datahist->GetMaximum() + datahist->GetBinError(datahist->GetMaximumBin()));
      }else{
        hsk->SetMinimum(0);
      }

      if(debug) printf("set overlay\n");
      int ratio = 0;

      if(debug) printf("set hsk\n");
      hsk->SetMaximum(1.35*histmax);

      hsk->Draw("hist");
      hsk->GetXaxis()->SetTitle(unit[i] == "" ? titleX[i].Data() : (titleX[i] + " [" + unit[i] + "]").Data());
      hsk->GetXaxis()->SetLabelColor(kWhite);
      char str[30];
      sprintf(str,"Events / %4.2f %s",binwidth(i)*rebin[i], unit[i].Data());
      hsk->GetYaxis()->SetTitle(str);
      hsk->GetYaxis()->SetTitleOffset(1.6);
      hsk->GetYaxis()->SetLabelSize(hsk->GetYaxis()->GetLabelSize()*0.7);
      hsk->GetXaxis()->SetLabelSize(hsk->GetXaxis()->GetLabelSize()*0.7);
      hsk->GetYaxis()->SetTitleSize(hsk->GetYaxis()->GetTitleSize()*0.7);
      if(debug) printf("set blinding\n");

      for(Int_t j=1; j<nbin[i]+1; j++) {
        hmcR.SetBinContent(j,1);
        hmcR.SetBinError(j,hmc.GetBinContent(j)>0 ? hmc.GetBinError(j)/hmc.GetBinContent(j) : 0);
        if(dataref) hdataR.SetBinContent(j, hmc.GetBinContent(j)>0 ? datahist->GetBinContent(j)/hmc.GetBinContent(j) : 1);
        if(dataref) hdataR.SetBinError(j, ( datahist->GetBinContent(j)>0 && hmc.GetBinContent(j)>0 )? datahist->GetBinError(j)/hmc.GetBinContent(j) : 0);
      }

      if(debug) printf("setting hmcR\n");
      hmc.SetFillColor(1);
      hmc.SetLineColor(0);
      hmc.SetMarkerStyle(1);
      hmc.SetMarkerSize(0);
      hmc.SetMarkerColor(1);
      hmc.SetFillStyle(3004);
      hmc.Draw("E2 same");
      lg1->Draw("same");

      if(debug) printf("atlas label\n");
      ATLASLabel(0.15,0.900,workflow.Data(),kBlack,lumi.Data(), analysis.Data(), region.Data());

//===============================lower pad===============================
      padlow->SetFillStyle(4000);
      padlow->SetGrid(1,1);
      padlow->SetTopMargin(0.03);
      padlow->SetBottomMargin(0.35);
      padlow->SetRightMargin(0.08);
      padlow->SetLeftMargin(0.12);
      padlow->cd();

      if(debug) printf("plot data ratio\n");
      if(dataref) {
        hdataR.SetMarkerStyle(20);
        hdataR.SetMarkerSize(0.8);
        hdataR.SetMaximum(1.5);
        hdataR.SetMinimum(0.5);
        hdataR.GetYaxis()->SetNdivisions(504,false);
        hdataR.GetYaxis()->SetTitle("Data/Bkg");
        hdataR.GetYaxis()->SetTitleOffset(hdataR.GetYaxis()->GetTitleOffset()*1.08);
        hdataR.GetYaxis()->CenterTitle();
        hdataR.GetXaxis()->SetTitle(unit[i] == "" ? titleX[i].Data() : (titleX[i] + " [" + unit[i] + "]").Data());
        hdataR.GetXaxis()->SetTitleSize(hdataR.GetXaxis()->GetTitleSize()*0.7);
        hdataR.GetYaxis()->SetTitleSize(hdataR.GetYaxis()->GetTitleSize()*0.7);
      }
      hmcR.SetFillColor(1);
      hmcR.SetLineColor(0);
      hmcR.SetMarkerStyle(1);
      hmcR.SetMarkerSize(0);
      hmcR.SetMarkerColor(1);
      hmcR.SetFillStyle(3004);
      if(debug) printf("plot data ratio\n");
      hdataR.Draw("E same");
      hdataR.GetXaxis()->SetTitleOffset(3.4);
      hdataR.GetXaxis()->SetLabelSize(hdataR.GetXaxis()->GetLabelSize()*0.7); 
      hdataR.GetYaxis()->SetLabelSize(hdataR.GetYaxis()->GetLabelSize()*0.7); 
      if(debug) printf("plot mc ratio\n");
      hmcR.Draw("E2 same");
      TLine line;
      line.SetLineColor(2);
      line.DrawLine(hdataR.GetBinLowEdge(1), 1., hdataR.GetBinLowEdge(hdataR.GetNbinsX()+1), 1.);
      cv.cd();
      if(debug) printf("draw low pad\n");
      padlow->Draw();
      if(debug) printf("printing\n");

//===============================upper pad signal and blinded data===============================

      padhi->cd();
      if(!overlaysamples.size()) {
        if(dataref) datahist->Draw("E same");
        cv.SaveAs(outputdir + "/" + region + "/" + name[i] + ".pdf");
      }
      if(sensitivevariable == name[i]) printf("region %s, background yield: %4.2f\n", region.Data(), hmc.Integral());

      for(auto overlaysample: overlaysamples){
        
        TLegend *lgsig = (TLegend*) lg1->Clone();
        if(debug) { printf("overlay: %s\n", overlaysample.Data()); }
        if(grabhist(overlaysample,region,i)) histoverlay = (TH1D*)grabhist(overlaysample,region,i)->Clone();
        if(doROC && sensitivevariable == name[i]) ROC_sig = (TH1D*) histoverlay->Clone();
        if(!histoverlay) continue;
        if(rebin[i] != 1) histoverlay->Rebin(rebin[i]);
        histoverlay->SetLineStyle(9);
        histoverlay->SetLineWidth(3);
        histoverlay->SetLineColor(kRed);
        histoverlay->SetFillColor(0);
        histoverlay->SetMinimum(0);
        ratio = histmax/histoverlay->GetMaximum()/3;
        if(ratio>10) ratio -= ratio%10;
        if(ratio>100) ratio -= ratio%100;
        if(ratio>1000) ratio -= ratio%1000;
        lgsig->AddEntry(histoverlay,(histoverlay->GetTitle() + (ratio > 0? "#times" + to_string(ratio) : "")).c_str(),"LP");

        if(sensitivevariable == name[i]){
          double _significance = 0;
          for(Int_t j=1; j<nbin[i]+1; j++) {
            if(histoverlay->GetBinContent(j) && hmc.GetBinContent(j)) {
              _significance += pow(significance(hmc.GetBinContent(j), histoverlay->GetBinContent(j)),2);
            }
          }
          if(doROC && sensitivevariable == name[i]){
            double bkgintegral = ROC_bkg->Integral();
            double sigintegral = ROC_sig->Integral();
            double sigeff = 1;
            double bkgrej = 0;
            ROC->SetPoint(0,sigeff,bkgrej);
            for (int i = 1; i < ROC_sig->GetNbinsX()+1; ++i)
            {
              sigeff -= ROC_sig->GetBinContent(i)/sigintegral;
              bkgrej += ROC_bkg->GetBinContent(i)/bkgintegral;
              ROC->SetPoint(i,sigeff,bkgrej);
            }
            outputfile->cd();
            ROC->Write(overlaysample + "_ROC");
            ROC_sig->Write(overlaysample + "_ROC_sig");
            ROC_bkg->Write(overlaysample + "_ROC_bkg");
            deletepointer(ROC);
            deletepointer(ROC_sig);
            deletepointer(ROC_bkg);
          }
          printf("signal %s yield: %4.2f, significance: %4.2f\n",overlaysample.Data(), histoverlay->Integral(), sqrt(_significance));
        }
        if(dataref) datahistblinded = (TH1D*)datahist->Clone();

        if(blinding && dataref){
          for(Int_t j=1; j<nbin[i]+1; j++) {
            if(histoverlay->GetBinContent(j)/sqrt(datahistblinded->GetBinContent(j)) > blinding) {
              datahistblinded->SetBinContent(j,0);
              datahistblinded->SetBinError(j,0);
            }
          }
          if(sensitivevariable == name[i]){
            for(int j = nbin[i]*3/4/rebin[i] ; j <= nbin[i] ; j++){
              datahistblinded->SetBinContent(j,0);
              datahistblinded->SetBinError(j,0);
            }
          }
        }
        if(ratio > 0) histoverlay->Scale(ratio);
        if(overlaysample != "") histoverlay->Draw("hist same");
        if(dataref) datahistblinded->Draw("E same");
        lgsig->SetBorderSize(0);
        lgsig->Draw();
        padhi->Update();
        cv.SaveAs(outputdir + "/" + region + "/" + name[i] + ".pdf");
        deletepointer(histoverlay);
        deletepointer(lgsig);
        if(dataref) deletepointer(datahistblinded);
      }
      deletepointer(hsk);
      deletepointer(lg1);
      deletepointer(padlow );
      deletepointer(padhi  );
      deletepointer(datahist);
      for(auto &iter : buffer) deletepointer(iter);
      if(debug) printf("end region %s\n",region.Data());
      cv.SaveAs(outputdir + "/" + region + "/" + name[i] + ".pdf]");
      cv.Clear();
    }
    if(debug) printf("end loop region\n");
  }
}
