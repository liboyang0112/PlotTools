#include "histSaver.h"
#include "fcnc_include.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "AtlasStyle.h"
#include "AtlasLabels.h"
#include "HISTFITTER.h"
#include "LatexChart.h"

using namespace std;
histSaver::histSaver(TString _outputfilename) {
  outputfilename = _outputfilename;
  trexdir = "trexinputs";
  inputfilename = "hists";
  nominalfilename = "";
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
  debug = 1;
  sensitivevariable = "";
}

histSaver::~histSaver() {
  if(debug) printf("histSaver::~histSaver()\n");
  for(auto& samp : plot_lib){
    for(auto &reg: samp.second) {
      for(auto &variation: reg.second) {
        for (int i = 0; i < v.size(); ++i){
          TH1D *target = variation.second[i];
            deletepointer(target);
          if(debug) cout<<"\rdone deleting histogram"<<std::endl<<std::flush;
        }
      }
    }
  }
  if(debug) std::cout<<"plot_lib destructed"<<std::endl;
  deletepointer(inputfile);
  if(debug) std::cout<<"inputfile destructed"<<std::endl;
  for(auto &file : outputfile)
    deletepointer(file.second);
  outputfile.clear();
  printf("histSaver::~histSaver() destructed\n");
}

void histSaver::printyield(TString region){
  printf("Print Yeild: %s\n", region.Data());
  double er;
  for(auto iter: plot_lib){
    TH1D* target = grabhist_int(iter.first,region,0,0);
    if(target){
      printf("%s: %4.3f \\pm %4.3f\n", iter.first.Data(), target->IntegralAndError(1,target->GetNbinsX(), er), er);
    }else{
      printf("Warning: histogram not found: %s, %s, %s\n", iter.first.Data(), region.Data(), v[0]->name.Data());
    }
  }
}

observable histSaver::calculateYield(TString region, string formula, TString variation){
  observable yield;
  istringstream iss(formula);
  vector<string> tokens{istream_iterator<string>{iss},
    istream_iterator<string>{}};
  if(tokens.size()%2) printf("Error: Wrong formula format: %s\nShould be like: 1 data -1 real -1 zll ...", formula.c_str());
  vector<TH1D*> newvec;
  observable scaleto(0,0);

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
    TString sample=tokens[icompon+1].c_str();
    if(grabhist(sample,region, sample == "data" ? "NOMINAL" : variation,0)){
      TH1D *target=grabhist_int(sample,region,0,0);
      double err;
      observable thisyield(target->IntegralAndError(1,target->GetNbinsX(),err),err);
      yield+=thisyield*numb;
    }
  }
  printf("histSaver::calculateYield() : Calcuated yield: %f+/-%f in %s\n", yield.nominal,yield.error,region.Data());
  return yield;
}

int histSaver::findvar(TString varname){
  for (int i = 0; i < v.size(); ++i)
  {
    if(v.at(i)->name == varname) return i;
  }
  printf("varname not found: %s\n", varname.Data());
  exit(0);
}

TH1D* histSaver::grabhist_int(TString sample, TString region, int ivar, bool vital){
  return grabhist(sample, region, "NOMINAL", ivar, vital);
}

TH1D* histSaver::grabhist(TString sample, TString region, TString variation, TString varname, bool vital){
  int ivar = -1;
  for (int i = 0; i < v.size(); ++i)
  {
    if(varname == v.at(i)->name){
      ivar = i;
      break;
    }
  }
  return grabhist(sample, region, variation, ivar, vital);
}

TH1D* histSaver::grabhist(TString sample, TString region, TString variation, int ivar, bool vital){
  if(sample == "data") variation = "NOMINAL";
  auto samp = plot_lib.find(sample);
  if(samp == plot_lib.end()) {
    if(debug) printf("histSaver:grabhist  Warning: sample %s not found\n", sample.Data());
    if(debug) show();
    if(vital) exit(0);
    return 0;
  }
  auto reg = samp->second.find(region);
  if(reg == samp->second.end()){
    if(debug) printf("histSaver:grabhist  Warning: region %s for sample %s not found\n", region.Data(), sample.Data());
    if(debug) show();
    if(vital) exit(0);
    return 0;
  }
  auto vari = reg->second.find(variation);
  if(vari == reg->second.end()){
    if(debug) printf("histSaver:grabhist  Warning: region %s for sample %s with variation %s not found\n", region.Data(), sample.Data(), variation.Data());
    if(debug) show();
    if(vital) exit(0);
    return 0;
  }
  return vari->second.at(ivar);
}

TH1D* histSaver::grabhist(TString sample, TString region, TString varname, bool vital){
  int ivar = -1;
  for (int i = 0; i < v.size(); ++i)
  {
    if(varname == v.at(i)->name){
      ivar = i;
      break;
    }
  }
  return grabhist(sample, region, ivar, vital);
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

Float_t histSaver::getVal(Int_t i) {
  Float_t tmp = -999999;
  if(address1[i])      { if(debug) printf("fill address1\n"); tmp = *address1[i]*v.at(i)->scale; }
  else if(address3[i]) { if(debug) printf("fill address3\n"); tmp = *address3[i]*v.at(i)->scale; }
  else if(address2[i]) { if(debug) printf("fill address2\n"); tmp = *address2[i]; }
  else {
    printf("error: fill variable failed. no type available for var %s\n",v.at(i)->name.Data());
    exit(0);
  }
  if(debug == 1) printf("fill value: %4.2f\n", tmp);
  if (!v.at(i)->xbins){
    if(tmp >= v.at(i)->xhigh) tmp = v.at(i)->xhigh*0.999999;
    if(tmp < v.at(i)->xlow) tmp = v.at(i)->xlow;
  }else{
    double xhi = v.at(i)->xbins->at(v.at(i)->xbins->size()-1);
    double xlow = v.at(i)->xbins->at(0);
    if(tmp >= xhi) tmp = xhi*0.999999;
    if(tmp < xlow) tmp = xlow;
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
  if(address1.size() == v.size()) {
    for (int i = 0; i < v.size(); ++i)
    {
      if(address2[i]) printf("histSaver::show()\t%s = %d\n", v.at(i)->name.Data(), *address2[i]);
      else if(address1[i]) printf("histSaver::show()\t%s = %4.2f\n", v.at(i)->name.Data(), *address1[i]*v.at(i)->scale);
      else if(address3[i]) printf("histSaver::show()\t%s = %4.2f\n", v.at(i)->name.Data(), *address3[i]*v.at(i)->scale);
    }
  }
}

float histSaver::binwidth(int i){
  return (v.at(i)->xhigh-v.at(i)->xlow)/v.at(i)->nbins;
}


void histSaver::merge_regions(vector<TString> inputregions, TString outputregion){
  if(debug) printf("histSaver::merge_regions\t");
  bool exist = 0;
  vector<TString> existregions;
  for(auto& iter:plot_lib ){
    if(debug) printf("=====================start merging sample %s=====================\n", iter.first.Data());
    existregions.clear();
    for(auto region:inputregions){
      if(iter.second.find(region) != iter.second.end()){
        existregions.push_back(region);
      }else if(debug) printf("histSaver::merge_regions\t input region: %s not found for sample %s\n",region.Data(), iter.first.Data());
    }
    if(existregions.size() == 0) continue;
    bool outputexist = 0;
    auto outiter = iter.second.find(outputregion);
    if(outiter != iter.second.end()){
      if(debug) printf("histSaver::merge_regions\t outputregion %s exist, overwrite it\n",outputregion.Data());
      exist = 1;
      for(auto &variation: outiter->second){
        for(int i = 0; i < v.size(); ++i) deletepointer(variation.second[i]);
        variation.second.clear();
      }
    }
    for(auto &variation : iter.second[existregions[0]]){
      for (int i = 0; i < v.size(); ++i)
      {
        auto &tmpiter = iter.second[outputregion][variation.first];
        tmpiter.push_back(0);
        for(auto region:existregions){
          TH1D* addtarget = grabhist(iter.first,region,variation.first,i);
          if(addtarget){
            if(tmpiter[i] == 0) tmpiter[i] = (TH1D*)addtarget->Clone(iter.first + "_" + variation.first+"_"+outputregion+"_"+v.at(i)->name + "_buffer");
            else tmpiter[i]->Add(addtarget);
          }
        }
      }
    }
  }
  if(!exist) regions.push_back(outputregion);
}
void histSaver::merge_regions(TString inputregion1, TString inputregion2, TString outputregion){
  if(debug) printf("histSaver::merge_regions\t %s and %s into %s\n",inputregion1.Data(),inputregion2.Data(),outputregion.Data());
  bool exist = 0;
  bool input1exist = 1;
  bool input2exist = 1;

  for(auto& iter:plot_lib ){
    if(debug) printf("=====================start merging sample %s=====================\n", iter.first.Data());
    input1exist = 1;
    input2exist = 1;
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
      for(auto &variation: iter.second[outputregion]){
        for(int i = 0; i < v.size(); ++i) deletepointer(variation.second[i]);
        variation.second.clear();
      }
    }
    for(auto &variation : iter.second[inputregion1])
    for (int i = 0; i < v.size(); ++i)
    {
      TH1D* addtarget1 = grabhist(iter.first,inputregion1,variation.first,i);
      TH1D* addtarget2 = grabhist(iter.first,inputregion2,variation.first,i);
      auto &tmpiter = iter.second[outputregion][variation.first];
      if(input1exist == 1 && addtarget1) tmpiter.push_back((TH1D*)addtarget1->Clone(iter.first + "_" + variation.first+"_"+outputregion+"_"+v.at(i)->name + "_buffer"));
      else if(addtarget2) tmpiter.push_back((TH1D*)addtarget2->Clone(iter.first + "_" + variation.first+"_"+outputregion+"_"+v.at(i)->name + +"_buffer"));
      else tmpiter.push_back(0);
      if(input1exist == 1 && input2exist == 1 && addtarget1 && addtarget2) {
        tmpiter[i]->Add(addtarget2);
        if(debug)
          printf("add %s to %s as %s\n", iter.second[inputregion2][variation.first][i]->GetName(),iter.second[inputregion1][variation.first][i]->GetName(),tmpiter[i]->GetName());
      }
    }
  }
  if(!exist) regions.push_back(outputregion);
}

void histSaver::add_sample(TString samplename, TString sampletitle, enum EColor color){
  samples.emplace_back(samplename,sampletitle,color,1);
  plot_lib[samplename];
}
void histSaver::init_hist(map<TString,map<TString,map<TString,vector<TH1D*>>>>::iterator sample_lib, TString region, TString variation){
  createdNP = variation;
  if(outputfile.find(variation) == outputfile.end()) outputfile[variation] = new TFile(outputfilename + "_" + variation + ".root","update");
  else outputfile[variation]->cd();
  auto sampleiter = find_if(samples.begin(),samples.end(),[sample_lib](fcncSample const& tmp){return sample_lib->first == tmp.name;});
  if(debug) {
    printf("histSaver::init_hist() : add new sample: %s\n", sample_lib->first);
    printf("histSaver::init_hist() : variation: %s\n",variation.Data());
  }
  if(debug) printf("histSaver::init_hist(): Region: %s\n", region.Data());
  pair<TString,map<TString,vector<TH1D*>>> regionpair;
  regionpair.first = region;
  auto regioniter = sample_lib->second.insert(regionpair).first;
  regionpair.first = region;
  for (int i = 0; i < v.size(); ++i){
    if(debug) printf("histSaver::init_hist(): var %s\n",v.at(i)->name.Data());
    TString histname = sample_lib->first + "_" + variation  + "_" +  region + "_" + v.at(i)->name + "_buffer";
    TH1D *created = new TH1D(histname,sampleiter->title,v.at(i)->nbins,v.at(i)->xlow,v.at(i)->xhigh);
    created->SetDirectory(0);
    regioniter->second[variation].push_back(created);
    if (sample_lib->first != "data")
    {
      created->Sumw2();
      created->SetFillColor(sampleiter->color);
      created->SetLineWidth(1);
      created->SetLineColor(kBlack);
      created->SetMarkerSize(0);
    }
  }
  
  if(debug == 1) printf("plot_lib[%s][%s][%s]\n", sample_lib->first.Data(), region.Data(), variation.Data());
  
  if (sample_lib->first == "data") dataref = 1;

  if(debug) printf("finished initializing %s\n", sample_lib->first.Data() );
}

vector<observable> histSaver::scale_to_data(TString scaleregion, string formula, TString scaleVariable, vector<double> slices, TString variation){
  int nslice = slices.size();
  int ivar = 0;
  for (; ivar < v.size(); ++ivar)
  {
    if(v[ivar]->name == scaleVariable) break;
  }

  if(!nslice) {
    nslice = v[ivar]->nbins+1;
    slices.push_back(v[ivar]->xlow);
    for (int i = 0; i < v[ivar]->nbins; ++i)
    {
      slices.push_back(slices[i]+binwidth(ivar));
    }
  }

  if(outputfile.find(variation) == outputfile.end()) outputfile[variation] = new TFile(outputfilename + "_" + variation + ".root", "update");
  else outputfile[variation]->cd();
  vector<TString> tokens = split(formula.c_str()," ");
  if(tokens.size()%2) printf("Error: Wrong formula format: %s\nShould be like: 1 real 1 zll ...", formula.c_str());
  vector<observable> scalefrom;
  vector<observable> scaleto;
  for(int i = 0 ; i < nslice ; i++)
    scalefrom.push_back(observable(0,0));
  scaleto = scalefrom;
  for(auto &sample: plot_lib){
    TH1D *target = grabhist(sample.first,scaleregion,variation,scaleVariable);
    int islice = 0;
    if(target){
      auto iter = find(tokens.begin(),tokens.end(), sample.first.Data());
      if(iter != tokens.end())
      {
        double numb = 0;
        try{
          numb = stof(string((*(iter-1)).Data()));
        } 
        catch(const std::invalid_argument& e){
          printf("Error: Wrong formula format: %s\nShould be like: 1 real 1 zll ...", formula.c_str());
          exit(1);
        }
        
        if(target->GetBinLowEdge(0) > slices[0]) {
          printf("WARNING: slice 1 (%4.2f, %4.2f) is lower than the low edge of the histogram %4.2f, please check variable %s\n", slices[0], slices[1], target->GetBinLowEdge(0), scaleVariable.Data());
        }
        for (int i = 1; i <= v[ivar]->nbins; ++i)
        {
          if(target->GetBinLowEdge(i) < slices[0]) continue;
          if(islice == nslice-1) break;
          scalefrom[islice] += observable(target->GetBinContent(i),target->GetBinError(i))*numb;
          if(target->GetBinLowEdge(i) >= slices[islice+1]) islice+=1;
        }
      }else{
        for (int i = 1; i <= v[ivar]->nbins; ++i)
        {
          if(target->GetBinLowEdge(i) < slices[0]) continue;
          if(islice == nslice-1) break;
          if(sample.first == "data") 
            scaleto[islice] += observable(target->GetBinContent(i),target->GetBinError(i));
          else scaleto[islice] -= observable(target->GetBinContent(i),target->GetBinError(i));
          if(target->GetBinLowEdge(i) >= slices[islice+1]) islice+=1;
        }
      }
    }
  }
  vector<observable> scalefactor;
  for (int i = 0; i < nslice; ++i)
  {
    scalefactor.push_back(scaleto[i]/scalefrom[i]);
  }
  printf("region %s, scale variable %s in %d slices:\n", scaleregion.Data(), scaleVariable.Data(), nslice);
  for (int i = 0; i < nslice-1; ++i)
    printf("(%4.2f, %4.2f): %4.2f +/- %4.2f to %4.2f +/- %4.2f, ratio: %4.2f +/- %4.2f\n",slices[i], slices[i+1],scalefrom[i].nominal,scalefrom[i].error,scaleto[i].nominal,scaleto[i].error,scalefactor[i].nominal,scalefactor[i].error);

  return scalefactor;
}

void histSaver::scale_sample(TString scaleregion, string formula, TString scaleVariable, vector<observable> scalefactor, vector<double> slices, TString variation){
  int nslice = slices.size();
  int ivar = 0;
  for (; ivar < v.size(); ++ivar)
  {
    if(v[ivar]->name == scaleVariable) break;
  }
  if(!nslice) {
    nslice = v[ivar]->nbins+1;
    slices.push_back(v[ivar]->xlow);
    for (int i = 0; i < v[ivar]->nbins; ++i)
    {
      slices.push_back(slices[i]+binwidth(ivar));
    }
  }
  vector<TString> tokens = split(formula.c_str()," ");
  if(tokens.size()%2) printf("Error: Wrong formula format: %s\nShould be like: 1 real 1 zll ...", formula.c_str());
  for (int i = 0; i < tokens.size(); ++i){
    if(!(i%2)) continue;
      int islice = 0;
      TH1D *target = grabhist(tokens[i],scaleregion,variation,scaleVariable);
      if(!target) {
        printf("histSaver::scale_sample : WARNING: hist not found grabhist(%s,%s,%s,%s)\n", tokens[i].Data(),scaleregion.Data(),variation.Data(),scaleVariable.Data());
        continue;
      }
      for (int i = 1; i <= v[ivar]->nbins; ++i)
      {
        if(target->GetBinLowEdge(i) < slices[0]) continue;
        if(islice == nslice-1) break;
        double scaletmp = scalefactor[islice].nominal;
        if(scaletmp<0) scaletmp = 0;
        target->SetBinContent(i,target->GetBinContent(i)*scaletmp);
        if(target->GetBinLowEdge(i) >= slices[islice+1]) islice+=1;
      }
  }
}

map<TString,vector<observable>>* histSaver::fit_scale_factor(vector<TString> *fit_regions, TString *variable, vector<TString> *scalesamples, const vector<double> *slices, TString *variation, vector<TString> *postfit_regions){
  auto *_scalesamples = new map<TString,map<TString,vector<TString>>>();
  auto *_postfit_regions = new map<TString,map<TString,vector<TString>>>();
  for(auto sample: *scalesamples){
    (*_scalesamples)[sample];
    (*_postfit_regions)[sample][sample] = *postfit_regions;
  }
  auto ret = fit_scale_factor(fit_regions, variable, _scalesamples, slices, variation, _postfit_regions);
  delete _scalesamples;
  delete _postfit_regions;
  return ret;
}

map<TString,vector<observable>>* histSaver::fit_scale_factor(vector<TString> *fit_regions, TString *variable, map<TString,map<TString,vector<TString>>> *scalesamples, const vector<double> *slices, TString *_variation, map<TString,map<TString,vector<TString>>> *postfit_regions){
  if(!postfit_regions) postfit_regions = scalesamples;
  auto *scalefactors = new map<TString,vector<observable>>();
  TString variation = _variation? *_variation:"NOMINAL";
  vector<observable> iter;
  int nbins = v[findvar(*variable)]->nbins;
  int ihists = 0;
  vector<int> binslices;
  HISTFITTER* fitter = new HISTFITTER();
  vector<TString> params;
  for(auto sample : *scalesamples) {
    if(sample.second.size()){
      for(auto sfForReg : sample.second){
        fitter->setparam("sf_" + sample.first + "_" + sfForReg.first, 1, 0.1, 0.,2.);
        printf("fit region: ");
        for(auto regions : sfForReg.second) printf(" %s ", regions.Data());
        printf("\n");
        params.push_back("sf_" + sample.first + "_" + sfForReg.first);
      }
    }else{
      fitter->setparam("sf_" + sample.first, 1, 0.1, 0.,2.);
      printf("fit region: All\n");
      params.push_back("sf_" + sample.first);
    }
  }
  for (int i = 0; i < slices->size()-1; ++i)
  {
    auto fitsamples = stackorder;
    fitsamples.push_back("data");
    for(auto sample : fitsamples){
      for(auto reg : *fit_regions){
        TH1D *target = grabhist(sample,reg,variation,*variable);
        if(!target) continue;
        if(ihists == 0) {
          binslices = resolveslices(target,slices);
        }
        bool scale = 0;
        TString SFname = "";
        TString addsample = sample;
        for(auto ssample : *scalesamples) {
          if(ssample.first == sample) {
            if(ssample.second.size()){
              for(auto sfForReg: ssample.second){
                for(auto sfreg: sfForReg.second){
                  if (sfreg == reg)
                  {
                    addsample = sample + "_" + sfForReg.first;
                  }
                }
              }
            }
            SFname = "sf_" + addsample;
          }
        }
        fitter->addfithist(sample,target,binslices[i],binslices[i+1]-1,SFname);
        ihists++;
      }
    }
    Double_t val[10],err[10];
    //fitter->debug();
//============================ do fit here============================
    //fitter->asimovfit(100,nprong[iprong]+"ptbin"+char(ptbin+'0')+".root");
    double chi2 = fitter->fit(val,err,0);
    int ipar = 0;
    for (auto par: params)
    {
      (*scalefactors)[par].push_back(observable(val[ipar],err[ipar]));
      ipar++;
    }
    fitter->clear();
  }
  for(auto samp : *postfit_regions){
    TH1D *target;
    map<TString,vector<TString>> plotregions;
    if(samp.second.size() == 0){
      plotregions["sf_" + samp.first] = *fit_regions;
    }else if(samp.second.size() == 1){
      for(auto sfForReg : samp.second)
        plotregions["sf_" + samp.first] = sfForReg.second;
    }else{
      for(auto sfForReg : samp.second)
        plotregions["sf_" + samp.first + "_" + sfForReg.first] = sfForReg.second;
    }

    for(auto sf : plotregions){
      for(auto reg : sf.second){
        target = grabhist(samp.first,reg,variation,*variable);
        if(!target) continue;
        for (int islice = 0; islice < slices->size()-1; ++islice)
        {
          for (int i = binslices[islice]; i < binslices[islice+1]; ++i)
          {
            target->SetBinContent(i,target->GetBinContent(i) * (*scalefactors)[sf.first][islice].nominal);
            target->SetBinError(i,target->GetBinError(i) * (*scalefactors)[sf.first][islice].nominal);
          }
        }
      }
    }
  }
  printf("fit regions:");
  for(auto reg: *fit_regions){
    printf(" %s ", reg.Data());
  }
  printf("\n");
  printf("fit samples:");
  for(auto param: params){
    printf(" %s ", param.Data());
  }
  printf("\n");
  for (int i = 0; i < slices->size()-1; ++i){
    printf("(%4.2f, %4.2f): ",(*slices)[i], (*slices)[i+1]);
    for(auto par: params){
      printf("%4.2f +/- %4.2f, ", (*scalefactors)[par][i].nominal, (*scalefactors)[par][i].error);
    }
    printf("\n");
  }

  return scalefactors;
}

vector<int> histSaver::resolveslices(TH1D* target, const vector<double> *slices){
  printf("histSaver::resolveslices(): total %d bins, from %f, to %f.\n",target->GetNbinsX(),target->GetBinLowEdge(1), target->GetBinLowEdge(target->GetNbinsX())+target->GetBinWidth(target->GetNbinsX()));
  if(debug){
    printf("slices ( ");
    for (int i = 0; i < slices->size(); ++i)
    {
      printf("%f ", slices->at(i));
    }
    printf(")\n");
  }
  vector<int> ret;
  if(target->GetBinLowEdge(1) > slices->at(0) ){
    printf("WARNING: slice 1 (%4.2f, %4.2f) is lower than the low edge of the histogram %4.2f, please check histogram %s\n", slices->at(0), slices->at(1), target->GetBinLowEdge(0), target->GetName());
  }
  if(target->GetXaxis()->GetXmax() < slices->at(slices->size()-1)) {
          printf("WARNING: last slice (%4.2f, %4.2f) is lower than the low edge of the histogram %4.2f, please check histogram %s\n", slices->at(slices->size()-2), slices->at(slices->size()-1), target->GetXaxis()->GetXmax(), target->GetName());
  }
  int islice = 0;
  for (int i = 1; i <= target->GetNbinsX(); ++i)
  {
    if(target->GetBinLowEdge(i) >= slices->at(islice)) {
      ret.push_back(i);
      islice+=1;
    }
    if(islice == slices->size()) break;
  }
  if(target->GetXaxis()->GetXmax() <= slices->at(slices->size()-1)) ret.push_back(target->GetNbinsX()+1);
  if(debug){
    printf("resolved slices ( ");
    for (int i = 0; i < ret.size(); ++i)
    {
      printf("%d ", ret.at(i));
    }
    printf(")\n");
  }
  return ret;
}

void histSaver::read_sample(TString samplename, TString savehistname, TString variation, TString sampleTitle, enum EColor color, double norm, TFile *_inputfile){

  TFile *readfromfile;

  if(_inputfile) readfromfile = _inputfile;
  else{
    if(!inputfile) inputfile = new TFile(inputfilename + ".root", "read");
    if(!inputfile) inputfile = new TFile(nominalfilename + ".root", "read");
    readfromfile = inputfile;
  }
  TString filename(readfromfile->GetName());
  if (debug == 1) printf("read from file: %s\n", filename.Data());
  if (samplename == "data") dataref = 1;
  TString histnameorig(savehistname + "_");
  for(auto const& region: regions) {
    TString histname;
    if(filename.Contains("NOMINAL") && variation.Contains("Xsec")){
      histname = histnameorig + "NOMINAL_" + region + "_";
    }else{
      histname = histnameorig + variation + "_" + region + "_";
    }
    if (debug == 1)
    {
      printf("read sample %s from %s region\n", samplename.Data(), region.Data());
    }
    auto &samplib = plot_lib[samplename];
    bool newRegion = samplib.find(region) == samplib.end();
    for (int i = 0; i < v.size(); ++i)
    {
      if(debug) printf("histSaver::read_sample() : Read file %s to get %s\n",readfromfile->GetName(), (histname + v.at(i)->name).Data());
      TH1D *readhist = (TH1D*)readfromfile->Get(histname + v.at(i)->name);
      if(!readhist) {
        if(debug) printf("histogram name not found: %s\n", (histname + v.at(i)->name).Data());
        if(i == samplib[region][variation].size()) samplib[region][variation].push_back(0);
        continue;
      }
      double tmp = readhist->Integral();
      if(tmp!=tmp){
        printf("Warning: %s->Integral() is nan, skip\n", (histname + v.at(i)->name).Data());
        if(i == samplib[region][variation].size()) samplib[region][variation].push_back(0);
        continue;
      }
      if(tmp==0){
        printf("Warning: %s->Integral() is 0, skip\n", (histname + v.at(i)->name).Data());
        if(i == samplib[region][variation].size()) samplib[region][variation].push_back(0);
        continue;
      }
      if(checkread){
        if(samplename == checkread_sample && region == checkread_region && variation == checkread_variation && i == checkread_variable){
          printf("read histogram %s, + %f\n", (histname + v.at(i)->name).Data(), readhist->GetBinContent(checkread_ibin)*norm);
        }
      }
      bool newhist = 0;
      TH1D* target;
      if(newRegion) {
        target = (TH1D*)(readfromfile->Get(histname + v.at(i)->name)->Clone());
        samplib[region][variation].push_back(target);
      }else {
        auto &tmp = samplib[region][variation];
        target = tmp[i];
        if(!target) {
          target = (TH1D*)(readfromfile->Get(histname + v.at(i)->name)->Clone());
          tmp[i] = target;
        }else{
          target->Add(readhist,norm);
          continue;
        }
      }
      target->SetName(samplename + "_" + variation + "_" + region + "_" + v.at(i)->name + "_buffer");
      target->Scale(norm);
      target->SetTitle(sampleTitle);
      target->SetFillColorAlpha(color,1);
      target->SetLineWidth(1);
      target->SetLineColor(kBlack);
      target->SetMarkerSize(0);
      target->SetDirectory(0);
    }
    if(debug) printf("histSaver::read_sample : finish read plot_lib[%s][%s][%s][%d]", samplename.Data(),region.Data(),variation.Data(),samplib[region][variation].size());
  }
}

void histSaver::add_region(TString region){
  regions.push_back(region);
  nregion += 1;
}

void histSaver::fill_hist(TString sample, TString region, TString variation){
  auto sampleiter = plot_lib.find(sample);
  if(sampleiter == plot_lib.end()) {
    printf("histSaver::fill_hist() ERROR: sample %s not found\n", sample.Data());
    show();
    exit(0);
  }
  if (weight_type == 0)
  {
    printf("ERROR: weight not set\n");
  }
  if(sampleiter->second.find(region) == sampleiter->second.end()){
    init_hist(sampleiter, region, variation);
    if(find(regions.begin(), regions.end(), region) == regions.end()) {
      regions.push_back(region);
      nregion += 1;
    }
  }


  for (int i = 0; i < v.size(); ++i){
    double fillval = getVal(i);
    if(fillval!=fillval) {
      printf("Warning: fill val is nan: \n");
      printf("plot_lib[%s][%s][%d]->Fill(%4.2f,%4.2f)\n", sample.Data(), region.Data(), i, fillval, weight_type == 1? *fweight : *dweight);
    }
    if(debug == 1) printf("plot_lib[%s][%s][%s][%d]->Fill(%4.2f,%4.2f)\n", sample.Data(), region.Data(), variation.Data(), i, fillval, weight_type == 1? *fweight : *dweight);
    TH1D *target = grabhist(sample,region,variation,i);
    if(target) target->Fill(fillval,weight_type == 1? *fweight : *dweight);
    else {
      if(!add_variation(sample,region,variation)) printf("add variation %s failed, sample %s doesnt exist\n", variation.Data(), sample.Data());
      TH1D *target = grabhist(sample,region,variation,i);
      if(target) target->Fill(fillval,weight_type == 1? *fweight : *dweight);
      else {
        printf("add_variation didnt work in filling sample %s, region %s, variation %s, variable %s\n",sample.Data(),region.Data(),variation.Data(),v[i]->name.Data());
        show();
        exit(0);
      }
    }
  }
}

void histSaver::fill_hist(TString sample, TString region){
  fill_hist(sample, region, "NOMINAL");
}


bool histSaver::find_sample(TString sample){
  if(plot_lib.find(sample) == plot_lib.end()) return 0;
  return 1;
}

bool histSaver::add_variation(TString sample,TString reg,TString variation){
  if(!find_sample(sample)) return 0;
  if(outputfile.find(variation) == outputfile.end()) outputfile[variation] = new TFile(outputfilename + "_" + variation + ".root", "update");
  else outputfile[variation]->cd();
  for (int i = 0; i < v.size(); ++i){
    if(plot_lib[sample][reg].begin() == plot_lib[sample][reg].end()) {
      printf("histSaver::add_variation() ERROR: No variation defined yet, cant add new variation\n");
      exit(0);
    }
    TH1D *created = (TH1D*) plot_lib[sample][reg][createdNP].at(i);
    if(!created){
      printf("histSaver::add_variation() ERROR: hist doesn't exist: plot_lib[%s][%s][%s][%d]\n",sample.Data(), reg.Data(), plot_lib[sample][reg].begin()->first.Data(),i);
      exit(0);
    }
    created = (TH1D*) created->Clone(sample + "_" + variation + "_" + reg + "_" + v.at(i)->name + "_buffer");
    created->Reset();
    created->SetDirectory(0);
    plot_lib[sample][reg][variation].push_back(created);
  }
  return 1;
}

void histSaver::write(){
  for(auto& iter: outputfile){
    for(auto& sample : plot_lib){
      for(auto& region: sample.second) {
        for(auto& variation : region.second){
          if(variation.first != iter.first) {
            continue;
          }
          double sum = variation.second[0]->Integral();
          if(sum == 0) continue;
          if(sum != sum) {
            printf("Warning: hist integral is nan, skip writing for %s\n", variation.second[0]->GetName());
            continue;
          }
          outputfile[variation.first]->cd();
          for (int i = 0; i < v.size(); ++i){
            if(variation.second[i]->GetMaximum() == sum) {
              continue;
            }
            //if(grabhist(iter.first,region,i)->Integral() == 0) {
            //  printf("Warning: histogram is empty: %s, %s, %d\n", iter.first.Data(),region.Data(),i);
            //}
            TString writename = variation.second[i]->GetName();
            writename.Remove(writename.Sizeof()-8,7); //remove "_buffer"
            printf("write histogram: %s\n", writename.Data());
            variation.second[i]->Write(writename,TObject::kWriteDelete);
          }
        }
      }
    }
    iter.second->Close();
    printf("histSaver::write() Written to file %s\n", iter.second->GetName());
  }
}

void histSaver::write_trexinput(TString NPname, TString writename, TString writeoption){
  if (writename == "")
  {
    writename = NPname;
  }
  gSystem->mkdir(trexdir);
  for (int i = 0; i < v.size(); ++i){
    gSystem->mkdir(trexdir + "/" + v.at(i)->name);
    for(auto const& region: regions) {
      bool muted = 0;
      for (auto const& mutedregion: mutedregions)
      {
        if(region.Contains(mutedregion))
          muted = 1;
      }
      if(muted) continue;

      gSystem->mkdir(trexdir + "/" + v.at(i)->name + "/" + region);
      for(auto& iter : plot_lib){
        if(NPname != "NOMINAL" && iter.first.Contains("data")) continue;
        TString filename = trexdir + "/" + v.at(i)->name + "/" + region + "/" + iter.first + ".root";
        TFile outputfile(filename, writeoption);
        if(debug) printf("Writing to file: %s, histoname: %s\n", filename.Data(), NPname.Data());
        TH1D *target = grabhist(iter.first,region,NPname,i);
        if(target) {
          target->Write(writename,TObject::kWriteDelete);
          if(!target->Integral()) printf("Warinig: plot_lib[%s][%s][%d] is empty\n", iter.first.Data(),region.Data(),i);
        }
        else if(debug) printf("Warning: histogram plot_lib[%s][%s][%d] not found\n", iter.first.Data(),region.Data(),i);
        outputfile.Close();
      }
    }
  }
}
void histSaver::clearhist(){
  if(debug) printf("histSaver::clearhist()\n");
  for(auto& sample : plot_lib){
    for(auto& region: sample.second) {
      for(auto& variation : region.second){
        for (int i = 0; i < v.size(); ++i){
          if(variation.second[i]){
            variation.second[i]->Reset();
          }else{
            printf("histSaver::Reset() Error: histogram not found: sample: %s, variable: %s, region: %s\n",sample.first.Data(), v.at(i)->name.Data(),region.first.Data());
          }
        }
      }
    }
  }
}

void histSaver::overlay(TString _overlaysample){
  overlaysamples.push_back(_overlaysample);
}

observable histSaver::templatesample(TString fromregion, TString variation,string formula,TString toregion,TString newsamplename,TString newsampletitle,enum EColor color, bool scaletogap, double SF){

  if(outputfile.find(variation) == outputfile.end()) outputfile[variation] = new TFile(outputfilename + "_" + variation + ".root", "update");
  else outputfile[variation]->cd();
  istringstream iss(formula);
  vector<string> tokens{istream_iterator<string>{iss},
    istream_iterator<string>{}};
  if(tokens.size()%2) printf("Error: Wrong formula format: %s\nShould be like: 1 data -1 real -1 zll ...", formula.c_str());
  vector<TH1D*> newvec;
  observable scaleto(0,0);
  for (int ivar = 0; ivar < v.size(); ++ivar)
  {
    TH1D *target = grabhist(tokens[1],fromregion, tokens[1] == "data" ? "NOMINAL" : variation,ivar);
    if(target){
      newvec.push_back((TH1D*)target->Clone(newsamplename+"_"+toregion+v[ivar]->name));
      newvec[ivar]->Reset();
      newvec[ivar]->SetNameTitle(newsamplename,newsampletitle);
      newvec[ivar]->SetFillColor(color);
    }else{
      newvec.push_back(0);
    }
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
    if(grabhist(tokens[icompon+1],toregion, tokens[icompon+1] == "data" ? "NOMINAL" : variation,0)){
      if(scaletogap) {
        double error = 0;
        observable tmp(grabhist(tokens[icompon+1],toregion, tokens[icompon+1] == "data" ? "NOMINAL" : variation,0)->Integral(),gethisterror(grabhist(tokens[icompon+1],toregion, tokens[icompon+1] == "data" ? "NOMINAL" : variation,0)));
        scaleto += tmp*numb;
      }
      for (int ivar = 0; ivar < v.size(); ++ivar)
      {
      	TH1D *target = grabhist(tokens[icompon+1],fromregion, tokens[icompon+1] == "data" ? "NOMINAL" : variation,ivar);
        if(target && newvec[ivar]) newvec[ivar]->Add(target,numb);
      }
    }
  }
  observable scalefactor;
  if(scaletogap) {
    observable scalefrom(newvec[0]->Integral(),gethisterror(newvec[0]));
    scalefactor = scaleto/scalefrom;
    printf("scale from %s: %4.2f +/- %4.2f\nto %s: %4.2f +/- %4.2f\nratio: %4.2f +/- %4.2f\n\n",
      fromregion.Data(), scalefrom.nominal, scalefrom.error,
      toregion.Data(), scaleto.nominal, scaleto.error,
      scalefactor.nominal, scalefactor.error);
    for(auto & hists : newvec){
      if(hists) hists->Scale(scalefactor.nominal);
    }
  }else{
    for(auto & hists : newvec){
      if(hists) hists->Scale(SF);
    }
  }
  for(int ivar = 0; ivar < v.size(); ivar++){
    plot_lib[newsamplename][toregion][variation] = newvec;
  }
  return scalefactor;
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
void histSaver::plot_stack(TString NPname, TString outdir, TString outputchartdir){
  SetAtlasStyle();
  TGaxis::SetMaxDigits(3);
  LatexChart* yield_chart = new LatexChart("yield");
  LatexChart* sgnf_chart = new LatexChart("significance");
  gSystem->mkdir(outdir);
  gSystem->mkdir(outputchartdir);
  TCanvas cv("cv","cv",600,600);
  TGraph* ROC;
  TH1D *ROC_sig = 0;
  TH1D *ROC_bkg = 0;
  TH1D* datahistorig = 0;
  vector<TH1D*> buffer;
  if(doROC){
    ROC = new TGraph();
    ROC -> SetName("ROC");
    ROC -> SetTitle("ROC");
  }
  TFile *outputrocfile = new TFile (outputfilename + "_roc.root", "recreate");
  for(auto const& region: regions) {
    bool muted = 0;
    for (auto const& mutedregion: mutedregions)
    {
      if(region.Contains(mutedregion))
        muted = 1;
    }
    if(muted) continue;
    gSystem->mkdir(outdir + "/" + region);
    for (int i = 0; i < v.size(); ++i){
      TPad *padlow = new TPad("lowpad","lowpad",0,0,1,0.3);
      TPad *padhi  = new TPad("hipad","hipad",0,0.3,1,1);
      TH1D hmc("hmc","hmc",v[i]->nbins/v[i]->rebin,v[i]->xlow,v[i]->xhigh);
      TH1D hmcR("hmcR","hmcR",v[i]->nbins/v[i]->rebin,v[i]->xlow,v[i]->xhigh);
      TH1D hdataR("hdataR","hdataR",v[i]->nbins/v[i]->rebin,v[i]->xlow,v[i]->xhigh);
//===============================upper pad bkg and unblinded data===============================
      hmc.Sumw2();
      THStack *hsk = new THStack(v.at(i)->name.Data(),v.at(i)->name.Data());
      TLegend* lg1 = 0;
      lg1 = new TLegend(0.45,0.7,0.97,0.95,"");
      lg1->SetNColumns(2);
      TH1D *histoverlay;
      if(debug) printf("set hists\n");
      for(auto& iter:stackorder ){
        if(iter == "data") continue;
        if(debug) {
          printf("plot_lib[%s][%s][%d]\n", iter.Data(), region.Data(), i);
        }
        TH1D *tmp = grabhist(iter,region,NPname,i);
        if(tmp) buffer.push_back((TH1D*)tmp->Clone());
        else continue;
        if(doROC && sensitivevariable == v.at(i)->name)
        {
          if(!ROC_bkg) ROC_bkg = (TH1D*) buffer.back()->Clone();
          else ROC_bkg->Add(buffer.back());
        }
        if(v.at(i)->rebin != 1) buffer.back()->Rebin(v.at(i)->rebin);
        hsk->Add(buffer.back());
        hmc.Add(buffer.back());
        lg1->AddEntry(buffer.back(),buffer.back()->GetTitle(),"F");
      }
      if(!hsk->GetMaximum()){
        printf("histSaver::plot_stack(): ERROR: stack has no entry for region %s, var %s, continue\n", region.Data(), v.at(i)->name.Data());
        continue;
      }
      double histmax = hmc.GetMaximum() + hmc.GetBinError(hmc.GetMaximumBin());

      TH1D * datahist = 0;
      if(debug) printf("set data\n");
      if (dataref) {
        datahistorig = grabhist("data",region,"NOMINAL",i);
        if(datahistorig) datahist = (TH1D*)datahistorig->Clone("dataClone");
        if(!datahist) {
          printf("histSaver::plot_stack(): WARNING: clone data histogram failed: region %s, variable %s\n", region.Data(), v.at(i)->name.Data());
          continue;
        } 
        if(v.at(i)->rebin != 1)
          datahist->Rebin(v.at(i)->rebin);
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

      cv.SaveAs(outdir + "/" + region + "/" + v[i]->name + ".pdf[");
      cv.cd();
      padhi->Draw();
      padhi->SetBottomMargin(0.017);
      padhi->SetRightMargin(0.08);
      padhi->SetLeftMargin(0.12);
      padhi->cd();
      
      hsk->Draw("hist");
      hsk->GetXaxis()->SetTitle(v.at(i)->unit == "" ? v.at(i)->title.Data() : (v.at(i)->title + " [" + v.at(i)->unit + "]").Data());
      hsk->GetXaxis()->SetLabelColor(kWhite);
      char str[30];
      sprintf(str,"Events / %4.2f %s",binwidth(i)*v.at(i)->rebin, v.at(i)->unit.Data());
      hsk->GetYaxis()->SetTitle(str);
      hsk->GetYaxis()->SetTitleOffset(1.6);
      hsk->GetYaxis()->SetLabelSize(hsk->GetYaxis()->GetLabelSize()*0.7);
      hsk->GetXaxis()->SetLabelSize(hsk->GetXaxis()->GetLabelSize()*0.7);
      hsk->GetYaxis()->SetTitleSize(hsk->GetYaxis()->GetTitleSize()*0.7);
      if(debug) printf("set blinding\n");

      for(Int_t j=1; j<v.at(i)->nbins+1; j++) {
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
      std::string regtitle = region.Data();
      findAndReplaceAll(regtitle,"reg","");
      findAndReplaceAll(regtitle,"vetobtagwp70","");
      findAndReplaceAll(regtitle,"_"," ");

      ATLASLabel(0.15,0.900,workflow.Data(),kBlack,lumi.Data(), analysis.Data(), regtitle.c_str());
      findAndReplaceAll(regtitle,"1l1tau1b3j_","TTH $\\tlhad$ ");
      findAndReplaceAll(regtitle,"1l1tau1b2j_","STH $\\tlhad$ ");
      findAndReplaceAll(regtitle,"1l1tau2b3j_","TTH $\\tlhad$ 2b ");
      findAndReplaceAll(regtitle,"1l1tau2b2j_","STH $\\tlhad$ 2b ");
      findAndReplaceAll(regtitle,"1l2tau1bnj_","$l\\thadhad$ ");
      findAndReplaceAll(regtitle,"1l2tau2bnj_","$l\\thadhad$ 2b ");
      findAndReplaceAll(regtitle,"2lSS1tau1bnj_","$2lSS\\thad$ ");
      findAndReplaceAll(regtitle,"2lSS1tau2bnj_","$2lSS\\thad$ 2b ");
//===============================blinded data===============================
      std::vector<TString> activeoverlay;
      for(auto overlaysample: overlaysamples){
        TH1D* histoverlaytmp = (TH1D*)grabhist(overlaysample,region,NPname,i);
        if(!histoverlaytmp){
          if(debug) printf("histSaver::plot_stack(): Warning: signal hist %s not found\n", overlaysample.Data());
          continue;
        }
        activeoverlay.push_back(overlaysample);
        if(blinding && dataref){
          for(Int_t j=1; j<v.at(i)->nbins+1; j++) {
            if(histoverlaytmp->GetBinContent(j)/sqrt(datahist->GetBinContent(j)) > blinding) {
              datahist->SetBinContent(j,0);
              datahist->SetBinError(j,0);
              hdataR.SetBinContent(j,0);
              hdataR.SetBinError(j,0);
            }
          }
          if(sensitivevariable == v.at(i)->name){
            for(int j = v.at(i)->nbins*3/4/v.at(i)->rebin ; j <= v.at(i)->nbins ; j++){
              datahist->SetBinContent(j,0);
              datahist->SetBinError(j,0);
              hdataR.SetBinContent(j,0);
              hdataR.SetBinError(j,0);
            }
          }
        }
      }
      if(dataref) {
        datahist->Draw("E same");
      }

//===============================lower pad===============================
      padlow->SetFillStyle(4000);
      padlow->SetGrid(1,1);
      padlow->SetTopMargin(0.03);
      padlow->SetBottomMargin(0.35);
      padlow->SetRightMargin(0.08);
      padlow->SetLeftMargin(0.12);
      padlow->cd();

      if(debug) printf("plot data ratio\n");
      if(dataref){
        hdataR.SetMarkerStyle(20);
        hdataR.SetMarkerSize(0.8);
      }
      hmcR.SetMaximum(1.5);
      hmcR.SetMinimum(0.5);
      hmcR.GetYaxis()->SetNdivisions(504,false);
      hmcR.GetYaxis()->SetTitle("Data/Bkg");
      hmcR.GetYaxis()->SetTitleOffset(hdataR.GetYaxis()->GetTitleOffset()*1.08);
      hmcR.GetYaxis()->CenterTitle();
      hmcR.GetXaxis()->SetTitle(v.at(i)->unit == "" ? v.at(i)->title.Data() : (v.at(i)->title + " [" + v.at(i)->unit + "]").Data());
      hmcR.GetXaxis()->SetTitleSize(hdataR.GetXaxis()->GetTitleSize()*0.7);
      hmcR.GetYaxis()->SetTitleSize(hdataR.GetYaxis()->GetTitleSize()*0.7);
      hmcR.SetFillColor(1);
      hmcR.SetLineColor(0);
      hmcR.SetMarkerStyle(1);
      hmcR.SetMarkerSize(0);
      hmcR.SetMarkerColor(1);
      hmcR.SetFillStyle(3004);
      hmcR.GetXaxis()->SetTitleOffset(3.4);
      hmcR.GetXaxis()->SetLabelSize(hdataR.GetXaxis()->GetLabelSize()*0.7); 
      hmcR.GetYaxis()->SetLabelSize(hdataR.GetYaxis()->GetLabelSize()*0.7); 
      if(debug) printf("plot mc ratio\n");
      hmcR.Draw("E2 same");
      if(debug) printf("plot data ratio\n");
      if(dataref){
        hdataR.Draw("E same");
      }
      TLine line;
      line.SetLineColor(2);
      line.DrawLine(hdataR.GetBinLowEdge(1), 1., hdataR.GetBinLowEdge(hdataR.GetNbinsX()+1), 1.);
      cv.cd();
      if(debug) printf("draw low pad\n");
      padlow->Draw();
      if(debug) printf("printing\n");

//===============================upper pad signal===============================

      padhi->cd();
      if(!activeoverlay.size()) {
        cv.SaveAs(outdir + "/" + region + "/" + v.at(i)->name + ".pdf");
      }

      if(sensitivevariable == v.at(i)->name) {
        if(dataref){
          yield_chart->set("data",regtitle,integral(datahistorig));
        }
        yield_chart->set("background",regtitle,integral(&hmc));
        printf("Region %s, Background yield: %f\n", region.Data(), hmc.Integral());
      }
      for(auto overlaysample: activeoverlay){
        
        TLegend *lgsig = (TLegend*) lg1->Clone();
        if(debug) { printf("overlay: %s\n", overlaysample.Data()); }
        histoverlay = grabhist(overlaysample,region,NPname,i);
        if(histoverlay) histoverlay = (TH1D*)histoverlay->Clone();
        if(doROC && sensitivevariable == v.at(i)->name) ROC_sig = (TH1D*) histoverlay->Clone();
        if(!histoverlay) continue;
        if(v.at(i)->rebin != 1) histoverlay->Rebin(v.at(i)->rebin);
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

        if(sensitivevariable == v.at(i)->name){
          double _significance = 0;
          for(Int_t j=1; j<v.at(i)->nbins+1; j++) {
            if(histoverlay->GetBinContent(j) && hmc.GetBinContent(j)) {
              if(hmc.GetBinContent(j) > 0 && histoverlay->GetBinContent(j) > 0)
                _significance += pow(significance(hmc.GetBinContent(j), histoverlay->GetBinContent(j)),2);
            }
          }
          if(doROC && sensitivevariable == v.at(i)->name){
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
            outputrocfile->cd();
            ROC->Write(overlaysample + "_ROC");
            ROC_sig->Write(overlaysample + "_ROC_sig");
            ROC_bkg->Write(overlaysample + "_ROC_bkg");
            deletepointer(ROC);
            deletepointer(ROC_sig);
            deletepointer(ROC_bkg);
          }
          std::string samptitle = histoverlay->GetTitle();
          findAndReplaceAll(samptitle," ","~");
          if(samptitle.find("#") != string::npos) samptitle = "$"+samptitle+"$";
          findAndReplaceAll(samptitle,"#","\\");
          findAndReplaceAll(samptitle,"%","\\%");
          findAndReplaceAll(samptitle,"rightarrow","to ");
          yield_chart->set(samptitle,regtitle,integral(histoverlay));
          sgnf_chart->set(samptitle,regtitle,sqrt(_significance));
          printf("signal %s yield: %4.2f, significance: %4.2f\n",overlaysample.Data(), histoverlay->Integral(), sqrt(_significance));
        }

        if(ratio > 0) histoverlay->Scale(ratio);
        if(overlaysample != "") histoverlay->Draw("hist same");
        lgsig->SetBorderSize(0);
        lgsig->Draw();
        padhi->Update();
        cv.SaveAs(outdir + "/" + region + "/" + v.at(i)->name + ".pdf");
        deletepointer(histoverlay);
        deletepointer(lgsig);
      }
      deletepointer(hsk);
      deletepointer(lg1);
      deletepointer(padlow );
      deletepointer(padhi  );
      deletepointer(datahist);
      for(auto &iter : buffer) deletepointer(iter);
      if(debug) printf("end region %s\n",region.Data());
      cv.SaveAs(outdir + "/" + region + "/" + v.at(i)->name + ".pdf]");
      cv.Clear();
    }
    if(debug) printf("end loop region\n");
  }
  if(yield_chart->rows.size()){
    yield_chart->caption = "The sample and data yield before the fit.";
    yield_chart->print((outputchartdir + "/" + "yield_chart").Data());
    sgnf_chart->caption = "The stat. only significance of the signal in each regions with the benchmark $\\mu$ value.";
    sgnf_chart->print((outputchartdir + "/" + "significance_chart").Data());
  }
  deletepointer(yield_chart);
  deletepointer(sgnf_chart);
  outputrocfile->Close();
  deletepointer(outputrocfile);
}
